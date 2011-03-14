#ifndef __FC_PAIRING_HEAP__
#define __FC_PAIRING_HEAP__

////////////////////////////////////////////////////////////////////////////////
// File    : FCPairHeap.h
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com
//           Ms.Moran Tzafrir  email: morantza@gmail.com
// Written : 16 February 2011, 27 October 2009
// 
// Copyright (C) 2011 Jonathan Eastep, 2009 Moran Tzafrir
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License 
// along with this program; if not, write to the Free Software Foundation
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
////////////////////////////////////////////////////////////////////////////////
// TODO:
//
////////////////////////////////////////////////////////////////////////////////

#include "cpp_framework.h"
#include "FCBase.h"
#include "PairingHeap.h"

//#define _FC_CAS_STATS

using namespace CCP;

template <class T>
class FCPairHeap : public FCBase<T> {
private:

        //constants -----------------------------------

        //inner classes -------------------------------

        // Private variables --------------------------

        //fields --------------------------------------
        AtomicInteger           _fc_lock;
        char                    _pad1[CACHE_LINE_SIZE];
        final int               _NUM_REP;
        final int               _REP_THRESHOLD;
        PairHeap<T>             _heap;

        //helper function -----------------------------
        inline_ void flat_combining() {
                final int repThresh = 0;
                for (int iTry=0;iTry<_NUM_REP; ++iTry) {
                        Memory::read_barrier();
                        int num_changes=0;
                        SlotInfo* curr_slot = FCBase<T>::_tail_slot.get();
                        while(null != curr_slot->_next) {
                                final FCIntPtr curr_value = curr_slot->_req_ans;
                                if(curr_value > FCBase<T>::_NULL_VALUE) {
                                        ++num_changes;
                                        _heap.insert((PtrNode<T>*) curr_value);
                                        curr_slot->_req_ans = FCBase<T>::_NULL_VALUE;
                                        curr_slot->_time_stamp = FCBase<T>::_NULL_VALUE;
                                } else if(FCBase<T>::_DEQ_VALUE == curr_value) {
                                        ++num_changes;
                                        curr_slot->_req_ans = -((FCIntPtr) _heap.deleteMin());
                                        curr_slot->_time_stamp = FCBase<T>::_NULL_VALUE;
                                }
                                curr_slot = curr_slot->_next;
                        }//while on slots
                        if(num_changes < repThresh)
                                break;
                        //Memory::write_barrier();
                }//for repetition
        }       
        
public:

        FCPairHeap()
        :       _NUM_REP(FCBase<T>::_NUM_THREADS),
                _REP_THRESHOLD((int)(Math::ceil(FCBase<T>::_NUM_THREADS/(1.7))))
        { }

        //enq ......................................................
        boolean add(final int iThread, PtrNode<T>* final inPtr) {
                final FCIntPtr inValue = (FCIntPtr) inPtr;
                CasInfo& my_cas_info = FCBase<T>::_cas_info_ary[iThread];

                SlotInfo* my_slot = FCBase<T>::_tls_slot_info.get();
                if(null == my_slot)
                        my_slot = FCBase<T>::get_new_slot();

                SlotInfo* volatile&     my_next = my_slot->_next;
                FCIntPtr volatile& my_re_ans = my_slot->_req_ans;
                my_re_ans = inValue;

                do {
                        //this is needed because the combiner may remove you
                        if (null == my_next)
                                FCBase<T>::enq_slot(my_slot);

                        boolean is_cas = false;
                        if(lock_fc(_fc_lock, is_cas)) {
#ifdef _FC_CAS_STATS
                                ++(my_cas_info._succ);
#endif
                                ++(my_cas_info._locks);
                                FCBase<T>::machine_start_fc(iThread);
                                flat_combining();
                                _fc_lock.set(0);
                                FCBase<T>::machine_end_fc(iThread);
#ifdef _FC_CAS_STATS
                                ++(my_cas_info._ops);
#endif
                                return true;
                        } else {
                                Memory::write_barrier();
#ifdef _FC_CAS_STATS
                                if(!is_cas)
                                        ++(my_cas_info._failed);
#endif
                                while(FCBase<T>::_NULL_VALUE != my_re_ans && 0 != _fc_lock.getNotSafe()) {
                                        FCBase<T>::thread_wait(iThread);
                                } 
                                Memory::read_barrier();
                                if(FCBase<T>::_NULL_VALUE == my_re_ans) {
#ifdef _FC_CAS_STATS
                                        ++(my_cas_info._ops);
#endif
                                        return true;
                                }
                        }
                } while(true);
        }

        //deq ......................................................
        PtrNode<T>* remove(final int iThread, PtrNode<T>* final inPtr) {
                final FCIntPtr inValue = (FCIntPtr) inPtr;
                CasInfo& my_cas_info = FCBase<T>::_cas_info_ary[iThread];

                SlotInfo* my_slot = FCBase<T>::_tls_slot_info.get();
                if(null == my_slot)
                        my_slot = FCBase<T>::get_new_slot();

                SlotInfo* volatile&     my_next = my_slot->_next;
                FCIntPtr volatile& my_re_ans = my_slot->_req_ans;
                my_re_ans = FCBase<T>::_DEQ_VALUE;

                do {
                        //this is needed because the combiner may remove you
                        if(null == my_next)
                                FCBase<T>::enq_slot(my_slot);

                        boolean is_cas = false;
                        if(lock_fc(_fc_lock, is_cas)) {
#ifdef _FC_CAS_STATS
                                ++(my_cas_info._succ);
#endif
                                ++(my_cas_info._locks);
                                FCBase<T>::machine_start_fc(iThread);
                                flat_combining();
                                _fc_lock.set(0);
                                FCBase<T>::machine_end_fc(iThread);
#ifdef _FC_CAS_STATS
                                ++(my_cas_info._ops);
#endif
                                return (PtrNode<T>*) -(my_re_ans);
                        } else {
                                Memory::write_barrier();
#ifdef _FC_CAS_STATS
                                if(!is_cas)
                                        ++(my_cas_info._failed);
#endif
                                while(FCBase<T>::_DEQ_VALUE == my_re_ans && 0 != _fc_lock.getNotSafe()) {
                                        FCBase<T>::thread_wait(iThread);
                                }
                                Memory::read_barrier();
                                if(FCBase<T>::_DEQ_VALUE != my_re_ans) {
#ifdef _FC_CAS_STATS
                                        ++(my_cas_info._ops);
#endif
                                        return (PtrNode<T>*) -(my_re_ans);
                                }
                        }
                } while(true);
        }

        //peek .....................................................
        PtrNode<T>* contain(final int iThread, PtrNode<T>* final inPtr) {
                final FCIntPtr inValue = (FCIntPtr) inPtr;
                return FCBase<T>::_NULL_VALUE;
        }

        //general .....................................................
        int size() {
                return 0;
        }

        final char* name() {
                return "FCPairHeap";
        }

};

#endif
