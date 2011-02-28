#ifndef __SMART_PAIRING_HEAP__
#define __SMART_PAIRING_HEAP__

////////////////////////////////////////////////////////////////////////////////
// File    : SmartPairHeap.h
// Author  : Jonathan Eastep  email: jonathan.eastep@gmail.com
// Written : 16 February 2011
//
// Copyright (C) 2011 Jonathan Eastep
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

#include "PairingHeap.h"
#include "FCBase.h"
#include "LearningEngine.h"
#include "SmartLockLite.h"
#include "Heartbeat.h"

using namespace CCP;

template <class T>
class SmartPairHeap : public FCBase<T> {
private:

        //constants -----------------------------------

        //inner classes -------------------------------

        // Private variables --------------------------

        //fields --------------------------------------
        SmartLockLite<FCIntPtr>*  _fc_lock;
        char                      _pad1[CACHE_LINE_SIZE];
        final int                 _NUM_REP;
        final int                 _REP_THRESHOLD;
        char                      _pad2[CACHE_LINE_SIZE];
        PairHeap<T>               _heap;
        char                      _pad3[CACHE_LINE_SIZE];
        Hb*                       _hbmon;
        LearningEngine*           _learner;
        int                       _sc_tune_id;

        //helper function -----------------------------
        inline_ void flat_combining(final int iThread) {

                int maxPasses;
                if ( 0 == FCBase<T>::_enable_scancount_tuning )
                        maxPasses = FCBase<T>::_num_passes;
                else
		        maxPasses = 1 + 2*_learner->getdiscval(_sc_tune_id, iThread);

                int total_changes = 0;
                for (int iTry=0;iTry<maxPasses; ++iTry) {
		        int num_changes = 0;
                        //Memory::read_barrier();

                        SlotInfo* curr_slot = FCBase<T>::_tail_slot.get();
                        while(null != curr_slot->_next) {
                                final FCIntPtr curr_value = curr_slot->_req_ans;
                                if(curr_value > FCBase<T>::_NULL_VALUE) {
                                        if ( 0 == _gIsDedicatedMode )
                                                ++num_changes;
                                        _heap.insert((PtrNode<T>*) curr_value);
                                        curr_slot->_req_ans = FCBase<T>::_NULL_VALUE;
                                        curr_slot->_time_stamp = FCBase<T>::_NULL_VALUE;
                                } else if(FCBase<T>::_DEQ_VALUE == curr_value) {
                                        final FCIntPtr rv = -((FCIntPtr) _heap.deleteMin());
                                        if ( (rv != FCBase<T>::_NULL_VALUE) || (0 == _gIsDedicatedMode) )
                                                ++num_changes;
                                        curr_slot->_req_ans = rv;
                                        curr_slot->_time_stamp = FCBase<T>::_NULL_VALUE;
                                }
                                curr_slot = curr_slot->_next;
                        }//while on slots

                        total_changes += num_changes;

                }//for repetition

                if ( total_changes && (FCBase<T>::_enable_scancount_tuning || FCBase<T>::_enable_lock_scheduling) )
                        _hbmon->heartbeat(total_changes);

        }       
        
public:

        SmartPairHeap(Hb* hbmon, LearningEngine* learner)
        :       _NUM_REP(FCBase<T>::_NUM_THREADS),
                _REP_THRESHOLD((int)(Math::ceil(FCBase<T>::_NUM_THREADS/(1.7)))),
	        _hbmon(hbmon),
	        _learner(learner)
        {
                _sc_tune_id = 0;
                if ( 0 != FCBase<T>::_enable_scancount_tuning )
                        _sc_tune_id = _learner->register_sc_tune_id();

                _fc_lock = new SmartLockLite<FCIntPtr>(FCBase<T>::_NUM_THREADS, _learner);
                Memory::read_write_barrier();
        }

        virtual ~SmartPairHeap() 
        {
                delete _fc_lock;
        }

        //enq ......................................................
        boolean add(final int iThread, PtrNode<T>* final inPtr) {
                final FCIntPtr inValue = (FCIntPtr) inPtr;
 
                SlotInfo* my_slot = FCBase<T>::_tls_slot_info.get();
                if(null == my_slot)
                        my_slot = FCBase<T>::get_new_slot();

                SlotInfo* volatile&     my_next = my_slot->_next;
                FCIntPtr volatile* my_re_ans = &my_slot->_req_ans;
                *my_re_ans = inValue;

                //this is needed because the combiner may remove you
                if (null == my_next)
                   FCBase<T>::enq_slot(my_slot);

                //Memory::write_barrier();
                boolean is_cas = _fc_lock->lock(my_re_ans, inValue, iThread);
                // when we get here, we either aborted or succeeded
                // abort happens when we got our answer
                if ( is_cas )
                {
                   // got the lock so we should do flat combining
                   CasInfo& my_cas_info = FCBase<T>::_cas_info_ary[iThread];
                   ++(my_cas_info._locks);
                   flat_combining(iThread);
                   _fc_lock->unlock(iThread);
                }
                return true;
        }

        //deq ......................................................
        PtrNode<T>* remove(final int iThread, PtrNode<T>* final inPtr) {
                final FCIntPtr inValue = (FCIntPtr) inPtr;

                SlotInfo* my_slot = FCBase<T>::_tls_slot_info.get();
                if(null == my_slot)
                        my_slot = FCBase<T>::get_new_slot();

                SlotInfo* volatile&     my_next = my_slot->_next;
                FCIntPtr volatile* my_re_ans = &my_slot->_req_ans;
                *my_re_ans = FCBase<T>::_DEQ_VALUE;

                //this is needed because the combiner may remove you
                if(null == my_next)
                   FCBase<T>::enq_slot(my_slot);

                //Memory::write_barrier();
                boolean is_cas = _fc_lock->lock(my_re_ans, FCBase<T>::_DEQ_VALUE, iThread);
                if( is_cas ) 
                {
                   CasInfo& my_cas_info = FCBase<T>::_cas_info_ary[iThread];
                   ++(my_cas_info._locks);
                   flat_combining(iThread);
                   _fc_lock->unlock(iThread);
                }
                return (PtrNode<T>*) -(*my_re_ans);
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
                return "SmartPairHeap";
        }

        void cas_reset(final int iThread) {
                _fc_lock->resetcasops(iThread);
                FCBase<T>::_cas_info_ary[iThread].reset();
        }

        void print_custom() {
                int failed = 0;
                int succ = 0;
                int ops = 0;
                int locks = 0;

                for (int i=0; i<FCBase<T>::_NUM_THREADS; ++i) {
                        failed += FCBase<T>::_cas_info_ary[i]._failed;
                        succ += FCBase<T>::_cas_info_ary[i]._succ;
                        ops += FCBase<T>::_cas_info_ary[i]._ops;
                        locks += FCBase<T>::_cas_info_ary[i]._locks;
                }
                int tmp1 = _fc_lock->getcasops();
                int tmp2 = _fc_lock->getcasfails();
                succ += tmp1 - tmp2;
                failed += tmp2;
                printf(" 0 0 0 0 0 0 ( %d, %d, %d, %d, %d )", ops, locks, succ, failed, failed+succ);
        }

};

#endif
