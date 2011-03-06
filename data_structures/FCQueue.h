#ifndef __FC_QUEUE__
#define __FC_QUEUE__

////////////////////////////////////////////////////////////////////////////////
// File    : FCQueue.h
// Author  : Jonathan Eastep   email: jonathan.eastep@gmail.com
//           Ms.Moran Tzafrir  email: morantza@gmail.com
// Written : 16 February 2011, 27 October 2009
//
// Copyright (C) 2011 Jonathan Eastep, 2009 Moran Tzafrir.
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

//#define FCCASSTATS


using namespace CCP;

template <class T>
class FCQueue : public FCBase<T> {
private:

        //constants -----------------------------------

        //inner classes -------------------------------
        struct Node {
                Node* volatile          _next;
                FCIntPtr volatile       _values[256];

                static Node* get_new(final int in_num_values) {
                        final size_t new_size = (sizeof(Node) + (in_num_values + 2 - 256) * sizeof(FCIntPtr));

                        Node* final new_node = (Node*) malloc(new_size);
                        new_node->_next = null;
                        return new_node;
                }
        };

        //fields --------------------------------------
        AtomicInteger   _fc_lock;
        char            _pad1[CACHE_LINE_SIZE];
        final int       _NUM_REP;
        final int       _REP_THRESHOLD;
        Node* volatile  _head;
        Node* volatile  _tail;
        int volatile    _NODE_SIZE;
        Node* volatile  _new_node;

        //helper function -----------------------------
        inline_ void flat_combining() {
                
                // prepare for enq
                FCIntPtr volatile* enq_value_ary;
                if(null == _new_node) 
                        _new_node = Node::get_new(_NODE_SIZE);

                enq_value_ary = _new_node->_values;
                *enq_value_ary = 1;
                ++enq_value_ary;

                // prepare for deq
                FCIntPtr volatile * deq_value_ary = _tail->_values;
                deq_value_ary += deq_value_ary[0];

                int num_added = 0;
                //final int repThresh = _REP_THRESHOLD;
                final int repThresh = 0;

                for (int iTry=0;iTry<_NUM_REP; ++iTry) {
                        //Memory::read_barrier();

                        int num_changes=0;
                        SlotInfo* curr_slot = FCBase<T>::_tail_slot.get();
                        while(null != curr_slot->_next) {
                                final FCIntPtr curr_value = curr_slot->_req_ans;
                                if(curr_value > FCBase<T>::_NULL_VALUE) {
                                        ++num_changes;
                                        *enq_value_ary = curr_value;
                                        ++enq_value_ary;
                                        curr_slot->_req_ans = FCBase<T>::_NULL_VALUE;
                                        curr_slot->_time_stamp = FCBase<T>::_NULL_VALUE;

                                        ++num_added;
                                        if(num_added >= _NODE_SIZE) {
                                                Node* final new_node2 = Node::get_new(_NODE_SIZE+4);
                                                memcpy((void*)(new_node2->_values), (void*)(_new_node->_values), (_NODE_SIZE+2)*sizeof(FCIntPtr) );
                                                free(_new_node);
                                                _new_node = new_node2; 
                                                enq_value_ary = _new_node->_values;
                                                *enq_value_ary = 1;
                                                ++enq_value_ary;
                                                enq_value_ary += _NODE_SIZE;
                                                _NODE_SIZE += 4;
                                        }
                                } else if(FCBase<T>::_DEQ_VALUE == curr_value) {
                                        ++num_changes;
                                        final FCIntPtr curr_deq = *deq_value_ary;
                                        if(0 != curr_deq) {
                                                curr_slot->_req_ans = -curr_deq;
                                                curr_slot->_time_stamp = FCBase<T>::_NULL_VALUE;
                                                ++deq_value_ary;
                                        } else if(null != _tail->_next) {
                                                Node* tmp = _tail;
                                                _tail = _tail->_next;
                                                free(tmp);
                                                deq_value_ary = _tail->_values;
                                                deq_value_ary += deq_value_ary[0];
                                                continue;
                                        } else {
                                                curr_slot->_req_ans = FCBase<T>::_NULL_VALUE;
                                                curr_slot->_time_stamp = FCBase<T>::_NULL_VALUE;
                                        } 
                                }
                                curr_slot = curr_slot->_next;
                        }//while on slots

                        //Memory::write_barrier();
                        if(num_changes < repThresh)
                                break;
                }//for repetition

                if(0 == *deq_value_ary && null != _tail->_next) {
                        Node* tmp = _tail;
                        _tail = _tail->_next;
                        free(tmp);
                } else {
                        _tail->_values[0] = (deq_value_ary -  _tail->_values);
                }

                if(enq_value_ary != (_new_node->_values + 1)) {
                        *enq_value_ary = 0;
                        _head->_next = _new_node;
                        _head = _new_node;
                        _new_node  = null;
                } 
        }

public:
        //public operations ---------------------------
        FCQueue() 
	:       _NUM_REP(FCBase<T>::_NUM_THREADS),
                _REP_THRESHOLD((int)(Math::ceil(FCBase<T>::_NUM_THREADS/(1.7))))
        {
                _head = Node::get_new(FCBase<T>::_NUM_THREADS);
                _tail = _head;
                _head->_values[0] = 1;
                _head->_values[1] = 0;

                FCBase<T>::_timestamp = 0;
                _NODE_SIZE = 4;
                _new_node = null;
        }

        virtual ~FCQueue() { }

        //enq ......................................................
        boolean add(final int iThread, PtrNode<T>* final inPtr) {
                final FCIntPtr inValue = (FCIntPtr) inPtr;
                CasInfo& my_cas_info = FCBase<T>::_cas_info_ary[iThread];

                //SlotInfo* my_slot = _tls_slot_info;
                SlotInfo* my_slot = FCBase<T>::_tls_slot_info.get();
                if(null == my_slot)
                        my_slot = FCBase<T>::get_new_slot();

                SlotInfo* volatile&  my_next   = my_slot->_next;
                FCIntPtr volatile*   my_re_ans = &my_slot->_req_ans;
                *my_re_ans = inValue;

                do {
                        //this is needed because the combiner may remove you
                        if (null == my_next)
                                FCBase<T>::enq_slot(my_slot);

                        boolean is_cas = true;
                        if(lock_fc(_fc_lock, is_cas)) {
#ifdef FCCASSTATS
                                ++(my_cas_info._succ);
#endif
                                ++(my_cas_info._locks);
                                FCBase<T>::machine_start_fc(iThread);
                                flat_combining();
                                _fc_lock.set(0);
                                FCBase<T>::machine_end_fc(iThread);
#ifdef FCCASSTATS
                                ++(my_cas_info._ops);
#endif
                                return true;
                        } else {
                                //Memory::write_barrier();
#ifdef FCCASSTATS
                                if(!is_cas)
                                        ++(my_cas_info._failed);
#endif
                                while(FCBase<T>::_NULL_VALUE != *my_re_ans && 0 != _fc_lock.getNotSafe()) {
                                        FCBase<T>::thread_wait(iThread);
                                } 
                                //Memory::read_barrier();
                                if(FCBase<T>::_NULL_VALUE == *my_re_ans) {
#ifdef FCCASSTATS
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

                //SlotInfo* my_slot = _tls_slot_info;
                SlotInfo* my_slot = FCBase<T>::_tls_slot_info.get();
                if(null == my_slot)
                        my_slot = FCBase<T>::get_new_slot();

                SlotInfo* volatile&     my_next = my_slot->_next;
                FCIntPtr volatile* my_re_ans = &my_slot->_req_ans;
                *my_re_ans = FCBase<T>::_DEQ_VALUE;

                do {
                        //this is needed because the combiner may remove you
                        if(null == my_next)
                                FCBase<T>::enq_slot(my_slot);

                        boolean is_cas = true;
                        if(lock_fc(_fc_lock, is_cas)) {
#ifdef FCCASSTATS
                                ++(my_cas_info._succ);
#endif
                                ++(my_cas_info._locks);
                                FCBase<T>::machine_start_fc(iThread);
                                flat_combining();
                                _fc_lock.set(0);
                                FCBase<T>::machine_end_fc(iThread);
#ifdef FCCASSTATS
                                ++(my_cas_info._ops);
#endif
                                return (PtrNode<T>*) -(*my_re_ans);
                        } else {
                                //Memory::write_barrier();
#ifdef FCCASSTATS
                                if(!is_cas)
                                        ++(my_cas_info._failed);
#endif
                                while(FCBase<T>::_DEQ_VALUE == *my_re_ans && 0 != _fc_lock.getNotSafe()) {
                                        FCBase<T>::thread_wait(iThread);
                                }
                                //Memory::read_barrier();
                                if(FCBase<T>::_DEQ_VALUE != *my_re_ans) {
#ifdef FCCASSTATS
                                        ++(my_cas_info._ops);
#endif
                                        return (PtrNode<T>*) -(*my_re_ans);
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
                return "FCQueue";
        }

};


#endif
