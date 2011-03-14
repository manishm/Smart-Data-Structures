#ifndef __FC_SKIP_LIST__
#define __FC_SKIP_LIST__

////////////////////////////////////////////////////////////////////////////////
// File    : FCSkipList.h
// Author  : Jonathan Eastep   email: jonathan.eastep@gmail.com
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

#include "FCBase.h"
#include "cpp_framework.h"

//#define _FC_CAS_STATS

using namespace CCP;

template <class T>
class FCSkipList : public FCBase<T> { 
protected://consts
        static final int _MAX_LEVEL     = 20;

protected://types

        class Node {
        public:
                FCIntPtr        _key;
                PtrNode<T>*     _element;
                int             _top_level;
                int             _counter;
                Node*           _next[_MAX_LEVEL+1];

        public:

                static Node* getNewNode(PtrNode<T>* final theElement) {
                        Node* final new_node  = (Node*) malloc(sizeof(Node));
                        new_node->_element    = theElement;
                        new_node->_key        = theElement->getkey();
                        new_node->_top_level  = _MAX_LEVEL;
                        new_node->_counter    = 1;
                        return new_node;
                }

                static Node* getNewNode(PtrNode<T>* final theElement, final int height) {
                        Node* final new_node  = (Node*) malloc(sizeof(Node) + height + 1 - _MAX_LEVEL);
                        new_node->_element    = theElement;
                        new_node->_key        = theElement->getkey();
                        new_node->_top_level  = height;
                        new_node->_counter    = 1;
                        return new_node;
                }
        };

protected://SkipList fields
        VolatileType<_u64>      _random_seed;
        Node*   final           _head;
        Node* final             _tail;

protected://Flat Combining fields

        AtomicInteger           _fc_lock;
        char                    _pad1[CACHE_LINE_SIZE];
        final int               _NUM_REP;
        final int               _REP_THRESHOLD;
        Node*                   preds[_MAX_LEVEL + 1];
        Node*                   succs[_MAX_LEVEL + 1];
        SlotInfo*               _saved_remove_node[1024];
	Node*                   _saved_node_ptr[1024];

protected://methods
        inline_ int randomLevel() {
                int x = (int)(_random_seed.get()  & U64(0xFFFFFF));
                x ^= x << 13;
                x ^= x >> 17;
                _random_seed.set( x ^= x << 5 );
                if ((x & 0x80000001) != 0) {// test highest and lowest bits
                        //printf("0\n");  fflush(stdout);
                        return 1;
                }
                int level = 2;
                while (((x >>= 1) & 1) != 0) 
                        ++level;
                //printf("%d\n",level);  fflush(stdout);
                if(level > (_MAX_LEVEL-1))
                        return (_MAX_LEVEL-1);
                else
                        return level;
        }

        inline_ Node* find(PtrNode<T>* final inPtr) {
                final FCIntPtr key = inPtr->getkey();
                Node* pPred;
                Node* pCurr;
                pPred = _head;
                Node* found_node = null;

                for (int iLevel = _MAX_LEVEL-1; iLevel >= 0; --iLevel) {

                        pCurr = pPred->_next[iLevel];

                        while (key > pCurr->_key) {
                                pPred = pCurr; 
                                pCurr = pPred->_next[iLevel];
                        }

                        if (null == found_node && key == pCurr->_key) {
                                found_node = pCurr;
                        }

                        preds[iLevel] = pPred;
                        succs[iLevel] = pCurr;
                }

                return found_node;
        }


        inline_ void flat_combining() {

                int num_removed = 0;
                final int top_level = randomLevel();
                final int repThresh = 0;
                for (int iTry=0;iTry<_NUM_REP; ++iTry) {

                        int num_changes=0;
                        SlotInfo* curr_slot = FCBase<T>::_tail_slot.get();
                        Memory::read_barrier();

                        while(null != curr_slot->_next) {

                                if ( curr_slot->_deq_pending ) {
                                        curr_slot = curr_slot->_next;
                                        continue;
                                }

                                final FCIntPtr inValue = curr_slot->_req_ans;
                                if(inValue > FCBase<T>::_NULL_VALUE) {
                                        ++num_changes;
                                        curr_slot->_req_ans = FCBase<T>::_NULL_VALUE;
                                        curr_slot->_time_stamp = FCBase<T>::_NULL_VALUE;

                                        //ADD ......................................................
                                        Node* node_found = find((PtrNode<T>*) inValue);
                                        if (null != node_found) {
                                                ++(node_found->_counter);
                                                curr_slot = curr_slot->_next;
                                                continue;
                                        }

                                        // first link succs ........................................
                                        // then link next fields of preds ..........................
                                        Node* new_node = Node::getNewNode((PtrNode<T>*) inValue, top_level);
                                        Node** new_node_next = new_node->_next;
                                        Node** curr_succ = succs;
                                        Node** curr_preds = preds;

                                        for (int level = 0; level < top_level; ++level) {
                                                *new_node_next = *curr_succ;
                                                (*curr_preds)->_next[level] = new_node;
                                                ++new_node_next;
                                                ++curr_succ;
                                                ++curr_preds;
                                        }

                                        //..........................................................
                                        curr_slot = curr_slot->_next;
                                        continue;

                                } else if(FCBase<T>::_DEQ_VALUE == inValue) {
                                        ++num_changes;
                                        curr_slot->_deq_pending = true;
                                        //REMOVE ...................................................
                                        _saved_remove_node[num_removed] = curr_slot;
                                        ++num_removed;
					assert(num_removed < 1024);
                                        curr_slot = curr_slot->_next;
                                        continue;
                                } else {
                                        curr_slot = curr_slot->_next;
                                        continue;
                                }

                        } //while on slots

                        if(num_changes < repThresh)
                                break;
                }

                //..................................................................
                Node* remove_node = (_head->_next[0]);
                int max_level = -1;
                int iSaved = 0;
                for (int iRemove = 0; iRemove<num_removed; ++iRemove) {

                        if ( _tail != remove_node ) {
                                SlotInfo* dequeuer = _saved_remove_node[iRemove];
                                dequeuer->_req_ans = -((FCIntPtr) remove_node->_element);
                                dequeuer->_time_stamp = FCBase<T>::_NULL_VALUE;
                                dequeuer->_deq_pending = false;

                                --(remove_node->_counter);
                                if(0 == remove_node->_counter) {
                                        if(remove_node->_top_level > max_level) {
                                                max_level = remove_node->_top_level;
                                        }
                                        _saved_node_ptr[iSaved++] = remove_node;
                                        remove_node = remove_node->_next[0];
                                } 
                        }
                        else
                        {
                                SlotInfo* dequeuer = _saved_remove_node[iRemove];
                                dequeuer->_req_ans = FCBase<T>::_NULL_VALUE;
                                dequeuer->_time_stamp = FCBase<T>::_NULL_VALUE;
                                dequeuer->_deq_pending = false;
                        }
                }

                if(-1 != max_level) {
                        Node* pred = _head;
                        Node* curr;

                        for (int iLevel = (max_level-1); iLevel >= 0; ) {
                                curr = pred->_next[iLevel];

                                if(0 != curr->_counter) {
                                        _head->_next[iLevel] = curr;
                                        --iLevel;
                                } else {
                                        pred = curr; 
                                        curr = pred->_next[iLevel];
                                }
                        }
                }

                for(int i = 0; i < iSaved; i++)
                        free(_saved_node_ptr[i]);

        }

public://methods

        FCSkipList()
        : _head( Node::getNewNode(new PtrNode<T>(FCBase<T>::_MIN_INT, null)) ),
          _tail( Node::getNewNode(new PtrNode<T>(FCBase<T>::_MAX_INT, null)) ),
          _NUM_REP( Math::Min(2, FCBase<T>::_NUM_THREADS)),
          _REP_THRESHOLD((int)(Math::ceil(FCBase<T>::_NUM_THREADS/(1.7))))
        {
                //initialize head to point to tail .....................................
                for (int iLevel = 0; iLevel < _head->_top_level; ++iLevel)
                        _head->_next[iLevel] = _tail;
        }

        ~FCSkipList() { }

public://methods

        //deq ......................................................
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

        //peek ......................................................................
        PtrNode<T>* contain(final int iThread, PtrNode<T>* final inPtr) {
                final FCIntPtr inValue = (FCIntPtr) inPtr;
                return FCBase<T>::_NULL_VALUE;
        }

public://methods

        int size() {
                return 0;
        }

        final char* name() {
                return "FCSkipList";
        }

};

#endif
