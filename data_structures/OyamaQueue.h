#ifndef __OYAMA_QUEUE__
#define __OYAMA_QUEUE__

////////////////////////////////////////////////////////////////////////////////
// File    : OyamaQueue.h
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com
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
// MERCHANTABILITY99 or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
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
#include "QueueForLogSync.h"

using namespace CCP;

template <class T> class OyamaQueue : public FCBase<T> {
private:

        struct Node {
                const FCIntPtr      _value;
                Node* volatile      _next;

                Node(final FCIntPtr x) : _value(x) {
                        _next = null;
                }
        };

        struct Req {
                VolatileType<FCIntPtr> _req_ans;
                int                         _d1,d2,d3,d4;

                Req(final FCIntPtr in_req) : _req_ans(in_req) { }
        };

        AtomicInteger           _log_lock;
        QueueForLogSync<Req>    _log;
        Node* volatile          _head;
        Node* volatile          _tail;

        void execute_log(CasInfo& my_cas_info) {
                for (int i=0; i<FCBase<T>::_NUM_THREADS; ++i) {
                        Req* final curr_req = _log.deq(my_cas_info);
                        if(null == curr_req)
                                break;

                        final FCIntPtr req_ans = curr_req->_req_ans;
                        if(req_ans > FCBase<T>::_NULL_VALUE) {
                                Node* const new_node = new Node(req_ans);       // Allocate a new node from the free list
                                new_node->_next = null;                         // Set next pointer of node to NULL

                                _tail->_next = new_node;                        // Link node at the end of the linked list
                                _tail = new_node;                               // Swing Tail to node

                                 curr_req->_req_ans = FCBase<T>::_NULL_VALUE;
                                 //delete curr_req;

                        } else if(FCBase<T>::_DEQ_VALUE == req_ans) {

                                Node* const old_head = _head;                   // Read Head
                                Node* const new_head = old_head->_next;         // Read next pointer

                                if(null == new_head) {                                          // Is queue empty?
                                        curr_req->_req_ans = FCBase<T>::_NULL_VALUE;
                                } else {
                                        final FCIntPtr curr_value = new_head->_value;      // Queue not empty.  Read value before release
                                        _head = new_head;                                       // Swing Head to next node
                                        curr_req->_req_ans = -(curr_value);
                                }
                                //delete curr_req;
                        }
                }
        }

public:
        OyamaQueue() {
                Node* const new_node = new Node(FCBase<T>::_NULL_VALUE);
                new_node->_next = null;
                _head = _tail = new_node;
        }

        boolean add(final int iThread, PtrNode<T>* final inPtr) { 
                final FCIntPtr inValue = (FCIntPtr) inPtr;

                Req* final my_req( new Req(inValue) );
                CasInfo& my_cas_info = FCBase<T>::_cas_info_ary[iThread];
                _log.enq(my_req, my_cas_info);

                do {
                        boolean is_cas = true;
                        if(lock_fc(_log_lock, is_cas)) {
                                ++(my_cas_info._succ);
                                execute_log(my_cas_info);
                                _log_lock.set(0);
                                ++(my_cas_info._ops);
                                return true;
                        } else {
                                if(!is_cas)
                                        ++(my_cas_info._failed);
                                Memory::write_barrier();
                                while(FCBase<T>::_NULL_VALUE != my_req->_req_ans && 0 != _log_lock.getNotSafe()) {
                                        FCBase<T>::thread_wait(iThread);
                                } 
                                Memory::read_barrier();
                                if(FCBase<T>::_NULL_VALUE !=  my_req->_req_ans) {
                                        ++(my_cas_info._ops);
                                        return true;
                                }
                        }
                } while(true);
        }

        PtrNode<T>* remove(final int iThread, PtrNode<T>* final inPtr) { 
                final FCIntPtr inValue = (FCIntPtr) inPtr;

                Req* final my_req( new Req(FCBase<T>::_DEQ_VALUE) );
                CasInfo& my_cas_info = FCBase<T>::_cas_info_ary[iThread];
                _log.enq(my_req, my_cas_info);

                do {
                        boolean is_cas = true;
                        if(lock_fc(_log_lock, is_cas)) {
                                ++(my_cas_info._succ);
                                execute_log(my_cas_info);
                                _log_lock.set(0);
                                ++(my_cas_info._ops);
                                return (PtrNode<T>*) -(my_req->_req_ans);
                        } else {
                                if(!is_cas)
                                        ++(my_cas_info._failed);
                                Memory::write_barrier();
                                while(FCBase<T>::_DEQ_VALUE == (my_req->_req_ans) && 0 != _log_lock.getNotSafe()) {
                                        FCBase<T>::thread_wait(iThread);
                                } 
                                Memory::read_barrier();
                                if(FCBase<T>::_DEQ_VALUE !=  my_req->_req_ans) {
                                        ++(my_cas_info._ops);
                                        return (PtrNode<T>*) -(my_req->_req_ans);
                                }
                        }
                } while(true);
        }

        PtrNode<T>* contain(final int iThread, PtrNode<T>* final inPtr) { 
                final FCIntPtr inValue = (FCIntPtr) inPtr;
                return FCBase<T>::_NULL_VALUE;
        }

        int size() {
                return 0;
        }

        const char* name() {
                return "OyamaQueue";
        }

};


#endif

