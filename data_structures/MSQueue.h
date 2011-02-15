#ifndef __MS_QUEUE__
#define __MS_QUEUE__

////////////////////////////////////////////////////////////////////////////////
// File    : MSQueue.h
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

#include "FCBase.h"
#include "cpp_framework.h"

using namespace CCP;

template <class T> class MSQueue  : public FCBase<T> {
private:

        struct Node {
                AtomicStampedReference<Node>    _next;
                FCIntPtr volatile               _value;

                Node(final FCIntPtr x) : _value(x) {}
        };

        AtomicStampedReference<Node>    _head;
        AtomicStampedReference<Node>    _tail;

public:
        MSQueue() {
                // Allocate a free node
                Node* const new_node = new Node(FCBase<T>::_NULL_VALUE);

                // Make it the only node in the linked list
                new_node->_next.set(null,0);
                _head.set(new_node,0);
                _tail.set(new_node,0);
        }

        boolean add(final int iThread, PtrNode<T>* final inPtr) { 
                final FCIntPtr inValue = (FCIntPtr) inPtr;

                CasInfo& my_cas_info = FCBase<T>::_cas_info_ary[iThread];

                // Allocate a new node from the free list
                // Set next pointer of node to NULL
                Node* new_node = new Node(inValue);                             
                new_node->_next.set(null, 0);          
                
                AtomicStampedReference<Node> tail;
                AtomicStampedReference<Node> next;

                // Keep trying until Enqueue is done
                do {                                
                        Memory::read_write_barrier();
                        // Read Tail.ptr and Tail.count together
                        // Read next ptr and count fields together
                        tail = _tail;                                                                   
                        next = tail.getReference()->_next;

                        // Are tail and next consistent?
                        if(tail == _tail) {
                                // Was Tail pointing to the last node?
                                if (null == next.getReference()) {
                                        // Try to link node at the end of the linked list
                                        if(tail.getReference()->_next.compareAndSet(next.getReference(), new_node, next.getStamp(), next.getStamp()+1)) {
                                                ++(my_cas_info._succ);
                                                // Enqueue is done.  Exit loop
                                                break;  
                                        } else {
                                                ++(my_cas_info._failed);
                                        }
                                } else {                
                                        // Tail was not pointing to the last node Try to swing Tail to the next node
                                        if(_tail.compareAndSet(tail.getReference(), next.getReference(), tail.getStamp(), next.getStamp()+1))
                                                ++(my_cas_info._succ);
                                        else
                                                ++(my_cas_info._failed);
                                }
                        }
                } while(true);

                // Enqueue is done.  Try to swing Tail to the inserted node
                if (_tail.compareAndSet(tail, new_node, tail.getStamp(), tail.getStamp()+1))
                        ++(my_cas_info._succ);
                else 
                        ++(my_cas_info._failed);

                ++(my_cas_info._ops);
                return true;
        }

        PtrNode<T>* remove(final int iThread, PtrNode<T>* final inPtr) { 
                final FCIntPtr inValue = (FCIntPtr) inPtr;

                CasInfo& my_cas_info = FCBase<T>::_cas_info_ary[iThread];

                AtomicStampedReference<Node> head;
                AtomicStampedReference<Node> tail;
                AtomicStampedReference<Node> next;

                // Keep trying until Dequeue is done
                do {                              
                        Memory::read_write_barrier();
                        // Read Head and Tail and Head.ptr-next
                        head = _head;
                        tail = _tail;
                        next = head.getReference()->_next;

                        // Are head, tail, and next consistent?
                        if(head == _head) { 
                                // Is queue empty or Tail falling behind?
                                if(head.getReference() == tail.getReference()) {
                                        // Is queue empty?
                                        if (null == next.getReference()) {                                       
                                                ++(my_cas_info._ops);
                                                // Queue is empty, couldn't dequeue
                                                return FCBase<T>::_NULL_VALUE; 
                                        }
                                        // Tail is falling behind.  Try to advance it
                                        if(_tail.compareAndSet(tail.getReference(), next.getReference(), tail.getStamp(), tail.getStamp()+1))
                                                ++(my_cas_info._succ);
                                        else
                                                ++(my_cas_info._failed);

                                } else {                     
                                        // No need to deal with Tail
                                        // Read value before CAS Otherwise, another dequeue might free the next node
                                        final FCIntPtr rtrn_value = next.getReference()->_value;

                                        // Try to swing Head to the next node
                                        if (_head.compareAndSet(head.getReference(), next.getReference(), head.getStamp(), head.getStamp()+1)) {
                                                ++(my_cas_info._succ);
                                                // It is safe now to free the old node
                                                //free(head.ptr)                     
                                                ++(my_cas_info._ops);
                                                // Queue was not empty, dequeue succeeded  
                                                return (PtrNode<T>*) rtrn_value;     
                                        } else {
                                                ++(my_cas_info._failed);
                                        }
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
                return "MSQueue";
        }

};

#endif

