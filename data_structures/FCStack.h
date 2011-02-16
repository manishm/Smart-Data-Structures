#ifndef __FC_STACK__
#define __FC_STACK__

////////////////////////////////////////////////////////////////////////////////
// File    : FCStack.h
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

template <class T>
class FCStack : public FCBase<T> {
private:

        //fields -------------------------------------- 
        FCIntPtr volatile        _req_ary[FCBase<T>::_MAX_THREADS];

        AtomicInteger            _fc_lock;

        int volatile             _top_indx;
        int volatile             _stack_ary_size;
        FCIntPtr volatile*       _stack_ary;            

        //helper function -----------------------------
        inline_ void flat_combining() {
                //
                FCIntPtr volatile* final stack_ary_end = _stack_ary + _stack_ary_size;             
                register FCIntPtr volatile* stack_ary = &(_stack_ary[_top_indx]);
                register FCIntPtr volatile* req_ary;

                for (register int iTry=0; iTry<3; ++iTry) {
                        req_ary = _req_ary;
                        for(register int iReq=0; iReq < FCBase<T>::_NUM_THREADS;) {
                                register final FCIntPtr curr_value = *req_ary;
                                if(curr_value > FCBase<T>::_NULL_VALUE) {
                                        if(stack_ary_end != stack_ary) {
                                                *stack_ary = curr_value;
                                                ++stack_ary;
                                                *req_ary = FCBase<T>::_NULL_VALUE;
                                        }
                                } else if(FCBase<T>::_DEQ_VALUE == curr_value) {
                                        if(stack_ary  != _stack_ary) {
                                                --stack_ary;
                                                *req_ary = -(*stack_ary);
                                        } else {
                                                *req_ary = FCBase<T>::_NULL_VALUE;
                                        }
                                } 
                                ++req_ary;
                                ++iReq;
                        }
                }

                //
                _top_indx = (stack_ary - _stack_ary);
        }

public:
        //public operations ---------------------------
        FCStack() 
        {
                for (int iReq=0; iReq<FCBase<T>::_MAX_THREADS; ++iReq) {
                        _req_ary[iReq] = 0;
                }

                _stack_ary_size = 1024;
                _stack_ary = new FCIntPtr[_stack_ary_size];
                _top_indx = 0;
        }

        ~FCStack()
        {
                delete[] _stack_ary;
        }

        //push ......................................................
        boolean add(final int iThread, PtrNode<T>* final inPtr) {
                final FCIntPtr inValue = (FCIntPtr) inPtr;

                FCIntPtr volatile* final my_value = &(_req_ary[iThread]);
                *my_value = inValue;

                do {
                        boolean is_cas = false;
                        if(lock_fc(_fc_lock, is_cas)) {
                                //sched_start(iThread);
                                flat_combining();
                                _fc_lock.set(0);
                                //sched_stop(iThread);
                                return true;
                        } else {
                                Memory::write_barrier();
                                while(FCBase<T>::_NULL_VALUE != *my_value && 0 != _fc_lock.getNotSafe()) {
                                        thread_wait(FCBase<T>::_NUM_THREADS, true);
                                } 
                                Memory::read_barrier();
                                if(FCBase<T>::_NULL_VALUE == *my_value)
                                        return true;
                        }
                } while(true);
        }

        //pop ......................................................
        PtrNode<T>* remove(final int iThread, PtrNode<T>* final inPtr) {
                final FCIntPtr inValue = (FCIntPtr) inPtr;

                FCIntPtr volatile* final my_value = &(_req_ary[iThread]);
                *my_value = FCBase<T>::_DEQ_VALUE;

                do {
                        boolean is_cas = false;
                        if(lock_fc(_fc_lock, is_cas)) {
                                //sched_start(iThread);
                                flat_combining();
                                _fc_lock.set(0);
                                //sched_stop(iThread);
                                return (PtrNode<T>*) -*my_value;
                        } else {
                                Memory::write_barrier();
                                while (FCBase<T>::_DEQ_VALUE == *my_value && 0 != _fc_lock.getNotSafe()) {
                                        thread_wait(FCBase<T>::_NUM_THREADS, true);
                                } 
                                Memory::read_barrier();
                                if(FCBase<T>::_DEQ_VALUE != *my_value)
                                        return (PtrNode<T>*) -*my_value;
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
                return "FCStack";
        }


};

#endif
