#ifndef __MUTEX_QUEUE__
#define __MUTEX_QUEUE__

////////////////////////////////////////////////////////////////////////////////
// File    : MutexQueue.h
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com
// Written : 20 April 2011
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

#include "FCBase.h"
#include "cpp_framework.h"
#include <pthread.h>

using namespace CCP;


template <class T> class MutexQueue  : public FCBase<T> {

private:
	PtrNode<T>* volatile  _head;
	PtrNode<T>* volatile  _tail;
        volatile int          _count;
        pthread_mutex_t       _mutex;

public:
        MutexQueue() {
	        _head = FCBase<T>::_NULL_VALUE;
		_tail = _head;
                _count = 0;
		pthread_mutex_init(&_mutex, null);
        }

        ~MutexQueue() {}

        boolean add(final int iThread, PtrNode<T>* final inPtr) { 
 
		pthread_mutex_lock(&_mutex);

                inPtr->next = FCBase<T>::_NULL_VALUE;
                if ( _count )
		        _tail->next = inPtr;
                else
		        _head = inPtr;
                _tail = inPtr;
                ++_count;
		pthread_mutex_unlock(&_mutex);

                return true;
        }

        PtrNode<T>* remove(final int iThread, PtrNode<T>* final inPtr) { 

	        PtrNode<T>* tmp = FCBase<T>::_NULL_VALUE; 

		pthread_mutex_lock(&_mutex);
                if ( _count ) {
		        tmp = _head;
                        _head = _head->next;
                        --_count;
		}
		pthread_mutex_unlock(&_mutex);

                return tmp;
        }

        PtrNode<T>* contain(final int iThread, PtrNode<T>* final inPtr) { 
                return FCBase<T>::_NULL_VALUE;
        }

        int size() {
                return _count;
        }

        bool empty() {
	        return ( 0 == _count );
	}

        const char* name() {
                return "MutexQueue";
        }

};


#endif

