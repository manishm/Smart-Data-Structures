#ifndef __MUTEX_PAIRINGHEAP__
#define __MUTEX_PAIRINGHEAP__

////////////////////////////////////////////////////////////////////////////////
// File    : MutexPairingHeap.h
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com
// Written : 21 April 2011
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
#include "PairingHeap.h"
#include <pthread.h>

using namespace CCP;


template <class T> class MutexPairHeap  : public FCBase<T> {

private:
        volatile int          _count;
        pthread_mutex_t       _mutex;
        PairHeap<T>           _heap;

public:
        MutexPairHeap() {
                _count = 0;
		pthread_mutex_init(&_mutex, null);
        }

        ~MutexPairHeap() {}

        boolean add(final int iThread, PtrNode<T>* final inPtr) { 
 
		pthread_mutex_lock(&_mutex);

                _heap.insert(inPtr);
                ++_count;

		pthread_mutex_unlock(&_mutex);

                return true;
        }

        PtrNode<T>* remove(final int iThread, PtrNode<T>* final inPtr) { 

	        PtrNode<T>* rv = FCBase<T>::_NULL_VALUE; 

		pthread_mutex_lock(&_mutex);

                if ( _count ) {
		        rv = _heap.deleteMin();
                        --_count;
		}

		pthread_mutex_unlock(&_mutex);

                return rv;
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
                return "MutexPairHeap";
        }

};

#endif

