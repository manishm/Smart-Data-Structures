
#ifndef __NODE_H__
#define __NODE_H__

////////////////////////////////////////////////////////////////////////////////
// File    : Node.h                                                             
// Author  : Jonathan Eastep   email: jonathan.eastep@gmail.com                 
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

#include <assert.h>
#include "cpp_framework.h"
#include "portable_defns.h"

template <typename T>
class PtrNode {

protected:
        template<class C> friend class WorkQueueFreeList;
        template<class C> friend class MutexQueue;

        FCIntPtr     key;
        T*           value;
        volatile bool          free;
        PtrNode<T>*  volatile  next;

public:
        PtrNode()
	: key(0), value(null), free(true)
	{
	}

        PtrNode(const FCIntPtr key, T* const value)
        : key(key), value(value), free(true) 
	{
	        assert( ((FCIntPtr) value) >= 0 );
        }

        virtual ~PtrNode() 
        {};

        bool operator<(const PtrNode<T>& rhs) const {
                return key < rhs.key;
        }

        FCIntPtr getkey() {
                return key;
        }

        T* getvalue() {
                return value;
        }
 
        void setkey(FCIntPtr k) {
	        key = k;
	}

        void setvalue(T* v) {
	        value = v;
	}

        void setfree(bool v) {
	        free = v;
	}

        bool getfree() {
	        return free;
	}

};

class FCIntPtrNode : public PtrNode<FCIntPtr> {
public:
        FCIntPtrNode(const FCIntPtr keyandval): PtrNode<FCIntPtr>(keyandval, &key) {
                //there are some restrictions on the integer range
                //this may need to change if additional ops beyond add/remove are supported; base it on allowable user-space pointers
                //FIXME: we need asserts here!!!
                //assert( keyandval > FCBase<FCIntPtr>::_NULL_VALUE );  //if not, will screw up distinguishing ops
                //assert( -keyandval != FCBase<FCIntPtr>::_DEQ_VALUE ); //if not, will screw up distinguishing ops
                //assert( keyandval < FCBase<FCIntPtr>::_MAX_INT );     //if not, may or may not screw up skiplist sentinal nodes
        }

        ~FCIntPtrNode() {}

        FCIntPtr getvalue() {
                return *value;
        }
};



template <class T>
class FreeList  {

private:

        T**                  _free_list;
        volatile int         _head;
        volatile int         _tail;
        volatile int         _size;
        int                  _max_size;
        T*                   _buff;

        static volatile int  _lck;

public:

        FreeList(int threads, int size_per_thread) 
	: _size(threads*size_per_thread)
        {
	        _max_size = _size;
	        _free_list = new T*[_max_size];
                _buff = new T[_max_size]();
                for(int i = 0; i < _max_size; i++)
		        _free_list[i] = &_buff[i];

                _head = 0;
	        _tail = 0;
		CCP::Memory::read_write_barrier();
	}

        FreeList(int threads, int size_per_thread, T* buff) 
	: _buff(NULL), _size(threads*size_per_thread)
        {
		CCP::Memory::read_write_barrier();
	        _max_size = _size;
	        _free_list = new T*[_max_size];
                for(int i = 0; i < _max_size; i++)
		        _free_list[i] = &buff[i];

                _head = 0;
	        _tail = 0;
		CCP::Memory::read_write_barrier();
	}

        ~FreeList() {
	        delete[] _buff;
	        delete[] _free_list;
	        CCP::Memory::read_write_barrier();
	}

        inline_ T* alloc(int iThread) {
	        T* e = null;

	        while( !CAS(&_lck, 0, 1) );

                if ( _size > 0 ) {
                        --_size;

			e = _free_list[_head];

			++_head;
			if ( _head == _max_size )
			        _head = 0;

		}

                FASTORE(&_lck, 0);
                return e;
	}

        inline_ void dealloc(int iThread, T* e) {
	        while( !CAS(&_lck, 0, 1) );

                assert( _size < _max_size );
                ++_size;

                _free_list[_tail] = e;

                ++_tail;
                if ( _tail == _max_size )
		       _tail = 0;
 
                FASTORE(&_lck, 0);
	}

};

template<class T> volatile int FreeList<T>::_lck = 0;




#if 1

//synchronization-free free list
//assumptions:
//  for each element
//    only one thread allocs it
//    only one thread (possibly the same) uses and deallocs it
//implementation:
//  per-thread free lists
//  special 2-way synchronization based on a hand-shake
//  bad programs can overfill the buffer and cause deadlock
template <class T>
class WorkQueueFreeList  {

private:
          
        class AlignedInt {
        public:
	        char  _pad1[CACHE_LINE_SIZE];
	        int   _val;
		char  _pad2[CACHE_LINE_SIZE-sizeof(int)];

	        AlignedInt() : _val(0) {}
	        AlignedInt(int v) :_val(v) {}
	};

        T*            _free_list;
        AlignedInt*   _heads;
        int           _sz_per;
        int           _threads;
        bool          _list_alloc; 

public:

        //static volatile int  _lck;

        WorkQueueFreeList(int threads, int size_per_thread) 
	: _sz_per(size_per_thread), _threads(threads), _list_alloc(true)
        {
	        //std::cout << "size_per_thread= " << size_per_thread << std::endl;
	        //std::cout << "threads= " << threads << std::endl;
	        _free_list = new T[threads*size_per_thread];
                _heads = new AlignedInt[threads];

                for(int i = 0; i < threads; i++)
		        _heads[i] = AlignedInt(i * size_per_thread);

		CCP::Memory::read_write_barrier();
	}

        WorkQueueFreeList(int threads, int size_per_thread, T* buff) 
	: _sz_per(size_per_thread), _threads(threads), _free_list(buff), _list_alloc(false)
        {
		CCP::Memory::read_write_barrier();
	        //std::cout << "size_per_thread= " << size_per_thread << std::endl;
	        //std::cout << "threads= " << threads << std::endl;

                _heads = new AlignedInt[threads];
                for(int i = 0; i < threads; i++)
		        _heads[i] = AlignedInt(i * size_per_thread);

		CCP::Memory::read_write_barrier();
	}

        ~WorkQueueFreeList() {
	        CCP::Memory::read_write_barrier();
                if ( _list_alloc )
	                delete[] _free_list;
	        delete[] _heads;
	}

        inline_ T* alloc(int iThread) {
	        //while( !CAS(&_lck, 0, 1) );

                //test
	        //CCP::Memory::read_write_barrier();

	        int& head = _heads[iThread]._val;
		//std::cout << "head = " << head << std::endl;
	        T* e = &_free_list[head];

                while( !e->free );

                ++head;
                if ( head == ((iThread + 1) * _sz_per) )
		        head -= _sz_per;

                e->free = false;

                //test
		//CCP::Memory::read_write_barrier();

                //FASTORE(&_lck, 0);
                return e;
	}

        inline_ void dealloc(int iThread, T* e) {
	        //while( !CAS(&_lck, 0, 1) );

                //test
	        //CCP::Memory::read_write_barrier();

	        while( e->free );
                e->free = true;

                //test
		//CCP::Memory::read_write_barrier();

                //FASTORE(&_lck, 0);
	}

};

//template<class T> volatile int WorkQueueFreeList<T>::_lck = 0;


#else

typedef FreeList WorkQueueFreeList;

#endif

#endif


