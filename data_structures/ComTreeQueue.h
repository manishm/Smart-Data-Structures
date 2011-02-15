#ifndef __COM_TREE_QUEUE__
#define __COM_TREE_QUEUE__

//FIXME: This seems to have a bug in it
//I don't understand the placement of threads in the array-based tree
//Also don't understand if two threads are supposed to share a node or not
//the LOG in the constructor looks wrong too


////////////////////////////////////////////////////////////////////////////////
// File    : ComTreeQueue.h
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com
// Author  : Ms.Moran Tzafrir  email: morantza@gmail.com
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

using namespace CCP;

template <class T> class ComTreeQueue : public FCBase<T> {
private:

        //-----------------------------------------------------
        class TTASLock {
        public:
                char               _pad1[CACHE_LINE_SIZE];
                VolatileType<_u32> _lock;
                char               _pad2[CACHE_LINE_SIZE];

                TTASLock() : _lock(0) {}
                ~TTASLock() {_lock=0;}

                inline_ void init() {_lock = 0;}

                inline_ void lock(CasInfo& my_cas_info) {
                        do {
                                if(0 == _lock && (0 == MUTEX_ENTER(&_lock))) {
                                        ++(my_cas_info._succ);
                                        return;
                                } else {
                                        ++(my_cas_info._failed);
                                }
                        } while(true);
                }

                inline_ bool tryLock(CasInfo& my_cas_info) {
                        if ( 0 == MUTEX_ENTER(&_lock) ) {
                                ++(my_cas_info._succ);
                                return true;
                        } else {
                                ++(my_cas_info._failed);
                                return false;
                        }
                }

                inline_ bool isLocked() {
                        return 0 != _lock;
                }

                inline_ void unlock() {
                        MUTEX_EXIT(&_lock);
                }
        };

        //inner classes queue ---------------------------------
        struct Node {
                char                 _pad1[CACHE_LINE_SIZE];
                final FCIntPtr       _value;
                Node* volatile       _next;
                int                  d1, d2, d3, d4, d5, d6, d7;
                int                  d21, d22, d23, d24, d25, d26, d27;
                char                 _pad2[CACHE_LINE_SIZE];

                Node(final FCIntPtr x) : _value(x) {
                        _next = null;
                }
        };
        //fields queue ----------------------------------------
        char           _pad1[CACHE_LINE_SIZE];
        Node* volatile _head;
        Node* volatile _tail;
        char           _pad2[CACHE_LINE_SIZE];

        //inner class -----------------------------------------
        #define REQUEST_NONE  (0)
        #define REQUEST_LEFT  (4)
        #define REQUEST_RIGHT (8)

        struct ThreadInfo {
                char                            _pad1[CACHE_LINE_SIZE];
                int                             _path_length;
                int                             _path_ary[128];
                int                             _add_req_ary[128];

                VolatileType<int>               _num_op;
                VolatileType<int>               _num_deq;
                VolatileType<FCIntPtr>          _op_ary[128];
                VolatileType<ThreadInfo*>       _removers[128];
                char                            _pad2[CACHE_LINE_SIZE];
                VolatileType<FCIntPtr>          _rv;
                char                            _pad3[CACHE_LINE_SIZE];

                ThreadInfo() {
                        init();
                }

                void init() {
                        _path_length    = 0;
                        _num_op         = 0;
                        _num_deq        = 0;
                        _rv             = FCBase<T>::_NULL_VALUE;

                        for (int i=0; i<128; ++i) {
                                _path_ary[i] = 0;
                                _add_req_ary[i] = REQUEST_NONE;
                                _op_ary[i] = 0;
                                _removers[i] = null;
                        }
                }

                void add_thread_info(ThreadInfo* threadInfo) {
                        final int copy_op = threadInfo->_num_op;
                        if((_num_op + copy_op) >= 128) {
                                System_err_println("*******************************************");
                        }
                        for (int i=0; i<copy_op; ++i) {
                                _op_ary[_num_op+i] = threadInfo->_op_ary[i];
                        }
                        _num_op += copy_op; 

                        final int copy_deq = threadInfo->_num_deq;
                        if((_num_deq + copy_deq) >= 128) {
                                System_err_println("*******************************************");
                        }
                        for (int i=0; i<copy_deq; ++i) {
                                _removers[_num_deq+i] = threadInfo->_removers[i];
                        }
                        _num_deq += copy_deq; 
                }
        };

        struct TreeNode {
                char pad1[CACHE_LINE_SIZE];
                TTASLock                        _lock;
                VolatileType<ThreadInfo*>       _left;
                VolatileType<ThreadInfo*>       _right;
                char pad2[CACHE_LINE_SIZE];

                TreeNode() : _left(null), _right(null)
                { }
        };

        //helper function -------------------------------------
        inline_ int NearestPowerOfTwo(final int x) {
                int mask = 1;
                while(mask < x) {
                        mask <<= 1;
                }
                return mask;
        }
        inline_ int left(final int i)                                   { return ((2*i) + 1); } 
        inline_ int right(final int i)                                  { return ((2*i) + 2); } 
        inline_ int parent(final int i)                                 { return ((i-1)/2);   } 
        inline_ int get_thread_indx(final int iThread)  { 
                return (_NUM_NODES - _MAX_THREADS + (iThread/2)); 
                //return (_NUM_NODES/2) + iThread;
        }
        inline_ bool get_thread_is_left(final int iThread) { return (0 != (iThread%2)); }
        inline_ bool is_parent_left(final int i) { return (0 != (i%2)); }
        //inline_ bool is_parent_left(final int i) { return (0 != (parent(i)%2));       }

        //fields ----------------------------------------------
        char pad3[CACHE_LINE_SIZE];
        final int               _MAX_THREADS;
        final int               _LOG_THREADS;
        final int               _NUM_NODES;
        TreeNode* final _tree;
        ThreadInfo*             _thread_info_ary;
        char pad4[CACHE_LINE_SIZE];

public:

        ComTreeQueue() 
          :     _MAX_THREADS( NearestPowerOfTwo(FCBase<T>::_NUM_THREADS/2) ),
                //_MAX_THREADS( NearestPowerOfTwo(FCBase<T>::_NUM_THREADS) ),   
                //_LOG_THREADS( (int) ceil(log((double)FCBase<T>::_NUM_THREADS)) ),
                _LOG_THREADS( (int) ceil( (double) Integer::log2(FCBase<T>::_NUM_THREADS) ) ),
                _NUM_NODES( (int) (pow((double)2, (double)(_LOG_THREADS+1))) ),
                _tree( new TreeNode[_NUM_NODES] )
        {
                //cerr << "_NUM_THREADS=" << FCBase<T>::_NUM_THREADS << " _MAX_THREADS=" 
                //     << _MAX_THREADS << " _LOG_THREADS=" << _LOG_THREADS << " _NUM_NODES=" 
                //     << _NUM_NODES << endl;
                //exit(0);
                //init queue ............................................
                // Allocate a free node
                // Make it the only node in the linked list
                // Both Head and Tail point to it
                Node* const new_node = new Node(FCBase<T>::_NULL_VALUE);        
                new_node->_next = null;                                         
                _head = _tail = new_node;                                       

                //.......................................................
                _thread_info_ary = new ThreadInfo[FCBase<T>::_NUM_THREADS];

        }

        void execute_requests(ThreadInfo* thread_info, CasInfo& my_cas_info) {
                // ......................................................
                int num_added = 0;
                Memory::read_barrier();

                int deq_indx = 0;

                for (int i=0; i<thread_info->_num_op; ++i) {
                        final FCIntPtr curr_value = thread_info->_op_ary[i];

                        //enq
                        if(curr_value > FCBase<T>::_NULL_VALUE) {                                               
                                Node* const new_node = new Node(curr_value);    // Allocate a new node from the free list
                                new_node->_next = null;                         // Set next pointer of node to NULL
                                _tail->_next = new_node;                        // Link node at the end of the linked list
                                _tail = new_node;                               // Swing Tail to node
                                ++(my_cas_info._ops);
                        } else { //deq

                                //if (curr_value != FCBase<T>::_NULL_VALUE) {
                                if (curr_value != FCBase<T>::_DEQ_VALUE) {
                                        cerr << "Oops. got weird op= " << curr_value << endl;
                                        exit(0);
                                }

                                if ( deq_indx >= thread_info->_num_deq )
                                {
                                        cerr << "Oops. too many deq. deq_indx= " << deq_indx << " num_deq= " << thread_info->_num_deq << endl;
                                        exit(0);
                                }

                                Node* const old_head = _head;                           // Read Head
                                Node* const new_head = old_head->_next;                 // Read next pointer
                                if(null != new_head) {                                          // Is queue empty?
                                        final FCIntPtr curr_value = new_head->_value;      // Queue not empty.  Read value before release
                                        if ( thread_info->_removers[deq_indx] == null ) {
                                          cerr << "oops..." << endl;
                                          exit(0);
                                        }
                                        thread_info->_removers[deq_indx++]->_rv = curr_value;
                                        _head = new_head;                                       // Swing Head to next node
                                } else {
                                        if ( thread_info->_removers[deq_indx] == null ) {
                                          cerr << "oops..." << endl;
                                          exit(0);
                                        }
                                        thread_info->_removers[deq_indx++]->_rv = FCBase<T>::_NULL_VALUE;
                                }
                                ++(my_cas_info._ops);
                                //}
                        }
                }
        }

        boolean add(final int iThread, PtrNode<T>* final inPtr) { 
                final FCIntPtr inValue = (FCIntPtr) inPtr;

                CasInfo&        my_cas_info = FCBase<T>::_cas_info_ary[iThread];
                ThreadInfo* thread_info = &(_thread_info_ary[iThread]);
                thread_info->init();

                //Add ENQ request ...........................................
                ++(thread_info->_num_op);
                thread_info->_op_ary[0] = inValue;
                if (inValue == 0 ) {
                        cerr << "oops. got bad inValue" << endl;
                        exit(0);
                }

                int i_node = get_thread_indx(iThread);
                bool is_left = get_thread_is_left(iThread);


                //lock all my path from leaf to root ........................
                Memory::read_write_barrier();
                bool is_exec;
                if (is_left) {
                        //if ( _tree[i_node]._left != null ) {
                        //        cerr << "wth? got " << _tree[i_node]._right << " btw, thread_info= " << thread_info << endl;
                        //      exit(0);
                        //}

                        _tree[i_node]._left = thread_info;
                }
                else {
                        //if ( _tree[i_node]._right != null ) {
                        //        cerr << "wth? got " << _tree[i_node]._right << " btw, thread_info= " << thread_info << endl;
                        //      exit(0);
                        //}

                        _tree[i_node]._right = thread_info;
                }

                do {
                        is_exec=true;

                        //climb the tree as much as we can ......................
                        if (_tree[i_node]._lock.tryLock(my_cas_info)) {
                                final int i_path = thread_info->_path_length;
                                ++(thread_info->_path_length);
                                thread_info->_path_ary[i_path] = i_node;

                                //we have the lock, take the request and climb-up....
                                if (is_left) {
                                        if (null != _tree[i_node]._right) {
                                                thread_info->add_thread_info(_tree[i_node]._right);
                                                thread_info->_add_req_ary[i_path] = REQUEST_RIGHT;
                                        }
                                } else {
                                        if (null != _tree[i_node]._left) {
                                                thread_info->add_thread_info(_tree[i_node]._left);
                                                thread_info->_add_req_ary[i_path] = REQUEST_LEFT;
                                        }
                                }
                                if(0==i_node) //if root
                                        break;
                                is_left = is_parent_left(i_node);
                                i_node = parent(i_node);
                        } else {
                                //wait for lock .....................................
                                do {
                                        //thread_wait(iThread);
                                } while (_tree[i_node]._lock.isLocked());

                                if (is_left) {
                                        if(null == _tree[i_node]._left) { //our requests been answered
                                                is_exec=false;
                                                break;
                                        }
                                } else {
                                        if(null == _tree[i_node]._right) { //our requests been answered
                                                is_exec=false;
                                                break;
                                        }
                                }
                        }
                } while(true);

                //got root execute requests .................................
                if(is_exec) {
                        execute_requests(thread_info, my_cas_info);
                }

                //unlock our path & mark answered requests ..................
                for (int i=(thread_info->_path_length-1); i>=0; --i) {
                        final int i_curr_node = thread_info->_path_ary[i];
                        if (REQUEST_LEFT == thread_info->_add_req_ary[i]) {
                                _tree[i_curr_node]._left = null;
                        } else if (REQUEST_RIGHT == thread_info->_add_req_ary[i]) {
                                _tree[i_curr_node]._right = null;
                        }
                        _tree[i_curr_node]._lock.unlock();
                }

                //...........................................................
                return true;
        }

        PtrNode<T>* remove(final int iThread, PtrNode<T>* final inPtr) { 
                final FCIntPtr inValue = (FCIntPtr) inPtr;

                CasInfo&        my_cas_info = FCBase<T>::_cas_info_ary[iThread];
                ThreadInfo* thread_info = &(_thread_info_ary[iThread]);
                thread_info->init();

                //Add DEQ request ...........................................
                int i_node = get_thread_indx(iThread);
                bool is_left = get_thread_is_left(iThread);
                ++(thread_info->_num_deq);
                thread_info->_removers[0] = thread_info;
                ++(thread_info->_num_op);
                thread_info->_op_ary[0] = FCBase<T>::_DEQ_VALUE;

                //lock all my path from leaf to root ........................
                Memory::write_barrier();
                bool is_exec;
                if (is_left) 
                        _tree[i_node]._left = thread_info;
                else 
                        _tree[i_node]._right = thread_info;
                do {
                        is_exec=true;

                        //climb the tree as much as we can ......................
                        if (_tree[i_node]._lock.tryLock(my_cas_info)) {
                                final int i_path = thread_info->_path_length;
                                ++(thread_info->_path_length);
                                thread_info->_path_ary[i_path] = i_node;

                                //we have the lock, take the request and climb-up....
                                if (is_left) {
                                        if (null != _tree[i_node]._right) {
                                                thread_info->add_thread_info(_tree[i_node]._right);
                                                thread_info->_add_req_ary[i_path] = REQUEST_RIGHT;
                                        }
                                } else {
                                        if (null != _tree[i_node]._left) {
                                                thread_info->add_thread_info(_tree[i_node]._left);
                                                thread_info->_add_req_ary[i_path] = REQUEST_LEFT;
                                        }
                                }
                                if(0==i_node) //if root
                                        break;
                                is_left = is_parent_left(i_node);
                                i_node = parent(i_node);
                        } else {
                                //wait for lock .....................................
                                do {
                                        //thread_wait(iThread);
                                } while (_tree[i_node]._lock.isLocked());

                                if (is_left) {
                                        if(null == _tree[i_node]._left) { //our requests been answered
                                                is_exec=false;
                                                break;
                                        }
                                } else {
                                        if(null == _tree[i_node]._right) { //our requests been answered
                                                is_exec=false;
                                                break;
                                        }
                                }
                        }
                } while(true);

                //got root execute requests .................................
                if(is_exec) {
                        execute_requests(thread_info, my_cas_info);
                }

                //unlock our path & mark answered requests ..................
                for (int i=(thread_info->_path_length-1); i>=0; --i) {
                        final int i_curr_node = thread_info->_path_ary[i];
                        if (REQUEST_LEFT == thread_info->_add_req_ary[i]) {
                                _tree[i_curr_node]._left = null;
                        } else if (REQUEST_RIGHT == thread_info->_add_req_ary[i]) {
                                _tree[i_curr_node]._right = null;
                        }
                        _tree[i_curr_node]._lock.unlock();
                }

                //...........................................................
                return (PtrNode<T>*) thread_info->_rv.get();
        }

        PtrNode<T>* contain(final int iThread, PtrNode<T>* final inPtr) { 
                final FCIntPtr inValue = (FCIntPtr) inPtr;
                return FCBase<T>::_NULL_VALUE;
        }

        int size() {
                return 0;
        }

        const char* name() {
                return "ComTreeQueue";
        }

};


#endif
