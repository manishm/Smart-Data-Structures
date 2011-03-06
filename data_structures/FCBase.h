#ifndef __FCBASE__
#define __FCBASE__

#define TIME_BASED_WAIT

////////////////////////////////////////////////////////////////////////////////
// File    : FCBase.h
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

#include "cpp_framework.h"
#include "Node.h"

#if defined(SPARC64) || defined(SPARC)          
        #include <schedctl.h>
#endif

extern int _gNumThreads;
extern int _gIsDedicatedMode;
extern volatile int _gIsStopThreads;


//these assume 64-bit FCIntPtr
#define _FC_NULL_VALUE (S64(0));
#define _FC_DEQ_VALUE  (S64(0x8000000000000000)+2);
#define _FC_MIN_INT    (S64(0x8000000000000000));
#define _FC_MAX_INT    (S64(0x7FFFFFFFFFFFFFFF));



struct CasInfo {
        int  _failed  ATTRIBUTE_CACHE_ALIGNED;
        int  _succ;
        int  _locks;
        int  _ops;
        char _pad2    ATTRIBUTE_CACHE_ALIGNED;

        CasInfo() {
                _failed = 0;
                _succ   = 0;
                _locks  = 0;
                _ops    = 0;
        }

        void reset() {
                _failed = 0;
                _succ   = 0;
                _locks  = 0;
                _ops    = 0;
        }
};


struct SlotInfo {
        //here 1 can post the request and wait for answer
        FCIntPtr volatile      _req_ans       ATTRIBUTE_CACHE_ALIGNED; 
        int volatile           _time_stamp; 
        SlotInfo* volatile     _next; 
        SlotInfo* volatile     _prev;
        void*                  _custem_info;
        bool                   _deq_pending;
        char                   _pad           ATTRIBUTE_CACHE_ALIGNED;

        SlotInfo() {
                _req_ans     = _FC_NULL_VALUE;
                _time_stamp  = 0;
                _next        = null;
                _prev        = null;
                _custem_info = null;
                _deq_pending = false;
        }
};

template <class T>
class FCBase {
public:

        static volatile int _num_post_read_write       ATTRIBUTE_CACHE_ALIGNED;
        static int          _num_passes;
        static int          _sync_interval;
        static int          _enable_scancount_tuning;
        static int          _enable_lock_scheduling;
        static int          _dynamic_work_size;
        static int          _dynamic_work_intervals;

        static final FCIntPtr _NULL_VALUE = _FC_NULL_VALUE;
        static final FCIntPtr _DEQ_VALUE  = _FC_DEQ_VALUE;
        static final FCIntPtr _MIN_INT    = _FC_MIN_INT;
        static final FCIntPtr _MAX_INT    = _FC_MAX_INT;

        CCP::AtomicInteger  _barr  ATTRIBUTE_CACHE_ALIGNED;

protected:

        //constants -----------------------------------
        final int         _NUM_THREADS  ATTRIBUTE_CACHE_ALIGNED;
        static final int  _MAX_THREADS  = 1024;
        final boolean     _IS_USE_CONDITION;
        int               _sync_count; 

        CasInfo           _cas_info_ary[_MAX_THREADS];
        int               _cpu_cash_contamination[8*1024*1024];
        int               _iTry[_MAX_THREADS*CACHE_LINE_SIZE];



        //helper function -----------------------------

        #if defined(SPARC) || defined(SPARC64)          
                schedctl_t* _schedctl_ary[_MAX_THREADS];

                inline void init_architecture_specific() {
                        for(int iReq=0; iReq < _MAX_THREADS; ++iReq)
                                _schedctl_ary[iReq] = null;
                }

                inline void machine_start_fc(final int iThread) {
                        if(null == _schedctl_ary[iThread])
                                _schedctl_ary[iThread] = schedctl_init();
                        schedctl_start(_schedctl_ary[iThread]);
                }

                inline void machine_end_fc(final int iThread) {
                        schedctl_stop(_schedctl_ary[iThread]);
                }

                inline void thread_wait(final int iThread, final boolean is_short=false) {
                        CCP::Memory::read_write_barrier();
                        CCP::Memory::read_write_barrier();
                        CCP::Memory::read_write_barrier();
                        CCP::Memory::read_write_barrier();
                        CCP::Memory::read_write_barrier();
                        CCP::Memory::read_write_barrier();
                }

                #define lock_fc(x,b) ( 0 == x.getNotSafe() && (b=x.compareAndSet(0, 0xFF)) )

        #else

                inline void init_architecture_specific() { }

                inline void machine_start_fc(final int iThread) { }

                inline void machine_end_fc(final int iThread) { }

                inline void thread_wait(final int iThread, final boolean is_short=false) { 
                        //JME: commented this out to compare apples to apples
                        //if(_NUM_THREADS < 8) {
                        //      CCP::Thread::sleep(0,1);
                        //} else {
                                //this seems to perform better
                                CCP::Thread::yield();
                        //}
                }

                #define lock_fc(x,b) ( 0 == x.getNotSafe() &&  0 == x.getNotSafe() &&  0 == x.getNotSafe() && (b=x.compareAndSet(0, 0xFF)) )
                //#define lock_fc(x,b) (b=x.compareAndSet(0, 0xFF))

        #endif

        //list inner types ------------------------------

        //list fields -----------------------------------
        CCP::ThreadLocal<SlotInfo*>    _tls_slot_info;
        CCP::AtomicReference<SlotInfo> _tail_slot;
        CCP::AtomicReference<SlotInfo> _head_slot;
        int volatile                   _timestamp;

        //list helper function --------------------------
        void init_slot_list() {
                SlotInfo* tmp = new SlotInfo();
                _tail_slot.set(tmp);
                _head_slot.set(tmp);
        }

        void deinit_slot_list() {
                //needed so we see write to _prev in get_new_slot or enq_slot
                CCP::Memory::read_write_barrier();

                SlotInfo *slot = _head_slot.get();
                do {
                        SlotInfo *tmp = slot;
                        slot = slot->_prev;
                        delete tmp;
                } while ( slot != null );
        }

        SlotInfo* get_new_slot() {
                SlotInfo* my_slot= new SlotInfo();
                _tls_slot_info.set(my_slot);

                SlotInfo* curr_tail;
                do {
                        curr_tail = _tail_slot.get();
                        my_slot->_next = curr_tail;
                } while(false == _tail_slot.compareAndSet(curr_tail, my_slot));
		curr_tail->_prev = my_slot;

                return my_slot;
        }

        void enq_slot(SlotInfo* p_slot) {
                SlotInfo* curr_tail;
                do {
                        curr_tail = _tail_slot.get();
                        p_slot->_next = curr_tail;
                } while(false == _tail_slot.compareAndSet(curr_tail, p_slot));
		curr_tail->_prev = p_slot;
        }

        void enq_slot_if_needed(SlotInfo* p_slot) {
                if(null == p_slot->_next) {
                        enq_slot(p_slot);
                }
        }

        int calc_parity(int val) {
           int parity = 0;
           for(int i = 0; i <= (32-1); i++)
           {
              parity ^= val & 1;
              val = val >> 1;
           }
           return parity;
        }

public:

        FCBase( final int num_threads = _gNumThreads, final boolean is_use_condition = false) 
        :       _NUM_THREADS(num_threads), 
                _IS_USE_CONDITION(is_use_condition),
                _sync_count(0)
        {
                init_architecture_specific();
                init_slot_list();
        }

        virtual ~FCBase() 
        {
                deinit_slot_list();
        }
        
        virtual void cas_reset(final int iThread) {
                _cas_info_ary[iThread].reset();;
        }

        virtual void print_cas() {
        }

        virtual void print_custom() {
                int failed = 0;
                int succ = 0;
                int ops = 0;
                int locks = 0;

                for (int i=0; i<_NUM_THREADS; ++i) {
                        failed += _cas_info_ary[i]._failed;
                        succ += _cas_info_ary[i]._succ;
                        ops += _cas_info_ary[i]._ops;
                        locks += _cas_info_ary[i]._locks;
                }
                printf(" 0 0 0 0 0 0 ( %d, %d, %d, %d, %d )", ops, locks, succ, failed, failed+succ);
        }

        virtual int post_computation(final int iThread) {

                /*
                if ( ++_sync_count == _sync_interval )
                {
                        int gen = _barr.getAndIncrement() / _NUM_THREADS;
                        while( (0 == _gIsStopThreads) && (gen <= (_barr.get() / _NUM_THREADS)) );
                        _sync_count = 0;
                }
                */

#ifdef TIME_BASED_WAIT
                if ( _num_post_read_write > 0 )
                {
#if 0
                        //this seems risky given nanosleeps weak guarantees
                        //but will be preferable if the guarantees get better
                        CCP::Thread::sleep(0,_num_post_read_write);
#else
                        //FIXME: this will need tuning for different machines
                        //however, it will be inaccurate in a not bad way if not tuned
                        //since all data structures use it. the time units will simply
                        //be of the wrong scale
                        //++_iTry[iThread*CACHE_LINE_SIZE];
                        _u64 its = (_num_post_read_write * 2 + 2) / 3; //factor
                        _u64 val = 0;
                        for(int i = 0; i < its; ++i)
                        {
                                val +=  val ^ U64(0xAAAA5555AAAA5555);
                        }
                        if ( val == 0 )
                                printf("val = %d", val);
#endif
                }
                return 0;
#else
                int sum=1;
                if(_num_post_read_write > 0) {
                        ++_iTry[iThread*CACHE_LINE_SIZE];
                        unsigned long start_indx = ((unsigned long)(_iTry[iThread*CACHE_LINE_SIZE] * (iThread+1) * 17777675))%(7*1024*1024);

                        for (unsigned long i=start_indx; i<start_indx+_num_post_read_write; ++i) {
                                sum += _cpu_cash_contamination[i];
                                _cpu_cash_contamination[i] =  sum;
                        }
                }
                return sum;
#endif
        }

        //..........................................................................
        virtual boolean add(final int iThread, PtrNode<T>* final inPtr) = 0;
        virtual PtrNode<T>* remove(final int iThread, PtrNode<T>* final inPtr) = 0;
        virtual PtrNode<T>* contain(final int iThread, PtrNode<T>* final inPtr) = 0;    

        //..........................................................................
        virtual int size() = 0;
        virtual final char* name() = 0;
};


#endif
