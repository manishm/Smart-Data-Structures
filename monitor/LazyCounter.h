#ifndef __LAZY_COUNTER__
#define __LAZY_COUNTER__

////////////////////////////////////////////////////////////////////////////////
// File    : LazyCounter.h                                                        
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com                
// Written : 4 March 2011 
//                                                                              
// Fast concurrent counter
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

#include "portable_defns.h"
#include "cpp_framework.h"
#include "Monitor.h"
#include "FCBase.h"


class LazyCounter : public Monitor, public FCBase<FCIntPtr> {

private:

        int              _nthreads        ATTRIBUTE_CACHE_ALIGNED;
        bool             _concurrent;
        _u64*            _buffer;
        volatile _u64*   _counters;

        //backoff settings in nanoseconds
        static final _u64  BACKOFF_START  = 100;
        static final _u64  BACKOFF_MAX    = 1600;

        char pad                          ATTRIBUTE_CACHE_ALIGNED;

public:

        //---------------------------------------------------------------------------
        // LazyCounter API
        //---------------------------------------------------------------------------

        LazyCounter(int nthreads, bool concurrent = true) 
	: _nthreads(nthreads),
          _concurrent(concurrent)
        {
	        if ( !concurrent )
		        _nthreads = 1;

	        //each thread gets a counter on a separate cache line
	        _buffer         = (_u64*)  malloc(CACHE_LINE_SIZE*_nthreads + CACHE_LINE_SIZE);
                _u64 intbuffer =  (((((_u64) _buffer) + CACHE_LINE_SIZE - 1) / CACHE_LINE_SIZE) * CACHE_LINE_SIZE);
                _counters = (volatile _u64*) intbuffer;

                for(int i = 0; i < _nthreads; i++)
		        _counters[i * CACHE_LINE_SIZE / sizeof(_u64)] = 0;

		CCP::Memory::read_write_barrier();
        }

        ~LazyCounter() 
        {
	        free(_buffer);
        }

        inline void increment(int tid, _u64 iAmt)
	{
	        //assume each counter is only written by one thread at a time
                //there is one counter per cacheline so atomic ops are not needed
	        tid = _concurrent ? tid : 0;
	        _counters[tid * CACHE_LINE_SIZE / sizeof(_u64)] += iAmt;
	}

        inline void reset(int tid)
	{
	        //assume each counter is only written by one thread at a time
                //there is one counter per cacheline so atomic ops are not needed
	        tid = _concurrent ? tid : 0;
	        _counters[tid * CACHE_LINE_SIZE / sizeof(_u64)] = 0;
	}

        inline _u64 get()
	{
                //TODO. next-line prefetch hints might be very useful
	        _u64 sum = 0;
	        for(int i = 0; i < _nthreads; i++) {
		        sum += _counters[i * CACHE_LINE_SIZE / sizeof(_u64)];
		}

                return sum;
	}

        inline _u64 waitchange(_u64 lastval)
	{
                //check if changed. if so return new val.
	        _u64 newval = get();
                if ( newval != lastval )
                        return newval;

                //spin with backoff until val changes
                _u64 backoff = BACKOFF_START; 
                do {
                        CCP::Thread::delay(backoff);
                        backoff <<= 1;
                        if ( backoff > BACKOFF_MAX )
			        return get();

                        newval = get();
                } while( newval == lastval );

                //return new val
                return newval;
	}

        inline _u64 spinchange(_u64 lastval)
	{
	        _u64 newval;
                do {
                        newval = get();
                } while( newval == lastval );

                return newval;
	}

        //---------------------------------------------------------------------------
        // Monitor API
        //---------------------------------------------------------------------------

        inline _u64 waitrewardnotsafe(_u64 lastval) { return waitchange(lastval); }

        inline _u64 getrewardnotsafe() { return get(); }

        inline _u64 getreward() { return get(); }

        inline void addrewardnotsafe(int tid, _u64 amt) { increment(tid, amt); }

        inline void addreward(int tid, _u64 amt) { increment(tid, amt); }


        //---------------------------------------------------------------------------
        // Hacky Benchmark API
        //---------------------------------------------------------------------------

        boolean add(final int iThread, PtrNode<FCIntPtr>* final inPtr) {
	        addreward(iThread, 1);
                return true;
	}

        PtrNode<FCIntPtr>* remove(final int iThread, PtrNode<FCIntPtr>* final inPtr) {
	        addreward(iThread, U64(-1));
                return (PtrNode<FCIntPtr>*) 1;
	}

        PtrNode<FCIntPtr>* contain(final int iThread, PtrNode<FCIntPtr>* final inPtr) {
	        _u64 rv = waitrewardnotsafe( (_u64) inPtr );
                return (PtrNode<FCIntPtr>*) rv;
	}    

        int size() {
                return 0;
        }

        final char* name() {
                return "lazycounter";
        }

};


#endif
