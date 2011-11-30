#ifndef __HEARTBEAT__
#define __HEARTBEAT__

////////////////////////////////////////////////////////////////////////////////
// File    : Heartbeat.h                                                        
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com                 
//           Henry Hoffman     email: hank@.mit.edu                             
//           David Wingate     email: wingated@mit.edu                          
// Written : 16 February 2011                                                   
//                                                                              
// A much simplified version of Heartbeats                                    
//                                                                              
// Copyright (C) 2011 Jonathan Eastep, Henry Hoffman, David Wingate             
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


class Hb: public Monitor, public FCBase<FCIntPtr> {

private:

        //mode. non-concurrent mode enables optimization
        final bool         concurrent     ATTRIBUTE_CACHE_ALIGNED;

        //written atomically: via write in !concurrent mode otherwise via atomic ops 
        volatile _u64      total_items    ATTRIBUTE_CACHE_ALIGNED; 
        volatile _u64      total_changes;

        //backoff settings in nanoseconds
        static final _u64  BACKOFF_START  ATTRIBUTE_CACHE_ALIGNED  = 100;
        static final _u64  BACKOFF_MAX                             = 6400;

        char               pad            ATTRIBUTE_CACHE_ALIGNED;

public:

        //---------------------------------------------------------------------------
        // Heartbeats API
        //--------------------------------------------------------------------------- 

        Hb(bool concurrent = true)
	: concurrent(concurrent),
	  total_items(0), total_changes(0)
        {
                CCP::Memory::read_write_barrier();
        }

        ~Hb() {}

        inline _u64 waitchange(_u64& changes)
	{
                //check if changed. if so return newval.
	        _u64 lastchanges = changes;
	        _u64 newval = readheartbeatsnotsafe(changes);
                if ( changes != lastchanges )
                        return newval;

                //spin with backoff until val changes
                _u64 backoff = BACKOFF_START; 
                do {
                        CCP::Thread::delay(backoff);
                        backoff <<= 1;
                        if ( backoff > BACKOFF_MAX )
			        backoff = BACKOFF_START;
                        newval = readheartbeatsnotsafe(changes);
                } while( changes == lastchanges );

                return newval;
	}

        inline _u64 readheartbeatsnotsafe() {
	        return total_items;
        }

        inline _u64 readheartbeats() {
	        return concurrent ? FAADD(&total_items, 0) : total_items;  
        }

        inline _u64 readheartbeatsnotsafe(_u64& changes) {
	        changes = total_changes;
	        return total_items;
        }

        inline _u64 readheartbeats(_u64& changes) {
	        changes = concurrent ? FAADD(&total_changes, 0) : total_changes;
	        return concurrent ? FAADD(&total_items, 0) : total_items;  
        }

        inline void heartbeatnotsafe(int num_beats = 1) {
	        if ( num_beats != 0 ) {
	                total_changes += 1;
	                total_items += num_beats;
		}
        }  

        inline void heartbeat(int num_beats = 1) {
	        if ( num_beats != 0 ) {
	                if ( concurrent ) {
		                FAADD(&total_changes, 1);
                                FAADD(&total_items, num_beats);     
                        } else {
		                total_changes += 1;
		                total_items += num_beats;
		        }
                }
        }
        

        //---------------------------------------------------------------------------
        // Monitor API
        //---------------------------------------------------------------------------

        inline _u64 waitchangenotsafe(_u64& changes) { return waitchange(changes); }

        inline _u64 getrewardnotsafe() { return readheartbeatsnotsafe(); }

        inline _u64 getreward() { return readheartbeats(); }

        inline _u64 getrewardnotsafe(_u64& changes) { return readheartbeatsnotsafe(changes); }

        inline _u64 getreward(_u64& changes) { return readheartbeats(changes); }

        inline void addrewardnotsafe(int tid, _u64 amt) { heartbeatnotsafe(amt); }

        inline void addreward(int tid, _u64 amt) { heartbeat(amt); }


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
	        _u64 rv = (_u64) inPtr;
	        waitchangenotsafe(rv);
                return (PtrNode<FCIntPtr>*) rv;
	}    

        int size() {
                return 0;
        }

        final char* name() {
                return "heartbeat";
        }


};


#endif
