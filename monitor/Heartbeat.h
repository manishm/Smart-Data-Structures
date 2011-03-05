#ifndef __HEARTBEAT__
#define __HEARTBEAT__

////////////////////////////////////////////////////////////////////////////////
// File    : Heartbeat.h                                                        
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com                 
//           Henry Hoffman     email: hank@.mit.edu                             
//           David Wingate     email: wingated@mit.edu                          
// Written : 16 February 2011                                                   
//                                                                              
// Based on the Application Heartbeats paper                                    
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

#include "Monitor.h"
#include "cpp_framework.h"
#include "portable_defns.h"

class Hb: public Monitor {

private:

        volatile _u64    total_items  ATTRIBUTE_CACHE_ALIGNED; 
        final bool       concurrent   ATTRIBUTE_CACHE_ALIGNED;
        char             pad          ATTRIBUTE_CACHE_ALIGNED;

public:

        Hb(bool concurrent = true)
	: concurrent(concurrent),
	  total_items(0) 
        {
                CCP::Memory::read_write_barrier();
		//std::cerr << "concurrent = " << concurrent << std::endl;
        }

        ~Hb() {}

        inline _u64 getrewardnotsafe() {
                return total_items;
        }

        inline _u64 getreward() {
	        return concurrent ? FAADD(&total_items, 0) : getrewardnotsafe();  
        }

        inline void heartbeatnotsafe( int num_beats = 1 ) {
                total_items += num_beats;     
        }

        inline void heartbeat( int num_beats = 1 ) {
	        if ( concurrent )
                        FAADD(&total_items, num_beats);     
                else
		        heartbeatnotsafe(num_beats);
        }

        inline void addrewardnotsafe( int amt ) {
                heartbeatnotsafe(amt);
        }  

        inline void addreward( int amt ) {
                heartbeat(amt);
        }
};


#endif
