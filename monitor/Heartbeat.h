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

        volatile _u64    total_items_since_last_read   ATTRIBUTE_CACHE_ALIGNED;
        volatile _u64    total_items                   ATTRIBUTE_CACHE_ALIGNED; //aligned?
        char             pad[CACHE_LINE_SIZE];

        void init( void ) {
                total_items_since_last_read = 0;
                total_items = 0;
                CCP::Memory::read_write_barrier();
        }

        void finalize( void )
        { }

        double get_lastread_count( void ) {
                _u64 val;
                val = FASTORE( &total_items_since_last_read, 0 );
                //cout << "reward is " << val << endl;
                return val;
        }

        double get_read_count( void ) {
                //fixme: atomic operation may not be required here
                return FAADD( &total_items, 0 );
        }

        double get_read_count_low_overhead( void ) {
                return total_items;
        }

public:

        Hb() {
                init();
        }

        ~Hb() {
                finalize();
        }

        double getreward() {
                return get_read_count();     
        }

        double getrewardlowoverhead() {
                return get_read_count_low_overhead();
        }

        // this returns the reward that has been registered since the
        // last time this function was called. remember this object 
        // can potentially be shared among multiple other objects
        double getrewardsincelast() {
                return get_lastread_count();     
        }

        void heartbeat( int num_beats ) {
                //cout << "adding heartbeat reward" << endl;
                FAADD( &total_items_since_last_read, num_beats );
                FAADD( &total_items, num_beats );     
        }

        void heartbeat() {
                heartbeat(1);     
        }

        void addreward( int rewardcode, double amt ) {
                int iamt = amt;
                heartbeat(iamt);
        }

};


#endif
