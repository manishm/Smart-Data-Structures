#ifndef MONITOR_H
#define MONITOR_H

////////////////////////////////////////////////////////////////////////////////
// File    : Monitor.h                                                          
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com                 
// Written : 16 February 2011
//                                                                              
// Generic wrapper for reward monitoring frameworks
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

class Monitor {

public:
        // this returns the total reward that has ever been registered
        virtual _u64 getreward() = 0;

        virtual _u64 getrewardnotsafe() = 0;

        virtual void addreward(int amt) = 0;

        virtual void addrewardnotsafe(int amt) = 0;
};

#endif
