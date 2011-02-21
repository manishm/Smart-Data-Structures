//////////////////////////////////////////////////////////////////////////////// 
// File    : FCBase.cpp  
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


#include "FCBase.h"

template<> volatile int FCBase<FCIntPtr>::_num_post_read_write = 800;
template<> int          FCBase<FCIntPtr>::_num_passes = 14;
template<> int          FCBase<FCIntPtr>::_sync_interval = 0;
template<> int          FCBase<FCIntPtr>::_enable_scancount_tuning = 1;
template<> int          FCBase<FCIntPtr>::_enable_lock_scheduling = 0;
template<> int          FCBase<FCIntPtr>::_dynamic_work_size = 0;
template<> int          FCBase<FCIntPtr>::_dynamic_work_intervals = 0;
template<> double       FCBase<FCIntPtr>::_rl_to_sleepidle_ratio = 1.0;
