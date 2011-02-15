////////////////////////////////////////////////////////////////////////////////             
// File    : LearningEngine.cpp                                                                      
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
// MERCHANTABILITY99 or FITNESS FOR A PARTICULAR PURPOSE. See the GNU                        
// General Public License for more details.                                                  
//                                                                                           
// You should have received a copy of the GNU General Public License                         
// along with this program; if not, write to the Free Software Foundation  
// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA                    
////////////////////////////////////////////////////////////////////////////////     

#include "LearningEngine.h"

volatile unsigned int         LearningEngine::refcount = 0;
volatile unsigned int         LearningEngine::htlock = 0;
volatile unsigned int         LearningEngine::signalquit = 0;
volatile unsigned int         LearningEngine::lelock = 0;
volatile unsigned int         LearningEngine::lelistpos = 0;
std::vector<LearningEngine*>  LearningEngine::lelist;
pthread_t                     LearningEngine::managerthread;
volatile _u64                 LearningEngine::objstatus = 0;


