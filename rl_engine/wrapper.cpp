                                                        
////////////////////////////////////////////////////////////////////////////////                
// File    : wrapper.cpp                                                               
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com 
//           David Wingate     email: wingated@mit.edu                                          
// Written : 16 February 2011                                                                   
//                                         
// Copyright (C) 2011 Jonathan Eastep, David Wingate 
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


int main() {



   while ( 1 ) {

      rl_nac_action_sample( r );

      for ( int i=0; i<num_threads; i++ ) {
         threads[i].setpriority( vals[i] );
      }

      reward = read_heartbeat_reward();

      rl_nac_update( r, reward, statefeats );
         
   }

}
