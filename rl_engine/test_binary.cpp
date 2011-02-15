                 
////////////////////////////////////////////////////////////////////////////////                
// File    : test_binary.cpp                                                             
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

#include <stdio.h>
#include <math.h>
#include "rl_agent_c.h"

int vals[1];

rl_act_entry_t raes[1] =
  {
    { RLA_BINARY, 1, 0, vals },
  };

rl_act_desc_t rad = { 1, raes };

int main() {
  rl_nac_t r;
  double statefeats[1], reward;
  int pcnt1;
  double *params;

  r = rl_nac_init( 1, &rad );
  statefeats[0] = 1;

  rl_nac_get_params( r, 0, &pcnt1, &params );

  for ( int i=0; i< 10000; i++ ) {
    rl_nac_action_sample( r );

    // the reward is proportional to the value.
    reward = vals[0];

    rl_nac_update( r, reward, statefeats );

    printf( "%.2f ", reward );
    for ( int j=0; j<pcnt1; j++ ) {
      printf( "%.2f ", params[j] );
    }
    printf( "\n" );

  }  

}
