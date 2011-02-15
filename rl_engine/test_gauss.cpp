                 
////////////////////////////////////////////////////////////////////////////////                
// File    : test_guass.cpp                                                            
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

double means[2];

rl_act_entry_t raes[2] =
  {
    { RLA_GAUSS, 1, 0, &means[0] },
    { RLA_GAUSS, 1, 0, &means[1] },
  };

rl_act_desc_t rad = { 2, raes };

int main() {
  rl_nac_t r;
  double statefeats[1], reward;
  int pcnt1, pcnt2;
  double *p1, *p2;

  r = rl_nac_init( 1, &rad );
  statefeats[0] = 1;

  rl_nac_get_params( r, 0, &pcnt1, &p1 );
  rl_nac_get_params( r, 1, &pcnt2, &p2 );

  for ( int i=0; i< 50000; i++ ) {
    rl_nac_action_sample( r );

    // you get a reward based on how closely you sampled from 5,5
    reward = -sqrt( (means[0] - 5)*(means[0] - 5) + 
		    (means[1] - 5)*(means[1] - 5) );

    rl_nac_update( r, reward, statefeats );
    printf( "%.2f %.2f %.2f\n", reward, *p1, *p2 );

  }  

}
