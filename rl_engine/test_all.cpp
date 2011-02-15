                
////////////////////////////////////////////////////////////////////////////////                
// File    : test_all.cpp                                                          
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

#include "distance.h"

#define NUM_ACTS 5

double means[2];
int perm_vals[5];
int disc_vals[1];
int bin_vals[1];

rl_act_entry_t raes[NUM_ACTS] =
  {
    { RLA_PERM,     5, 0,  perm_vals },
    { RLA_GAUSS,    1, 0,  &means[0] },
    { RLA_GAUSS,    1, 0,  &means[1] },
    { RLA_DISCRETE, 1, 10, disc_vals },
    { RLA_BINARY,   1, 0,  bin_vals },
  };

rl_act_desc_t rad = { NUM_ACTS, raes };

int main() {
  rl_nac_t r;
  double statefeats[1], reward;
  int pcnt1;
  double *params;
  char true_answer[6] = { '0', '1', '2', '3', '4', 0 };
  char cur_answer[100];

  r = rl_nac_init( 1, &rad );
  statefeats[0] = 1;

  for ( int i=0; i< 10000; i++ ) {
    rl_nac_action_sample( r );

    sprintf( cur_answer, "%d%d%d%d%d", 
	     perm_vals[0], perm_vals[1], perm_vals[2],
	     perm_vals[3], perm_vals[4] );

    // you get a reward based on how closely the sampled vector is to
    // the reference one, according to Levenshtein distance.

    reward = 0;

    // perm
    reward += -distance( true_answer, cur_answer );

    // gauss
    reward += -sqrt( (means[0] - 5)*(means[0] - 5) + 
		    (means[1] - 5)*(means[1] - 5) );

    // discrete
    reward += disc_vals[0];

    // binary
    reward += bin_vals[0];

    rl_nac_update( r, reward, statefeats );

    // print out all of the parameters
    printf( "%.2f ", reward );
    for ( int a=0; a<NUM_ACTS; a++ ) {
      rl_nac_get_params( r, a, &pcnt1, &params );
      for ( int j=0; j<pcnt1; j++ ) {
	printf( "%.2f ", params[j] );
      }
    }
    printf( "\n" );

  }  

}
