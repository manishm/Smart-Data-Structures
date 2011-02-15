#ifndef _RL_AGENT_C_H
#define _RL_AGENT_C_H

                 
////////////////////////////////////////////////////////////////////////////////                
// File    : rl_agent_c.h                                                           
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


#ifdef __cplusplus
extern "C" {
#endif

/*
 * =============================================================================
 */

/* C-style declarations */

enum rl_act_type_t {
  RLA_PERM,
  RLA_GAUSS,
  RLA_BINARY,
  RLA_DISCRETE,
};

typedef void* rl_nac_t;

typedef struct rl_act_entry_t {
  rl_act_type_t type;
  int first_param; // sometimes the dimension; sometimes the number of samples
  int second_param; // only needed for RLA_DISCRETE
  void *vals;
} rl_act_entry_t;

typedef struct rl_act_desc_t {
  int act_cnt;
  rl_act_entry_t *acts;
} rl_act_desc_t;


rl_nac_t rl_nac_init( int num_state_feats, rl_act_desc_t *rad );
void rl_nac_deinit( rl_nac_t r );

void rl_nac_action_sample( rl_nac_t _r );
void rl_nac_update( rl_nac_t _r, double reward, double *statefeats );
void rl_nac_get_params( rl_nac_t _r, int act_num,
			int *param_cnt, double **params );

#ifdef __cplusplus
}
#endif

#endif
