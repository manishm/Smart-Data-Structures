
#ifndef RL_AGENT_H
#define RL_AGENT_H

                 
////////////////////////////////////////////////////////////////////////////////                
// File    : rl_agent.h                                                                 
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


/*

We want to be able to parametrically describe a set of actions.  The
idea is that we'll have a set of objects that we know how to sample
from, and that we can compute derivatives from.

Supported types right now are:

* Uniform over k elements
* Binary
* Permutation
* Multivariate (symmetrical) gaussian

 */

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// can only include this once...
#define SELDON_WITH_CBLAS
#define SELDON_WITH_LAPACK
#include "Seldon.hxx"
using namespace Seldon;

typedef Matrix<double> matrix_t;
typedef Vector<double> vector_t;

#include "rl_agent_c.h"

// =============================================================================

#include "rl_util.cpp"

// =============================================================================

unsigned int _unique_seed = 0;

class rl_act {
 public:
  rl_act() {
    srand48_r( time(NULL) + _unique_seed, &rng_state );
    _unique_seed++;
  };

  virtual void act_sample( void ) {};

  // thread-safe random number generator state
  struct drand48_data rng_state;

  // everything has a gradient vector.
  vector_t grad;
  // this is what we'll eventually use to take a gradient step.
  vector_t natural_grad;

  // everything has a real-valued parameter vector.
  vector_t params;

  // each action type has a different stepsize (!)
  double stepsize;
};

// =============================================================================

class rl_act_perm : public rl_act {
 public:
  rl_act_perm( int _dim, void *_vals ) {
    init( _dim, (unsigned int *)_vals );
  };
  rl_act_perm( int _dim, unsigned int *_vals ) {
    init( _dim, _vals );
  };
  void init( int _dim, unsigned int *_vals ) {
    dim = _dim;
    vals = _vals;

    stepsize = 2.0;
    //stepsize = 10.0;

    grad.Resize( _dim );
    natural_grad.Resize( _dim );
    params.Resize( _dim );
    log_probs.Resize( _dim );
    probs.Resize( _dim );
  };
  void act_sample( void );

  int dim;
  unsigned int *vals;
  vector_t log_probs, probs;
};

class rl_act_gauss : public rl_act {
 public:
  rl_act_gauss( int _num_samples, void *_vals ) {
    init( _num_samples, (double *)_vals );
  };
  rl_act_gauss( int _num_samples, double *_vals ) {
    init( _num_samples, _vals );
  };
  void init( int _num_samples, double *_vals ) {
    num_samples = _num_samples;
    vals = _vals;

    stepsize = 0.1;

    grad.Resize( 1 );
    natural_grad.Resize( 1 );
    params.Resize( 1 );
  };
  void act_sample( void );

  int num_samples;
  double *vals;
};

class rl_act_binary : public rl_act {
 public:
  rl_act_binary( int _num_samples, void *_vals ) {
    init( _num_samples, (unsigned char *)_vals );
  };
  rl_act_binary( int _num_samples, unsigned char *_vals ) {
    init( _num_samples, _vals );
  };
  void init( int _num_samples, unsigned char *_vals ) {
    num_samples = _num_samples;
    vals = _vals;

    stepsize = 2.0;

    grad.Resize( 1 );
    natural_grad.Resize( 1 );
    params.Resize( 1 );
  };
  void act_sample( void );

  int num_samples;
  unsigned char *vals;
};

class rl_act_discrete : public rl_act {
 public:
  rl_act_discrete( int _num_samples, int _num_opts, void *_vals ) {
    init( _num_samples, _num_opts, (unsigned int *)_vals );
  };
  rl_act_discrete( int _num_samples, int _num_opts, unsigned int *_vals ) {
    init( _num_samples, _num_opts, _vals );
  };
  void init( int _num_samples, int _num_opts, unsigned int *_vals ) {
    num_samples = _num_samples;
    num_opts = _num_opts;
    vals = _vals;

    stepsize = .1;
    //stepsize = .5;
    //stepsize = 6.0;

    grad.Resize( num_opts );
    natural_grad.Resize( num_opts );
    params.Resize( num_opts );
    probs.Resize( num_opts );
  };
  void act_sample( void );

  int num_samples;    // number of things to sample
  int num_opts;       // number of options 
  unsigned int *vals; // place to write the specific values sampled
  vector_t probs;
};

// =============================================================================

typedef std::vector< rl_act * > act_vec;

class rl_nac {

 public:

  // construction and setup
  rl_nac( int state_feat_cnt );
  void add_action( rl_act_type_t type,
		   int first_param,
		   int second_param,
		   void *vals );
  void update_sizes( void );

  // NAC algorithm methods
  void update( double reward, double *statefeats );
  void add_obs( double reward, double *statefeats );
  vector_t collect_gradient( void );
  void solve_for_grad( void );
  bool check_for_convergence( void );
  void take_gradient_step( void );
  bool time_to_check( void );
  
  /* description of state features */
  int state_feat_cnt;

  /* action space descriptors */  
  int act_feat_cnt;
  int act_cnt;
  act_vec *acts;

  /* nac parameters */
  double gamma;  /* rl discount factor */
  double beta;   /* forgetting factor */
  double lambda; /* eligibility trace parameter */
  
  /* current state of all needed variables */
  matrix_t A;
  vector_t b, z;
  vector_t phi_xt;

  // some hackish stuff
  int updates_since_last_step;
};

// =============================================================================

#endif
