
////////////////////////////////////////////////////////////////////////////////                
// File    : rl_agent.cpp                                                                          
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


#include "rl_agent.h"

//#define DEBUG

// after how many updates do we solve the linear system?
//#define UPDATE_THRESH 20
#define UPDATE_THRESH 200
//#define UPDATE_THRESH 100
//#define UPDATE_THRESH 20000
//#define UPDATE_THRESH 2000

// regularizer for when we solve the linear least squares
#define REGULARIZER 1e-6

// L1 regularizer for parameters -- drives everything towards zero in
// the absence of information
//#define PARAM_REGULARIZER 1e-4
#define PARAM_REGULARIZER 0

/*
 * =============================================================================
 */

rl_nac::rl_nac( int _state_feat_cnt ) {
  acts = new act_vec;

  state_feat_cnt = _state_feat_cnt;
  act_feat_cnt = 0;
  act_cnt = 0;

  gamma = 0.9;
  lambda = 0.1;
  beta = 0.0;  // for now

  updates_since_last_step = 0;
}

void rl_nac::add_action( rl_act_type_t type,
			 int first_param,
			 int second_param,
			 void *vals ) {
  rl_act *r;

  switch( type ) {

  case RLA_PERM:
    r = new rl_act_perm( first_param, vals );
    break;

  case RLA_GAUSS:
    r = new rl_act_gauss( first_param, vals );
    break;

  case RLA_BINARY:
    r = new rl_act_binary( first_param, vals );
    break;

  case RLA_DISCRETE:
    r = new rl_act_discrete( first_param, second_param, vals );
    break;
  }

  acts->push_back( r );
  act_feat_cnt += (*acts)[act_cnt]->params.GetM();;

  act_cnt++;

  update_sizes();

}

void rl_nac::update_sizes( void ) {
  // resize all arrays appropriately.

  int sz = state_feat_cnt + act_feat_cnt;

  A.Resize( sz, sz );
  b.Resize( sz );
  z.Resize( sz );

  phi_xt.Resize( state_feat_cnt );

  A.Fill( 0 );
  b.Fill( 0 );
  z.Fill( 0 );

  phi_xt.Fill( 0 );

}

/*
 * =============================================================================
 */

void rl_nac::update( double reward, double *statefeats ) {

  // incorporate new information
  add_obs( reward, statefeats );

  // we don't want to solve the linear system of equations every time.
  if ( time_to_check() ) {
    solve_for_grad();
    if ( check_for_convergence() ) {
      take_gradient_step();
    }
  }

  for ( int a=0; a<act_cnt; a++ ) {
    // drive parameters towards zero
    Add( -PARAM_REGULARIZER, (*acts)[a]->params, (*acts)[a]->params );
  }
}

void rl_nac::add_obs( double reward, double *statefeats ) {
  // the reward is r_t
  // the statefeats are phi(x_t+1)
  vector_t phi_xtp1;
  phi_xtp1.SetData( state_feat_cnt, statefeats );

#ifdef DEBUG
  printf( "phi_xtp1: " ); phi_xtp1.Print();
#endif

  // construct twid phi_t
  vector_t twid_phi_t, zeros;
  zeros.Resize( act_feat_cnt );
  zeros.Fill( 0 );
  vec_cat( phi_xtp1, zeros, twid_phi_t );
#ifdef DEBUG
  printf( "twid_phi_t: " ); twid_phi_t.Print();
#endif

  // construct hat phi_t 
  vector_t act_grad = collect_gradient();
#ifdef DEBUG
  printf( "Act_grad: " ); act_grad.Print();
#endif
  vector_t hat_phi_t;
  vec_cat( phi_xt, act_grad, hat_phi_t );
#ifdef DEBUG
  printf( "phi_xt: " ); phi_xt.Print();
#endif

  // update z_t
  Add( lambda, z, hat_phi_t );  // don't reverse these parameters...
  z = hat_phi_t;
#ifdef DEBUG
  printf( "z: " ); z.Print();
#endif

  // update A_t;  tmpv = hat_phi_t - gamma * twid_phi_t;
  vector_t tmpv = hat_phi_t;
  Add( -gamma, twid_phi_t, tmpv );
#ifdef DEBUG
  printf( "tmpv: " ); tmpv.Print();
#endif

  Rank1Update( 1.0, z, tmpv, A );
#ifdef DEBUG
  printf( "A: " ); A.Print();
#endif

  // update b
  Add( reward, z, b );
#ifdef DEBUG
  printf( "b: " ); b.Print();
#endif

  // rotate features
  phi_xt = phi_xtp1;

  // clean up
  phi_xtp1.Nullify();

  updates_since_last_step ++;

#ifdef DEBUG
  printf("\n\n");
#endif

}

vector_t rl_nac::collect_gradient( void ) {

  vector_t result( act_feat_cnt );

  int i = 0;
  for ( int a=0; a<act_cnt; a++ ) {
    int act_a_cnt = (*acts)[a]->grad.GetM();
    for ( int j=0; j<act_a_cnt; j++ ) {
      result(i) = (*acts)[a]->grad(j);
      i++;
    }   
  }

  return result;
}

void rl_nac::solve_for_grad( void ) {

  matrix_t copyA;
  vector_t copyb;

  LapackInfo info(0);

  copyA = A;
  copyb = b;
  
  // Regularize this: A = A + REGULARIZER*eye(N)
  int cnt = A.GetM();
  for ( int i=0; i<cnt; i++ ) {
    copyA(i,i) = copyA(i,i) + REGULARIZER;
  }

#ifdef DEBUG
  printf("=====================================\n");
  printf("Solving for gradient:\n");
  copyA.Print();
  copyb.Print();
  printf("=====================================\n");
#endif

  vector_t tau;
  GetQR( copyA, tau, info );

  if ( info.GetInfo() != 0 ) {
    fprintf( stderr, "Error in QR decomposition!\n" );
    for ( int a=0; a<act_cnt; a++ ) {
      (*acts)[a]->natural_grad.Fill(0);
    }
    return;
  }

  SolveQR( copyA, tau, copyb, info );

  if ( info.GetInfo() != 0 ) {
    fprintf( stderr, "Error solving linear system!\n" );
    for ( int a=0; a<act_cnt; a++ ) {
      (*acts)[a]->natural_grad.Fill(0);
    }
    return;
  }

  // results are stored in copyb
#ifdef DEBUG
  printf("Final answer:\n");
  copyb.Print();
#endif

  // now we need to farm these out to the actions!  we start i at
  // state_feat_cnt because copyb contains [w_t+1 v_t+1], but we only
  // want v_t+1 to add to the parameter vectors.
  int i = state_feat_cnt;
  for ( int a=0; a<act_cnt; a++ ) {
    int act_a_cnt = (*acts)[a]->grad.GetM();

    for ( int j=0; j<act_a_cnt; j++ ) {
      (*acts)[a]->natural_grad(j) = copyb(i);
      i++;
    }

    // make all vectors have unit length
    double n2 = Norm2( (*acts)[a]->natural_grad );

    if ( n2 > 1e-10 ) {
      // it's pretty easy to have a zero norm natural gradient, if you
      // received zero total reward, so to avoid nans, we only
      // normalize if we have a nonzero vector...
      Mlt( 1.0/n2, (*acts)[a]->natural_grad );
    }

#ifdef DEBUG
    printf( "Act %d natural grad:\n", a );
    (*acts)[a]->natural_grad.Print();
#endif
  }

}

bool rl_nac::check_for_convergence( void ) {
  return true;
}

void rl_nac::take_gradient_step( void ) {

#ifdef DEBUG
  printf("=====================================\n");
  printf("TAKING A STEP!!!\n");
  printf("=====================================\n");
#endif

  for ( int a=0; a<act_cnt; a++ ) {
#ifdef DEBUG
    printf( "Action %d:\n", a );
    printf( "  Parameters before: " );
    (*acts)[a]->params.Print();
    printf( "  Natural grad: " );
    (*acts)[a]->natural_grad.Print();
#endif

    Add( (*acts)[a]->stepsize,
	 (*acts)[a]->natural_grad,
	 (*acts)[a]->params );

#ifdef DEBUG
    printf( "  Parameters after: " );
    (*acts)[a]->params.Print();
#endif
  }

#ifdef DEBUG
  printf("=====================================\n");
#endif

  // XXX here, we should decay things with Beta.
  A.Fill( 0 );
  z.Fill( 0 );
  b.Fill( 0 );

  updates_since_last_step = 0;
}

bool rl_nac::time_to_check( void ) {
  if ( updates_since_last_step > UPDATE_THRESH ) {
    return true;
  }
  return false;
}

/*
 * =============================================================================
 */

/*

  In every case, we try to use softmax-style probability distributions
  to help alleviate boundary conditions.

 */

// -----------------------------------------------------------------------------

void rl_act_perm::act_sample( void ) {
  //  printf( "perm sample %d dim\n", dim );

  log_probs = params;
  grad = 1;

  for ( int i=0; i<dim; i++ ) {
    ll_to_prob( log_probs, probs );
    //    probs.Print();
    int a = sample( probs, &rng_state );

    vals[i] = a;

    // computation of grad = grad + -1*probs
    Add( -1, probs, grad );

    log_probs(a) = MFMIN;
  }

  //  printf( "Gradient: " ); grad.Print();
}

// -----------------------------------------------------------------------------

void rl_act_gauss::act_sample( void ) {
  //  printf( "gauss: %d samples\n", num_samples );

  // in this case, the parameter is the mean of the gaussian.
  // XXX we could also learn the variance, eh?

  for ( int i=0; i<num_samples; i++ ) {
    vals[i] = sample_gaussian( &rng_state ) + params(0);
    grad( 0 ) += vals[i] - params(0);
  }

  //  printf( "Gauss gradient: " ); grad.Print();
}

// -----------------------------------------------------------------------------

void rl_act_binary::act_sample( void ) {
  double tmp, p, ep, gradval;
  //  printf( "binary: %d samples\n", num_samples );

  ep = exp( -params(0) );
  p = 1.0 / ( 1.0 + ep );
  gradval = ep / ( 1.0 + ep );

  //  printf( "p=%.2f; gradval=%.2f\n", p, gradval );

  grad(0) = 0;

  for ( int i=0; i<num_samples; i++ ) {
    drand48_r( &rng_state, &tmp );
    vals[i] = ( tmp < p );
    if ( vals[i] ) {
      grad(0) += gradval;
    } else {
      grad(0) += -gradval;
    }
  }

  //  printf( "Gradient: " ); grad.Print();

}

// -----------------------------------------------------------------------------

void rl_act_discrete::act_sample( void ) {
  //  printf( "discrete sample: %d samples w/%d options\n", num_samples, num_opts );

  grad = 0;

  ll_to_prob( params, probs );

  for ( int i=0; i<num_samples; i++ ) {
    int a = sample( probs, &rng_state );
    vals[i] = a;
    // computation of grad = grad - probs
    Add( -1, probs, grad );
    grad( a ) += 1;
  }

  //  printf( "Gradient: " ); grad.Print();

}


/*
 * =============================================================================
 */

/* This is the C-style interface */

rl_nac_t rl_nac_init( int num_state_feats, rl_act_desc_t *rad ) {

  rl_nac *r = new rl_nac( num_state_feats );

  for ( int a=0; a<rad->act_cnt; a++ ) {
    r->add_action( rad->acts[a].type,
		   rad->acts[a].first_param,
		   rad->acts[a].second_param,
		   rad->acts[a].vals );
  }

  return r;
}

void rl_nac_deinit( rl_nac_t r ) {
   rl_nac *ther = (rl_nac*) r;
   delete ther;
}

void rl_nac_action_sample( rl_nac_t _r ) {
  int a;

  rl_nac *r = (rl_nac *)_r;

  /* this samples and populates the appropriate values and
     gradients */
  for ( a=0; a<r->act_cnt; a++ ) {
    (*(r->acts))[a]->act_sample();
  }

}

void rl_nac_update( rl_nac_t _r, double reward, double *statefeats ) {
  rl_nac *r = (rl_nac *)_r;
  r->update( reward, statefeats );
}

void rl_nac_get_params( rl_nac_t _r, int act_num,
			int *param_cnt, double **params ) {
  rl_nac *r = (rl_nac *)_r;
  
  (*param_cnt) = (*(r->acts))[act_num]->params.GetM();
  (*params) = (*(r->acts))[act_num]->params.GetData();
}


/*
 * =============================================================================
 */

