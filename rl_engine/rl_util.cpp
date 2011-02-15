                 
////////////////////////////////////////////////////////////////////////////////                
// File    : rl_util.cpp                                                                       
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


// minimum floating point value                                                 
#define MFMIN -1e300

//
// ===========================================================================
//
//  General utilities
//
// ===========================================================================
//

double max( vector_t &v ) {
  int cnt = v.GetM();
  double mymax = MFMIN;
  for ( int i=0; i<cnt; i++ ) {
    if ( v(i) > mymax ) mymax = v(i);
  }
  return mymax;
}

double sum( vector_t &v ) {
  int cnt = v.GetM();
  double total = 0;
  for ( int i=0; i<cnt; i++ ) {
    total += v(i);
  }
  return total;
}

void add( vector_t &v, double a ) {
  int cnt = v.GetM();
  for ( int i=0; i<cnt; i++ ) {
    v(i) += a;
  }
}

void vec_cat( vector_t &a, vector_t &b, vector_t &c ) {
  int i;

  c.Resize( a.GetM() + b.GetM() );

  for ( i=0; i<a.GetM(); i++ ) {
    c(i) = a(i);
  }

  for ( int j=0; j<b.GetM(); j++ ) {
    c(i) = b(j);
    i++;
  }
}

//
//
// ===========================================================================
//
//  Statistics utilities
//
// ===========================================================================
//

int randint( int maxval, struct drand48_data *rng_state ) {
  int rval;
  double tmpbob;

  drand48_r( rng_state, &tmpbob );
  rval = floor( ((double)maxval+1.0) * tmpbob );

  if ( rval > maxval ) // could happen if RV=1.0
    rval = maxval;

  return rval;
}

// given a multinomial probability vector, sample an index
int sample( vector_t &probs, struct drand48_data *rng_state ) {
  double p, cumsum;

  drand48_r( rng_state, &p );

  cumsum = 0;
  int cnt = probs.GetM();
  for ( unsigned int i=0; i<cnt; i++ ) {
    cumsum += probs(i);
    if ( p <= cumsum ) {
      return i;
    }
  }

  return cnt-1;
}

double sample_gaussian( struct drand48_data *rng_state ) {
  double u, v;
  drand48_r( rng_state, &u );
  drand48_r( rng_state, &v );
  // Box-Muller transform
  return sqrt( -2 * log(u) ) * cos( 2*M_PI*v );
}

int sample_unused_uniformly( vector_t &log_probs,
			     struct drand48_data *rng_state ) {
  double p, cumsum, num_unused;

  drand48_r( rng_state, &p );

  int cnt = log_probs.GetM();

  num_unused = 0;
  for ( unsigned int i=0; i<cnt; i++ ) {
    if ( log_probs(i) > MFMIN ) {
      num_unused++;
    }
  }

  cumsum = 0;
  for ( unsigned int i=0; i<cnt; i++ ) {
    if ( log_probs(i) > MFMIN ) {
      cumsum += (1.0/num_unused);
      if ( p <= cumsum ) {
	return i;
      }
    }
  }

  fprintf( stderr, "shouldn't get here!\n" );
  exit(0);
}

// convert a vector of log probabilities to probabilities,                      
// while avoiding underflow.                                                    
void ll_to_prob( vector_t &log_probs, vector_t &probs ) {

  int cnt = log_probs.GetM();

  probs = log_probs;
  add( probs, - max( log_probs ) );

  for ( unsigned int i=0; i<cnt; i++ ) {
    probs(i) = exp( probs(i) );
  }

  Mlt( 1.0 / sum( probs ), probs );

}
