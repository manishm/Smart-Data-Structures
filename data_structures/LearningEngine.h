
#ifndef LEARNING_ENGINE_H
#define LEARNING_ENGINE_H

////////////////////////////////////////////////////////////////////////////////
// File    : LearningEngine.h                                                  
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

#include <assert.h>
#include <pthread.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <math.h>
#include "Monitor.h"
#include "rl_agent_c.h"
#include "cpp_framework.h"
#include "portable_defns.h"

#include <unistd.h>

using namespace std;


static void* learningengine(void *);

class LearningEngine
{
public:

        enum learning_mode_t {
                disabled               = 0x0,
                random_lock_scheduling = 0x1,  //if set, use a random lock scheduling policy
                lock_scheduling        = 0x2,  //if set, use a learned lock scheduling policy
                scancount_tuning       = 0x4,  //if set, use a learned scancount setting
                manual_stepping        = 0x8,  //if set, no helper thread; require manual stepping
                inject_sleep           = 0x10, //if set, inject sleep each rl step
                inject_delay           = 0x20  //if set, inject delay time each rl step
        };

private:

        friend void* learningengine(void *);

        //general
        const learning_mode_t mode ATTRIBUTE_CACHE_ALIGNED;
        Monitor*              mon;        
        unsigned int          nthreads;

        //concurrency control
        static volatile unsigned int        htlock;
        static volatile unsigned int        refcount;
        static volatile unsigned int        signalquit;
        static volatile unsigned int        lelock;
        static volatile unsigned int        lelistpos;
        static std::vector<LearningEngine*> lelist;
        static volatile _u64                learnercount;
        volatile unsigned int               dellock  ATTRIBUTE_CACHE_ALIGNED;

        //threading
        pthread_attr_t   managerthreadattr  ATTRIBUTE_CACHE_ALIGNED;
        static pthread_t managerthread;

        //reward
        double   last_checkpointed_reward   ATTRIBUTE_CACHE_ALIGNED;
        double   total_reward_ever;
        timespec last_update_timestamp;

        //learning
        int                  num_lock_sched ATTRIBUTE_CACHE_ALIGNED;
        int                  num_sc_tune;
        volatile int         lock_sched_id;
        volatile int         sc_tune_id;
        int*                 perm_vals;
        int*                 disc_vals;
        int*                 ext_disc_vals;
        int*                 ext_perm_vals;
        double*              probs;
        rl_nac_t             r;
        struct drand48_data  rng_state      ATTRIBUTE_CACHE_ALIGNED;

        //slowdown 
        double               rl_to_sleepidle_ratio  ATTRIBUTE_CACHE_ALIGNED;
        _u64                 time_window[8];
        int                  time_indx;
        _u64                 injection_nanos;
        timespec             time1;
        timespec             time2;
        _u64                 time_sum;

        //padding
        char pad  ATTRIBUTE_CACHE_ALIGNED;


        //------------------------------------
        //Concurrency control and threading
        //------------------------------------

        static void lelistAdd(LearningEngine *s)
        {
                while( (0 != lelock) || !CAS(&lelock, 0, 1) );      
                lelist.push_back(s);
                CCP::Memory::read_write_barrier();
                lelock = 0;
        }

        static void lelistRemove(LearningEngine *s)
        {
                while( (0 != lelock) || !CAS(&lelock, 0, 1) );      
                std::vector<LearningEngine*>::iterator it = lelist.begin();
                while(*it != s)
                        it++;
                lelist.erase(it);
                CCP::Memory::read_write_barrier();
                lelock = 0;
        }
  
        static LearningEngine* lelistNext()
        {
               while( (0 != lelock) || !CAS(&lelock, 0, 1) );     
               if ( lelist.size() == 0 )
               {
                       CCP::Memory::read_write_barrier();
                       lelock = 0;
                       return NULL;
               }

               if ( lelistpos >= lelist.size() )
                       lelistpos = 0;

               LearningEngine *sl = lelist[lelistpos];
               lelistpos++;

               CCP::Memory::read_write_barrier();
               lelock = 0; 

               return sl;
        }

        void registerLearner(LearningEngine *l)
        {
               if (mode == disabled)
                       return; 

               //cerr << "In registerLearner with mode=" << mode << endl;
               
               initSlowdownMode();
               initAPI();
               initrl();

               while( (0 != htlock) || !CAS(&htlock, 0, 1) );

               learnercount += U64(0x0001000000000001);

               //CCP::Memory::read_write_barrier();

               if ( refcount++ == 0 )
               {
                       if ( 0 == (mode & manual_stepping) ) {
                               pthread_attr_init(&managerthreadattr);
                               pthread_attr_setdetachstate(&managerthreadattr, PTHREAD_CREATE_JOINABLE);
                               pthread_create(&managerthread, &managerthreadattr, learningengine, NULL);
                       }
               }

               CCP::Memory::read_write_barrier();
               htlock = 0;

               lelistAdd(l);
        }

        void unregisterLearner(LearningEngine *l)
        {
               if (mode == disabled)
                       return;

               //cerr << "In unregisterLearner with mode=" << mode << endl;

               lelistRemove(l);

               while( (0 != htlock) || !CAS(&htlock, 0, 1) );

               learnercount += U64(0xFFFF000000000001);     

               if ( --refcount == 0 )
               {
                       if ( 0 == (mode & manual_stepping) ) {
                              FAADD(&signalquit, 1);
                              pthread_join(managerthread, NULL);
                              FAADD(&signalquit, -1);
                       }
               }

               CCP::Memory::read_write_barrier();
               htlock = 0;

               while( (0 != dellock) || !CAS(&dellock, 0, 1) );

               deinitrl();
               deinitAPI();
               deinitSlowdownMode();

               CCP::Memory::read_write_barrier();
               dellock = 0;
        }


        //------------------------------------
        //Utilities
        //------------------------------------
        inline int clippriority(int pri)
        {
                //we have 64 discrete priority levels, 0-63, w/ 63 rsvd
                assert( pri >= 0 );
                return pri < 62 ? pri : 62;
        }

        //------------------------------------
        //Learning engine
        //------------------------------------

        void modeCheck() 
        {
                bool _is_disabled               = (mode == disabled);
                bool _is_random_lock_scheduling = (mode & random_lock_scheduling);
                bool _is_lock_scheduling        = (mode & lock_scheduling);
                bool _is_scancount_tuning       = (mode & scancount_tuning);
                bool _is_manual_stepping        = (mode & manual_stepping);
                bool _is_inject_sleep           = (mode & inject_sleep);
                bool _is_inject_delay           = (mode & inject_delay);

                bool err = false;
                err |= ( _is_random_lock_scheduling && (_is_lock_scheduling ||_is_scancount_tuning || _is_inject_sleep || _is_inject_delay ) );
                err |= ( _is_manual_stepping && ( _is_inject_sleep || _is_inject_delay ) );
                err |= ( _is_inject_sleep && _is_inject_delay );
                err |= ( _is_manual_stepping && !( _is_lock_scheduling || _is_scancount_tuning || _is_random_lock_scheduling) );
                err |= ( ( rl_to_sleepidle_ratio < .0001 ) || ( rl_to_sleepidle_ratio > 1.0 ) );
                err |= ( _is_lock_scheduling && (num_lock_sched < 1) );
                err |= ( _is_scancount_tuning && (num_sc_tune < 1) );

                if ( err ) {
                        cerr << "Sorry, unsupported mode: " << mode << endl;
                        exit(0);
                }
        }


        void initSlowdownMode()
        {
                const _u64 rlnanos = 1000;

                time_sum           = 8 * rlnanos;
                injection_nanos    = (_u64) ((((double) time_sum) / 8.) / rl_to_sleepidle_ratio);
                time_indx          = 0;

                for(int i = 0; i < 8; i++)
                        time_window[i] = rlnanos;
        }

        void deinitSlowdownMode()
        { }

        void initAPI()
        {
                //FIXME: the other cases aren't handled yet
                assert( (num_lock_sched == 0) || (num_lock_sched == 1) );
#if 1
                int permbytes = ((sizeof(int)*nthreads + CACHE_LINE_SIZE - 1) / CACHE_LINE_SIZE) * CACHE_LINE_SIZE;
                int discbytes = ((sizeof(int)*num_sc_tune + CACHE_LINE_SIZE - 1) / CACHE_LINE_SIZE) * CACHE_LINE_SIZE;

                perm_vals     = (int*) CCP::Memory::byte_aligned_malloc(permbytes, CACHE_LINE_SIZE);
                ext_perm_vals = (int*) CCP::Memory::byte_aligned_malloc(permbytes, CACHE_LINE_SIZE);
#else
                perm_vals     = new int[nthreads];
                ext_perm_vals = new int[nthreads];
#endif

                for (int i = 0; i<nthreads; ++i) {
                        int p = clippriority( nthreads - 1 - i );
                        perm_vals[i] = p;
                        ext_perm_vals[i] = p;
                }

                assert( num_sc_tune >= 0 );
#if 1
                disc_vals     = (int*) CCP::Memory::byte_aligned_malloc(discbytes, CACHE_LINE_SIZE);
                ext_disc_vals = (int*) CCP::Memory::byte_aligned_malloc(discbytes, CACHE_LINE_SIZE);
#else
                disc_vals     = new int[num_sc_tune];
                ext_disc_vals = new int[num_sc_tune];
#endif
                for (int i = 0; i<num_sc_tune; ++i) {
                        disc_vals[i] = nthreads;
                        ext_disc_vals[i] = nthreads;
                }
        }

        void deinitAPI()
        {
#if 1
	        CCP::Memory::byte_aligned_free(disc_vals);
		CCP::Memory::byte_aligned_free(ext_disc_vals);
		CCP::Memory::byte_aligned_free(perm_vals);
		CCP::Memory::byte_aligned_free(ext_perm_vals);
#else
                delete[] disc_vals;
                delete[] ext_disc_vals;
                delete[] perm_vals;
                delete[] ext_perm_vals;
#endif
        }

        void initrl()
        {
                if ( mode == disabled )
                        return;

                probs = new double[nthreads];
                srand48_r( 42, &rng_state );

                //cerr << "In initrl with mode=" << mode << endl;

                if ( (mode & lock_scheduling) && (mode & scancount_tuning) )
                {
                        //cerr << "Configuring ml for both" << endl;
                        //nthreads perm vals, 1 discrete (range 0 to 12)
                        rl_act_entry_t raes[] = 
                          {{ RLA_PERM, 0, 0, perm_vals },
                           { RLA_DISCRETE, 0, 13, disc_vals }};
                        raes[0].first_param = nthreads;
                        raes[1].first_param = num_sc_tune;
                        rl_act_desc_t rad = { 2, raes };

                        r = rl_nac_init( 1, &rad );
                }     
                else if ( mode & scancount_tuning ) {
                        //cerr << "Configuring ml for external discrete" << endl;
                        //0 perm vals, 1 discrete (range 0 to 12)
                        rl_act_entry_t raes[] = 
                          {{ RLA_DISCRETE, 0, 13, disc_vals }};
                        raes[0].first_param = num_sc_tune;
                        rl_act_desc_t rad = { 1, raes };

                        r = rl_nac_init( 1, &rad );
                } 
                else if ( mode & lock_scheduling ) {
                        //cerr << "Configuring ml for lock scheduling" << endl;
                        //nthreads perm vals
                        rl_act_entry_t raes[] = 
                          {{ RLA_PERM, 0, 0, perm_vals }};
                        raes[0].first_param = nthreads;
                        rl_act_desc_t rad = { 1, raes };

                        r = rl_nac_init( 1, &rad );
                }

                // initialize reward computation stuff
                last_checkpointed_reward = 0;
                total_reward_ever = 0;
                last_update_timestamp = rlgettime();
        }

        void deinitrl()
        {
                if ( mode == disabled )
                        return;

                //cerr << "In deinitrl mode=" << mode << endl;

                delete[] probs;

                if ( mode & (lock_scheduling | scancount_tuning) )
                        rl_nac_deinit( r );
        }
   
        timespec rlgettime()
        {
                timespec ts;
                clock_gettime( CLOCK_REALTIME, &ts );
                return ts;
        }

        double getmonitorsignal(_u64 lastval)
        {
	        //return mon ? mon->getrewardnotsafe() : 0;
                return mon ? mon->waitrewardnotsafe(lastval) : 0;
        }

        bool getreward( double *reward ) 
        {
                double total_reward_ever = getmonitorsignal(last_checkpointed_reward);
                double accumulated_heartbeats = total_reward_ever - last_checkpointed_reward;

                // all of the information we need to compute a rate.
                timespec tmp_time = rlgettime();
                _u64 elapsed = CCP::Thread::difftimes(last_update_timestamp, tmp_time);
                double time_elapsed = ((double) elapsed) / ((double) 1e9);

                // originals 10, .0001
                if ( (accumulated_heartbeats > 10) || (time_elapsed > 0.0001) ) {
                        last_update_timestamp = tmp_time;
                        double r = accumulated_heartbeats / time_elapsed;
                        *reward = r;
                        //cerr << "reward= " << r << endl;
                        last_checkpointed_reward = total_reward_ever;
                        return true;
                }
                
                return false;
        }

        void slowdown_part1()
        {
                //cerr << "injection_nanos= " << injection_nanos << endl;
                if ( mode & inject_sleep )
                        CCP::Thread::sleep(0, injection_nanos);
                else if ( mode & inject_delay ) {
                        CCP::Thread::delay(injection_nanos);
                }
                time1 = rlgettime();
        }

        void slowdown_part2()
        {
                time2 = rlgettime();
                _u64 delta = CCP::Thread::difftimes(time1,time2);
                time_sum = time_sum - time_window[time_indx] + delta;
                time_window[time_indx] = delta;
                time_indx = (time_indx + 1) % 8;
                double avg = ((double) time_sum) / 8.;
                injection_nanos = (_u64) (avg / rl_to_sleepidle_ratio);
        }

        void rlupdate() 
        {
                //cerr << "running rl update" << endl;
                //ML policy

	        if ( !CAS(&dellock, 0, 1) )
		        return;

                _u64 start_seq_no = (learnercount & U64(0xFFFFFFFFFFFF));
                _u64 status;
                _u64 num_live;
                _u64 seq_no;

                do
                {
                        if ( mode & (inject_sleep | inject_delay) )
                                slowdown_part1();

                        double reward;
                        if ( getreward( &reward ) )
                        {
                                rl_nac_action_sample( r );

                                //output the vals
                                memcpy(ext_disc_vals, disc_vals, num_sc_tune*sizeof(int));
                                for(int i = 0; i < nthreads; ++i)
                                        ext_perm_vals[i] = clippriority(nthreads - 1 - perm_vals[i]);
                                //cerr << "scancount = " << disc_vals[0] << endl;

                                double statefeats[] = {1.};
                                rl_nac_update( r, reward, statefeats );
                        }

                        if ( mode & (inject_sleep | inject_delay) )
                                slowdown_part2();

                        //cerr << "looping in rlupdate." << endl;

                        status = learnercount;
                        num_live = (status >> 48);
                        seq_no = (status & U64(0xFFFFFFFFFFFF));

                } while( (0 == (mode & manual_stepping)) && (num_live == 1) && (seq_no == start_seq_no) );

                CCP::Memory::read_write_barrier();
                dellock = 0;                
        }

        // given a multinomial probability vector, sample an index
        int sample( double *probs, struct drand48_data *rng_state ) {
                double p, cumsum;

                drand48_r( rng_state, &p );

                cumsum = 0;
                for (int i = 0; i < nthreads; i++ ) {
                        cumsum += probs[i];
                        if ( p <= cumsum )
                                return i;
                }

                return nthreads - 1;
        }

        void randupdate() 
        {
                //cerr << "calling rand update" << endl;

                // random policy

	        if ( !CAS(&dellock, 0, 1) )
		        return;

                // initialize them all equally likely (adjusted for initial renormalizing)
                for(int i = 0; i < nthreads; i++)
                        probs[i] = 1.;
                double sum = nthreads;

                for ( int priority=nthreads-1; priority>=0; priority-- ) {
                  
                        // renormalize
                        for(int i = 0; i < nthreads; i++)
                                probs[i] /= sum;

                        // sample and set priority
                        int a = sample( probs, &rng_state );
                        ext_perm_vals[a] = clippriority( priority );

                        // clear out the chosen one; set next sum
                        probs[a] = 0.;
                        sum = ((double) priority) / (priority + 1);        
                }

                CCP::Memory::read_write_barrier();
                dellock = 0; 
        }



 public:

        //---------------------
        //Interface
        //---------------------

         LearningEngine(unsigned int threads, Monitor *m, double rlfactor = 1.0, 
                        learning_mode_t mode = disabled, int num_lock_scheduling = 0, int num_scancount_tuning = 0 )
        :  nthreads(threads), 
           mon(m), 
           mode(mode),
           rl_to_sleepidle_ratio(rlfactor),
           num_lock_sched(num_lock_scheduling),
           num_sc_tune(num_scancount_tuning),
           lock_sched_id(0),
	   sc_tune_id(0),
	   dellock(0)
        {
	        //cerr << "instantiated a learner. num_sc_tune= " << num_sc_tune << endl;
                modeCheck();
                registerLearner(this);
        }

        ~LearningEngine() 
        {
                unregisterLearner(this);
        }

        inline int register_lock_sched_id()
	{ 
	        //FIXME: other cases aren't handled yet
	        assert( num_lock_sched == 1 ); 
                return 0;
	}

        inline int register_sc_tune_id()
	{
                int id = FAADD(&sc_tune_id, 1);
		return id;
	}

        inline int getdiscval(unsigned int sc_tune_id, unsigned int tid)
        {
                return ext_disc_vals[sc_tune_id];
        }

        inline int getpermval(unsigned int lock_sched_id, unsigned int tid)
        {
	        //FIXME. this will need to change when we support > 1 perm obj
                return ext_perm_vals[tid];
        }

        inline learning_mode_t getmode()
        {
                return mode;
        }

};



//------------------------------
//Learning Helper Thread
//------------------------------

static void* learningengine(void *)
{
        // get this valid, up-to-date in our cache
        LearningEngine::signalquit = FAADD(&LearningEngine::signalquit, 0);

        while( 0 == LearningEngine::signalquit )
        {
                // get next live lock (or retry)
                LearningEngine *le = LearningEngine::lelistNext();
                if ( le == NULL )
                        continue;
      
                if ( le->mode == LearningEngine::random_lock_scheduling )
                        le->randupdate();
                else
                        le->rlupdate();
        }

        return NULL;
}

#endif

