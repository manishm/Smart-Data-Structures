
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
                disabled = 0x0,
                random_lock_scheduling = 0x1,
                lock_scheduling = 0x2,
                scancount_tuning = 0x4
        };

private:

        const learning_mode_t mode ATTRIBUTE_CACHE_ALIGNED;
        
        friend void* learningengine(void *);

        //general
        unsigned int                        nthreads ATTRIBUTE_CACHE_ALIGNED;

        //concurrency control
        static volatile unsigned int        htlock        ATTRIBUTE_CACHE_ALIGNED;
        static volatile unsigned int        refcount      ATTRIBUTE_CACHE_ALIGNED;
        static volatile unsigned int        signalquit    ATTRIBUTE_CACHE_ALIGNED; //aligned?
        static volatile unsigned int        lelock        ATTRIBUTE_CACHE_ALIGNED;
        static volatile unsigned int        lelistpos     ATTRIBUTE_CACHE_ALIGNED;
        static std::vector<LearningEngine*> lelist        ATTRIBUTE_CACHE_ALIGNED; //aligned?
        static volatile _u64                learnercount  ATTRIBUTE_CACHE_ALIGNED;

        //threading
        pthread_attr_t   managerthreadattr  ATTRIBUTE_CACHE_ALIGNED;
        static pthread_t managerthread      ATTRIBUTE_CACHE_ALIGNED; //aligned?

        //reward
        double   last_checkpointed_reward   ATTRIBUTE_CACHE_ALIGNED;
        double   total_reward_ever          ATTRIBUTE_CACHE_ALIGNED; //aligned?
        double   last_update_timestamp      ATTRIBUTE_CACHE_ALIGNED; //aligned?
        Monitor  *mon                       ATTRIBUTE_CACHE_ALIGNED; //aligned?

        //learning
        rl_nac_t             r              ATTRIBUTE_CACHE_ALIGNED;
        int*                 perm_vals      ATTRIBUTE_CACHE_ALIGNED; //aligned?
        int*                 disc_vals      ATTRIBUTE_CACHE_ALIGNED; //aligned?
        int*                 ext_disc_vals  ATTRIBUTE_CACHE_ALIGNED; //aligned?
        int*                 ext_perm_vals  ATTRIBUTE_CACHE_ALIGNED; //aligned?
        double*              probs          ATTRIBUTE_CACHE_ALIGNED; //aligned?
        struct drand48_data  rng_state      ATTRIBUTE_CACHE_ALIGNED; //aligned?

        //padding
        char pad[CACHE_LINE_SIZE];


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

        void registerLearner(learning_mode_t mode, LearningEngine *l)
        {
               if (mode == disabled)
                       return; 

               //cerr << "In registerLearner with mode=" << mode << endl;
               
               while( (0 != htlock) || !CAS(&htlock, 0, 1) );

               learnercount += U64(0x0001000000000001);

               initrl(mode);
               //CCP::Memory::read_write_barrier();

               if ( refcount++ == 0 )
               {
                       pthread_attr_init(&managerthreadattr);
                       pthread_attr_setdetachstate(&managerthreadattr, PTHREAD_CREATE_JOINABLE);
                       pthread_create(&managerthread, &managerthreadattr, learningengine, NULL);
               }
               CCP::Memory::read_write_barrier();
               htlock = 0;

               lelistAdd(l);
        }

        void unregisterLearner(learning_mode_t mode, LearningEngine *l)
        {
               if (mode == disabled)
                       return;

               //cerr << "In unregisterLearner with mode=" << mode << endl;

               lelistRemove(l);

               while( (0 != htlock) || !CAS(&htlock, 0, 1) );

               learnercount += U64(0xFFFF000000000001);      

               if ( --refcount == 0 )
               {
                       FAADD(&signalquit, 1);
                       pthread_join(managerthread, NULL);
                       FAADD(&signalquit, -1);
                       //deinitrl();
               }

               deinitrl(mode);
               CCP::Memory::read_write_barrier();
               htlock = 0;
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

        void initAPI()
        {
                ext_disc_vals = new int[nthreads];
                ext_perm_vals = new int[nthreads];
                for (int i = 0; i<nthreads; ++i) {
                        ext_disc_vals[i] = nthreads;
                        ext_perm_vals[i] = clippriority( nthreads - 1 - i );
                }
        }

        void deinitAPI()
        {
                delete[] ext_disc_vals;
                delete[] ext_perm_vals;
        }

        void initrl(learning_mode_t mode)
        {
                if ( mode == disabled )
                        return;

                //cerr << "In initrl with mode=" << mode << endl;

                probs = new double[nthreads];
                srand48_r( 42, &rng_state );
                perm_vals = new int[nthreads];
                disc_vals = new int[nthreads];

                if ( (mode & lock_scheduling) && (mode & scancount_tuning) )
                {
                        //cerr << "Configuring ml for both" << endl;
                        //nthreads perm vals, 1 discrete (range 0 to 12)
                        rl_act_entry_t raes[] = 
                          {{ RLA_PERM, 0, 0, perm_vals },
                           { RLA_DISCRETE, 1, 13, disc_vals }};
                        raes[0].first_param = nthreads;
                        rl_act_desc_t rad = { 2, raes };

                        r = rl_nac_init( 1, &rad );
                }     
                else if ( mode & scancount_tuning ) {
                        //cerr << "Configuring ml for external discrete" << endl;
                        //0 perm vals, 1 discrete (range 0 to 12)
                        rl_act_entry_t raes[] = 
                          {{ RLA_DISCRETE, 1, 13, disc_vals }};
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

        void deinitrl(learning_mode_t mode)
        {
                if ( mode == disabled )
                        return;

                //cerr << "In deinitrl mode=" << mode << endl;

                delete[] probs;
                delete[] perm_vals;
                delete[] disc_vals;

                if ( mode & (lock_scheduling | scancount_tuning) )
                        rl_nac_deinit( r );
        }
   
        double rlgettime()
        {
                timespec ts;
                clock_gettime( CLOCK_REALTIME, &ts );
                return ((double)ts.tv_sec) + (((double)ts.tv_nsec) / ((double)1e9));
        }

        double getmonitorsignal()
        {
                return mon ? mon->getreward() : 0;  
        }

        int getreward( double *reward ) 
        {
                double total_reward_ever = getmonitorsignal();
                double accumulated_heartbeats = total_reward_ever - last_checkpointed_reward;

                // all of the information we need to compute a rate.
                double tmp_time = rlgettime();
                double time_elapsed = tmp_time - last_update_timestamp;

                // originals 10, .0001
                if ( accumulated_heartbeats > 10 || time_elapsed > 0.0001 ) {
                        last_update_timestamp = tmp_time;
                        *reward = accumulated_heartbeats / time_elapsed;
                        //cout << "using " << accumulated_heartbeats << " beats" << endl;
                        last_checkpointed_reward = total_reward_ever;
                        return 1;
                }
                
                return 0;
        }

        void rlupdate() 
        {
                //cerr << "running rl update" << endl;
                //ML policy

                _u64 start_seq_no = (learnercount & U64(0xFFFFFFFFFFFF));
                _u64 status;
                _u64 num_live;
                _u64 seq_no;

                do
                {
                        double reward;
                        if ( !getreward( &reward ) )
                        {
                                status = learnercount;
                                num_live = (status >> 48);
                                seq_no = (status & U64(0xFFFFFFFFFFFF));
                        }
                        else
                        {
                                rl_nac_action_sample( r );

                                //output the vals
                                memcpy(ext_disc_vals, disc_vals, nthreads*sizeof(int));
                                for(int i = 0; i < nthreads; ++i)
                                        ext_perm_vals[i] = clippriority(nthreads - 1 - perm_vals[i]);

                                double statefeats[] = {1.};
                                rl_nac_update( r, reward, statefeats );

                                status = learnercount;
                                num_live = (status >> 48);
                                seq_no = (status & U64(0xFFFFFFFFFFFF));
                        }

                        //cerr << "looping in rlupdate." << endl;

                } while( (num_live == 1) && (seq_no == start_seq_no) );

                CCP::Memory::read_write_barrier();

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
        }



 public:

        //---------------------
        //Interface
        //---------------------

        LearningEngine(unsigned int threads, Monitor *m, learning_mode_t mode = disabled)
        :  nthreads(threads), 
           mon(m), 
           mode(mode)
        {
                initAPI();
                registerLearner(mode, this);
        }

        ~LearningEngine() 
        {
                unregisterLearner(mode, this);
                deinitAPI();
        }

        inline int getdiscval(unsigned int id = 0)
        {
                return ext_disc_vals[id];
        }

        inline int getpermval(unsigned int id)
        {
                return ext_perm_vals[id];
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

