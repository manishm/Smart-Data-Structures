////////////////////////////////////////////////////////////////////////////////
// File    : main.cpp
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com
//           Ms.Moran Tzafrir  email: morantza@gmail.com
// Written : 16 February 2011, 13 October 2009
// 
// Multi-Platform C++ framework example benchamrk program
//
// Copyright (C) 2011 Jonathan Eastep, 2009 Moran Tzafrir.
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
// TODO:
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//INCLUDE DIRECTIVES
////////////////////////////////////////////////////////////////////////////////
#include <string>
#include <iostream>

//general includes .....................................
#include "cpp_framework.h"
#include "Configuration.h"
#include "Heartbeat.h"
#include "LazyCounter.h"

//FC research includes .................................
//queues
#include "FCQueue.h"
#include "SmartQueue.h"
#include "MSQueue.h"
#include "BasketsQueue.h"
//#include "ComTreeQueue.h"
#include "OyamaQueue.h"
//#include "OyamaQueueCom.h"
//skiplists
#include "FCSkipList.h"
#include "SmartSkipList.h"
#include "LFSkipList.h"
#include "LazySkipList.h"
//pairheaps
#include "FCPairingHeap.h"
#include "SmartPairingHeap.h"
//stacks
#include "FCStack.h"
//#include "SmartStack.h"  //doesn't exist yet
#include "LFStack.h"
#include "EliminationStack.h"


////////////////////////////////////////////////////////////////////////////////
//CONSTS
////////////////////////////////////////////////////////////////////////////////
static final int _ACTIONS_ARY_SIZE = 2*1024*1024;

////////////////////////////////////////////////////////////////////////////////
//GLOBALS
////////////////////////////////////////////////////////////////////////////////
static int                              _num_ds;
static FCBase<FCIntPtr>*                _gDS[1024];
static Configuration                    _gConfiguration;

static Random                           _gRand;
static int                              _gActionAry[_ACTIONS_ARY_SIZE];
static int                              _gTotalRandNum;

static int                              _g_thread_fill_table_size;
static int*                             _gRandNumAry;
static tick_t volatile *                _gThreadResultAry;

static int                              _gNumProcessors;
       int                              _gNumThreads;
       int                              _gIsDedicatedMode;
static int                              _gThroughputTime;

static Thread**                         _gThreads;
static AtomicInteger                    _gThreadStartCounter(0);
static AtomicInteger                    _gThreadEndCounter(0);
static VolatileType<tick_t>             _gStartTime(U64(0));
static VolatileType<tick_t>             _gEndTime(U64(0));

char                                    _pad1[64];
volatile int                            _gIsStopThreads(U64(0));
char                                    _pad2[64];

static tick_t                           _gResult = 0L;
static tick_t                           _gResultAdd = 0L;
static tick_t                           _gResultRemove = 0L;
static tick_t                           _gResultPeek = 0L;

static _u64 volatile                    _seed;
static boolean                          _is_tm=false;
static boolean                          _is_view=false;

static Monitor*                         _mon;

static final int                        _num_work_amts             = 10;
static int                            _work_amts[_num_work_amts] = {800, 6400, 200, 3200, 1600, 100, 400, 100, 400, 800};
//static int                            _work_amts[_num_work_amts] = {1600, 200, 400, 1600, 200, 1600, 3200, 100, 200, 800};
//static int                              _work_amts[_num_work_amts] = {800, 100, 6400, 200, 200, 100, 400, 800, 3200, 400};
//static int                            _work_amts[_num_work_amts] = {400, 200, 400, 200, 400, 200, 400, 200, 400, 200};

////////////////////////////////////////////////////////////////////////////////
//FORWARD DECLARETIONS
////////////////////////////////////////////////////////////////////////////////
void RunBenchmark();
void PrepareActions();
void PrepareRandomNumbers(final int size);
int NearestPowerOfTwo(final int x);
FCBase<FCIntPtr>* CreateDataStructure(char* final alg_name, LearningEngine* learner);


////////////////////////////////////////////////////////////////////////////////
//TYPEDEFS
////////////////////////////////////////////////////////////////////////////////

//this is used for benchmarking smart data structures in the independent-computation mode
class MixThread : public Thread {
public:
        final int _threadNo;

        MixThread (final int inThreadNo) : _threadNo(inThreadNo) {}

        void run() {

                PtrNode<FCIntPtr>* n = new FCIntPtrNode(0);

                //fill table ...........................................................
                for (int iNum=0; iNum < (_g_thread_fill_table_size/16); ++iNum) {

                                for (int iDb=0; iDb<_num_ds; ++iDb) {
                                        for (int i=0; i<16; ++i) {
                                                *n = FCIntPtrNode( Random::getRandom(_seed, _gConfiguration._capacity) + 2 );
                                                _gDS[iDb]->add(_threadNo, n);
                                                //_gDS[iDb]->add(_threadNo, Random::getRandom(_seed, _gConfiguration._capacity) + 2);
                                        }
                                }

                }
                for (int iDb=0; iDb<_num_ds; ++iDb) {
                        _gDS[iDb]->cas_reset(_threadNo);
                }

                //save start benchmark time ............................................
                final int start_counter = _gThreadStartCounter.getAndIncrement();
                if(start_counter == (_gNumThreads-1))
                        _gStartTime = System::currentTimeMillis();
                while((_gNumThreads) != _gThreadStartCounter.get()) {int i=_gThreadStartCounter.get();}

                //start thread benchmark ...............................................
                tick_t action_counter = 0;
                int iNumAdd = start_counter*1024;
                int iNumRemove = start_counter*1024;
                int iNumContain = start_counter*1024;
                int iOp = start_counter*128;
                do {
                        final int op = _gActionAry[iOp];

                        if(1==op) {
                                for (int iDb=0; iDb<_num_ds; ++iDb) {
                                        *n = FCIntPtrNode( _gRandNumAry[iNumAdd] );
                                        final int enq_value = _gDS[iDb]->add(_threadNo, n);
                                        ++iNumAdd;
                                        if(iNumAdd >= _gTotalRandNum)
                                                iNumAdd=0;
                                }
                                ++action_counter;
                                if ( 0 == _gConfiguration._internal_reward_mode )
				        _mon->addreward(_threadNo, _num_ds);
                        } else if(2==op) {
                                for (int iDb=0; iDb<_num_ds; ++iDb) {
                                        *n = FCIntPtrNode( _gRandNumAry[iNumRemove] );
                                        PtrNode<FCIntPtr>* final deq_value = _gDS[iDb]->remove(_threadNo, n);
                                        ++iNumRemove;
                                        if(iNumRemove >= _gTotalRandNum) { iNumRemove=0; }
                                }
                                ++action_counter;
                                if ( 0 == _gConfiguration._internal_reward_mode )
				        _mon->addreward(_threadNo, _num_ds);
                        } else {
                                for (int iDb=0; iDb<_num_ds; ++iDb) {
                                        *n = FCIntPtrNode( _gRandNumAry[iNumContain] );
                                        _gDS[iDb]->contain(_threadNo, n);
                                        ++iNumContain;
                                        if(iNumContain >= _gTotalRandNum) {     iNumContain=0; }
                                }
                                ++action_counter;
                                if ( 0 == _gConfiguration._internal_reward_mode )
				        _mon->addreward(_threadNo, _num_ds);
                        }

                        ++iOp;
                        if (iOp >= _ACTIONS_ARY_SIZE) {
                                iOp=0;
                        }

                        if (_gConfiguration._read_write_delay > 0) {
                                _gDS[0]->post_computation(_threadNo);
                        }

                        if (0 != _gIsStopThreads) {
                                //cerr << "stopping thread " << _threadNo << endl;
                                break;
                        }

                } while(true);

                //save end benchmark time ..............................................
                final int end_counter = _gThreadEndCounter.getAndIncrement();
                if(end_counter == (_gNumThreads-1)) {
                        _gEndTime = System::currentTimeMillis();
                }

                //save thread benchmark result .........................................
                _gThreadResultAry[_threadNo] = action_counter;
        }
};

//this is used for benchmarking smart data structures in producer-consumer mode
class AddThread : public Thread {
public:
        final int _threadNo;

        AddThread (final int inThreadNo) : _threadNo(inThreadNo) {}

        void run() {

                PtrNode<FCIntPtr>* n = new FCIntPtrNode(0);

                int iDb = _threadNo % _num_ds;

                for (int iDb=0; iDb<_num_ds; ++iDb) {
                        _gDS[iDb]->cas_reset(_threadNo);
		}

                //save start benchmark time ............................................
                final int start_counter = _gThreadStartCounter.getAndIncrement();
                if(start_counter == (_gNumThreads-1))
                        _gStartTime = System::currentTimeMillis();
                while((_gNumThreads) != _gThreadStartCounter.get()) {int i=_gThreadStartCounter.get();}

                //start thread benchmark ...............................................
                tick_t action_counter = 0;
                int iNumAdd = start_counter*1024;
                int iNumRemove = start_counter*1024;
                int iNumContain = start_counter*1024;
                int iOp = start_counter*128;
                do {

		        //for (int iDb=0; iDb<_num_ds; ++iDb) {
                                *n = FCIntPtrNode( _gRandNumAry[iNumAdd] );
                                final int enq_value = _gDS[iDb]->add(_threadNo, n);
                                ++iNumAdd;
                                if(iNumAdd >= _gTotalRandNum)
                                iNumAdd=0;
			//}
                        ++action_counter;
                        if ( 0 == _gConfiguration._internal_reward_mode )
				_mon->addreward(_threadNo, 1);

                        /*
                        if (_gConfiguration._read_write_delay > 0) {
                                _gDS[0]->post_computation(_threadNo);
                        }
                        */

                        if (0 != _gIsStopThreads) {
                                //cerr << "stopping thread " << _threadNo << endl;
                                break;
                        }

                } while(true);

                //save end benchmark time ..............................................
                final int end_counter = _gThreadEndCounter.getAndIncrement();
                if(end_counter == (_gNumThreads-1)) {
                        _gEndTime = System::currentTimeMillis();
                }

                //save thread benchmark result .........................................
                _gThreadResultAry[_threadNo] = action_counter;
        }
};

//this is used for benchmarking smart data structures and monitors in producer-consumer mode
//for monitors, remove is mapped to decrementing reward
class RemoveThread : public Thread {
public:
        final int _threadNo;

        RemoveThread (final int inThreadNo) : _threadNo(inThreadNo) {}

        void run() {

	        int iDb = _threadNo % _num_ds;

                for (int iDb=0; iDb<_num_ds; ++iDb) {
                        _gDS[iDb]->cas_reset(_threadNo);
		}

                //save start benchmark time ............................................
                final int start_counter = _gThreadStartCounter.getAndIncrement();
                if(start_counter == (_gNumThreads-1))
                        _gStartTime = System::currentTimeMillis();
                while((_gNumThreads) != _gThreadStartCounter.get()) {int i=_gThreadStartCounter.get();}

                //start thread benchmark ...............................................
                tick_t action_counter = 0;
                int iNumAdd = start_counter*1024;
                int iNumRemove = start_counter*1024;
                int iNumContain = start_counter*1024;
                int iOp = start_counter*128;

                do {
		        //for (int iDb=0; iDb<_num_ds; ++iDb) {
                                while ( (0 == _gIsStopThreads) && (0 == _gDS[iDb]->remove(_threadNo, NULL)) );  
                                ++iNumRemove;
                                if(iNumRemove >= _gTotalRandNum) { iNumRemove=0; }
			//}
                        ++action_counter;
                        if ( 0 == _gConfiguration._internal_reward_mode )
				_mon->addreward(_threadNo, 1);

                        if (_gConfiguration._read_write_delay > 0) {
                                _gDS[0]->post_computation(_threadNo);
                        }

                        if (0 != _gIsStopThreads) {
                                //cerr << "stopping thread " << _threadNo << endl;
                                break;
                        }

                } while(true);

                //save end benchmark time ..............................................
                final int end_counter = _gThreadEndCounter.getAndIncrement();
                if(end_counter == (_gNumThreads-1)) {
                        _gEndTime = System::currentTimeMillis();
                }

                //save thread benchmark result .........................................
                _gThreadResultAry[_threadNo] = action_counter;
        }
};


//this is used for benchmarking monitors in producer-consumer mode. 
//contains is mapped to waitrewardnotsafe
class ContainThread : public Thread {
public:
        final int _threadNo;

        ContainThread (final int inThreadNo) : _threadNo(inThreadNo) {}

        void run() {

	        //PtrNode<FCIntPtr>* n = new FCIntPtrNode(0);

                int iDb = _threadNo % _num_ds;

                for (int iDb=0; iDb<_num_ds; ++iDb) {
                        _gDS[iDb]->cas_reset(_threadNo);
		}

                //save start benchmark time ............................................
                final int start_counter = _gThreadStartCounter.getAndIncrement();
                if(start_counter == (_gNumThreads-1))
                        _gStartTime = System::currentTimeMillis();
                while((_gNumThreads) != _gThreadStartCounter.get()) {int i=_gThreadStartCounter.get();}

                //start thread benchmark ...............................................
                tick_t action_counter = 0;
                int iNumAdd = start_counter*1024;
                int iNumRemove = start_counter*1024;
                int iNumContain = start_counter*1024;
                int iOp = start_counter*128;
                _u64 last = 0;

                do {
		        //for (int iDb=0; iDb<_num_ds; ++iDb) {
		                //*n = FCIntPtrNode( _gRandNumAry[iNumContain] );
		                last = (_u64) _gDS[iDb]->contain(_threadNo, (FCIntPtrNode*) last);
                                ++iNumContain;
                                if(iNumContain >= _gTotalRandNum) {     iNumContain=0; }
		        //}
                        ++action_counter;
                        if ( 0 == _gConfiguration._internal_reward_mode )
				_mon->addreward(_threadNo, 1);

                        if (_gConfiguration._read_write_delay > 0) {
                                _gDS[0]->post_computation(_threadNo);
                        }

                        if (0 != _gIsStopThreads) {
                                //cerr << "stopping thread " << _threadNo << endl;
                                break;
                        }

                } while(true);

                //save end benchmark time ..............................................
                final int end_counter = _gThreadEndCounter.getAndIncrement();
                if(end_counter == (_gNumThreads-1)) {
                        _gEndTime = System::currentTimeMillis();
                }

                //save thread benchmark result .........................................
                _gThreadResultAry[_threadNo] = action_counter;
        }
};

////////////////////////////////////////////////////////////////////////////////
//MAIN
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv) {

        //..........................................................................
        _seed = Random::getSeed();

        //contaminate memory manager ...............................................
        for (int i=0; i< (1024*64); ++i) {
                void* final rand_mem = malloc( _gRand.nextInt(128)+1 );
                free(rand_mem);
        }
        
        //read benchmark configuration .............................................
        if(!_gConfiguration.read(argc, argv)) {
                if(!_gConfiguration.read()) {
                        System_out_println("USAGE: <algorithm> <testNum> <numThreads> <numActions> <maxKey> <insertOps> <deleteOps> <unrandom> <loadFactor> <badHash> <initialCount> <_throughput_time>");
                        System::exit(-1);
                }
        }

        //initialize global variables ..............................................
        int tmp = _gConfiguration._alg1_num + _gConfiguration._alg2_num + _gConfiguration._alg3_num + _gConfiguration._alg4_num;
        bool concurrent = (tmp != 1) || (0 == _gConfiguration._internal_reward_mode);

        _gNumProcessors     = 1; //Runtime.getRuntime().availableProcessors();
        _gNumThreads        = _gConfiguration._no_of_threads;
        _gIsDedicatedMode   = _gConfiguration._is_dedicated_mode;
        _gTotalRandNum      = Math::Max(_gConfiguration._capacity, 4*1024*1024);
        _gThroughputTime    = _gConfiguration._throughput_time;

        _mon = new Hb(concurrent);
        //_mon = new LazyCounter(_gNumThreads, concurrent);

        //prepare the random numbers ...............................................
        System_err_println("");
        System_err_println("    START create random numbers.");
        PrepareActions();
        PrepareRandomNumbers(_gTotalRandNum);
        System_err_println("    END   creating random numbers.");
        System_err_println("");
        //System.gc();

        //run the benchmark ........................................................
        RunBenchmark();

        //print results ............................................................
        if(0 == _gConfiguration._is_dedicated_mode) {
                System_out_format(" *%4d*", (unsigned int)_gResult);
        } else {
                System_out_format(" *%4d* *%4d* *%4d*", (unsigned int)_gResultAdd, (unsigned int)_gResultRemove, (unsigned int)_gResultPeek);
        }
        Thread::sleep(1*1000);
        for (int iDb=0; iDb<_num_ds; ++iDb) {
                _gDS[iDb]->print_custom();
                printf(" **%d %d**", (_gEndTime - _gStartTime), _gResult * (_gEndTime - _gStartTime));
                delete _gDS[iDb]; 
                _gDS[iDb] = null;
        }
        //can't do this unless we delete all learning engines too else le tries to read from mon after deleted segfault
	//delete _mon;

}



////////////////////////////////////////////////////////////////////////////////
//HELPER FUNCTIONS
////////////////////////////////////////////////////////////////////////////////
void RunBenchmark() {
        //print test information ...................................................
        System_err_println("Benchmark Curr: ");
        System_err_println("--------------");
        System_err_println("    numOfThreads:           " + Integer::toString( _gConfiguration._no_of_threads));

        System_err_println("    Algorithm1 Name:        " + std::string(_gConfiguration._alg1_name));
        System_err_println("    Algorithm1 Num:         " + Integer::toString(_gConfiguration._alg1_num));
        System_err_println("    Algorithm2 Name:        " + std::string(_gConfiguration._alg2_name));
        System_err_println("    Algorithm2 Num:         " + Integer::toString(_gConfiguration._alg2_num));
        System_err_println("    Algorithm3 Name:        " + std::string(_gConfiguration._alg3_name));
        System_err_println("    Algorithm3 Num:         " + Integer::toString(_gConfiguration._alg3_num));
        System_err_println("    Algorithm4 Name:        " + std::string(_gConfiguration._alg4_name));
        System_err_println("    Algorithm4 Num:         " + Integer::toString(_gConfiguration._alg4_num));

        System_err_println("    NumProcessors:          " + Integer::toString(_gNumProcessors));
        System_err_println("    testNo:                 " + Integer::toString(_gConfiguration._test_no));

        System_err_println("    addOps:                 " + Integer::toString(_gConfiguration._add_ops));
        System_err_println("    removeOps:              " + Integer::toString(_gConfiguration._remove_ops));
        System_err_println("    throughput_time:        " + Integer::toString(_gConfiguration._throughput_time));
        System_err_println("    is_dedicated_mode:      " + Integer::toString(_gConfiguration._is_dedicated_mode));
        System_err_println("    tm_status:              " + Integer::toString(_gConfiguration._tm_status) + (std::string)("   (0=Norm; else View)"));
        System_err_println("    read_write_delay:       " + Integer::toString(_gConfiguration._read_write_delay));
        System_err_println("    fc_passes:              " + Integer::toString(_gConfiguration._fc_passes));
        System_err_println("    barrier_interval:       " + Integer::toString(_gConfiguration._barrier_interval));
        System_err_println("    scancount_tuning:       " + Integer::toString(_gConfiguration._scancount_tuning));
        System_err_println("    lock_scheduling:        " + Integer::toString(_gConfiguration._lock_scheduling));
        System_err_println("    dynamic_work_size:      " + Integer::toString(_gConfiguration._dynamic_work_size));
        System_err_println("    dynamic_work_intervals: " + Integer::toString(_gConfiguration._dynamic_work_intervals));
        System_err_println("    internal_reward_mode:   " + Integer::toString(_gConfiguration._internal_reward_mode));

        _is_view = (0 != _gConfiguration._tm_status);

        char _sprintf_str[1024];
        sprintf(_sprintf_str, "%f",  _gConfiguration._load_factor);
        System_err_println("    loadFactor:             " + (std::string)(_sprintf_str));
        sprintf(_sprintf_str, "%f",  _gConfiguration._rl_to_sleepidle_ratio);
        System_err_println("    rl_to_sleepidle_ratio:  " + (std::string)(_sprintf_str));

        System_err_println("    initialCapacity:        " + Integer::toString(_gConfiguration._capacity));
        System_err_println("");
        FCBase<FCIntPtr>::_num_post_read_write = _gConfiguration._read_write_delay;
        FCBase<FCIntPtr>::_num_passes = _gConfiguration._fc_passes;
        FCBase<FCIntPtr>::_sync_interval = _gConfiguration._barrier_interval;
        FCBase<FCIntPtr>::_dynamic_work_size = _gConfiguration._dynamic_work_size;
        FCBase<FCIntPtr>::_dynamic_work_intervals = _gConfiguration._dynamic_work_intervals;

        //create appropriate data-structure ........................................
	LearningEngine *learner;
        _num_ds=0;

        int mode = LearningEngine::disabled;
        int num_lock_sched = 0;
        int num_sc_tune = 0;
        if ( 0   != _gConfiguration._lock_scheduling )  { mode |= LearningEngine::lock_scheduling;  num_lock_sched = 1; }
        if ( 0   != _gConfiguration._scancount_tuning ) { mode |= LearningEngine::scancount_tuning; num_sc_tune    = 1; }
        if ( 1.0 != _gConfiguration._rl_to_sleepidle_ratio )   { mode |= LearningEngine::inject_delay; }


	if ( 0 == strncmp(_gConfiguration._alg1_name, "smart", 5) )
                learner = new LearningEngine(_gNumThreads, _mon, _gConfiguration._rl_to_sleepidle_ratio,
					     (LearningEngine::learning_mode_t) mode, 
                                             num_lock_sched*_gConfiguration._alg1_num, num_sc_tune*_gConfiguration._alg1_num );
	else
	        learner = null;

        for (int i=0; i<(_gConfiguration._alg1_num); ++i) {
	        FCBase<FCIntPtr>* tmp = CreateDataStructure(_gConfiguration._alg1_name, learner);
                if ( tmp == null ) {
                   System_err_println("Invalid data structure type requested: " + std::string(_gConfiguration._alg1_name));
                   exit(0);
                }
                _gDS[_num_ds++] = tmp;
        }


	if ( 0 == strncmp(_gConfiguration._alg2_name, "smart", 5) )
                learner = new LearningEngine(_gNumThreads, _mon, _gConfiguration._rl_to_sleepidle_ratio,
					     (LearningEngine::learning_mode_t) mode, 
                                             num_lock_sched*_gConfiguration._alg2_num, num_sc_tune*_gConfiguration._alg2_num );
	else
	        learner = null;


        for (int i=0; i<(_gConfiguration._alg2_num); ++i) {
	        FCBase<FCIntPtr>* tmp = CreateDataStructure(_gConfiguration._alg2_name, learner);
                if ( tmp == null ) {
                   System_err_println("Invalid data structure type requested: " + std::string(_gConfiguration._alg2_name));
                   exit(0);
                }
                _gDS[_num_ds++] = tmp;
        }


	if ( 0 == strncmp(_gConfiguration._alg3_name, "smart", 5) )
                learner = new LearningEngine(_gNumThreads, _mon, _gConfiguration._rl_to_sleepidle_ratio,
					     (LearningEngine::learning_mode_t) mode, 
                                             num_lock_sched*_gConfiguration._alg3_num, num_sc_tune*_gConfiguration._alg3_num );
	else
	        learner = null;

        for (int i=0; i<(_gConfiguration._alg3_num); ++i) {
	        FCBase<FCIntPtr>* tmp = CreateDataStructure(_gConfiguration._alg3_name, learner);
                if ( tmp == null ) {
                   System_err_println("Invalid data structure type requested: " + std::string(_gConfiguration._alg3_name));
                   exit(0);
                }
                _gDS[_num_ds++] = tmp;
        }


	if ( 0 == strncmp(_gConfiguration._alg4_name, "smart", 5) )
                learner = new LearningEngine(_gNumThreads, _mon, _gConfiguration._rl_to_sleepidle_ratio,
					     (LearningEngine::learning_mode_t) mode, 
                                             num_lock_sched*_gConfiguration._alg4_num, num_sc_tune*_gConfiguration._alg4_num );
	else
	        learner = null;

        for (int i=0; i<(_gConfiguration._alg4_num); ++i) {
	        FCBase<FCIntPtr>* tmp = CreateDataStructure(_gConfiguration._alg4_name, learner);
                if ( tmp == null ) {
                   System_err_println("Invalid data structure type requested: " + std::string(_gConfiguration._alg4_name));
                   exit(0);
                }
                _gDS[_num_ds++] = tmp;
        }


        //install work amount
        if ( 0 != FCBase<FCIntPtr>::_dynamic_work_size )
        {
                assert( _num_work_amts >= 1 );
                for(int i=0; i<_num_ds; ++i)
                       _gDS[i]->_num_post_read_write = _work_amts[0];
        }

        //calculate how much each thread should add the data-structure initially ...
        final int table_size  = (int)((_gConfiguration._capacity) * (_gConfiguration._load_factor));
        _g_thread_fill_table_size = table_size / _gNumThreads;

        //create benchmark threads .................................................
        System_err_println("    START creating threads.");
        _gThreads               =  new Thread*[_gNumThreads];
        _gThreadResultAry       =  new tick_t[_gNumThreads];
        memset((void*)_gThreadResultAry, 0, sizeof(int)*_gNumThreads);

        final int num_add_threads       = (int) Math::ceil(_gNumThreads * (_gConfiguration._add_ops)/100.0);
        final int num_remove_threads    = (int) Math::floor(_gNumThreads * (_gConfiguration._remove_ops)/100.0);
        final int num_peek_threads      = _gNumThreads - num_add_threads - num_remove_threads;

        if(0 == _gConfiguration._is_dedicated_mode) {
                for(int iThread = 0; iThread < _gNumThreads; ++iThread) {
                        _gThreads[iThread] =  new MixThread(iThread);
                }
        } else {
                System_err_println("    num_add_threads:    " + Integer::toString(num_add_threads));
                System_err_println("    num_remove_threads: " + Integer::toString(num_remove_threads));
                System_err_println("    num_peek_threads:   " + Integer::toString(num_peek_threads));

                int curr_thread=0;
                for(int iThread = 0; iThread < num_add_threads; ++iThread) {
                        _gThreads[curr_thread] =  new AddThread(curr_thread);
                        ++curr_thread;
                }
                for(int iThread = 0; iThread < num_remove_threads; ++iThread) {
                        _gThreads[curr_thread] =  new RemoveThread(curr_thread);
                        ++curr_thread;
                }
                for(int iThread = 0; iThread < num_peek_threads; ++iThread) {
                        _gThreads[curr_thread] =  new ContainThread(curr_thread);
                        ++curr_thread;
                }
        }
        System_err_println("    END   creating threads.");
        System_err_println("");
        Thread::yield();

        //start the benchmark threads ..............................................
        System_err_println("    START threads.");
        for(int iThread = 0; iThread < _gNumThreads; ++iThread) {
                _gThreads[iThread]->start();
        }
        System_err_println("    END START  threads.");
        System_err_println("");

        //wait the throughput time, and then signal the threads to terminate ...
        if ( 0 == FCBase<FCIntPtr>::_dynamic_work_size ) {
                Thread::yield();
                Thread::sleep(_gThroughputTime*1000);
        } else {
	        assert( FCBase<FCIntPtr>::_dynamic_work_intervals > 10 );
                Thread::yield();
                //Thread::sleep(0, _gThroughputTime * (1000000000 / FCBase<FCIntPtr>::_dynamic_work_intervals));
		Thread::delay(_gThroughputTime * (1000000000 / FCBase<FCIntPtr>::_dynamic_work_intervals));

                for(int i=1; i<FCBase<FCIntPtr>::_dynamic_work_intervals; ++i) {
                        for(int j=0; j<_num_ds; ++j)
                                _gDS[j]->_num_post_read_write = _work_amts[i%_num_work_amts];
                        //Thread::sleep(0, _gThroughputTime * (1000000000 / FCBase<FCIntPtr>::_dynamic_work_intervals));
 		        Thread::delay(_gThroughputTime * (1000000000 / FCBase<FCIntPtr>::_dynamic_work_intervals));
                }
        }
        _gIsStopThreads = 1;

        //join the threads .........................................................
        for(int iThread = 0; iThread < _gNumThreads; ++iThread) {
                Thread::yield();
                Thread::yield();
                _gThreads[iThread]->join();
        }
        System_err_println("    ALL threads terminated.");
        System_err_println("");

        //calculate threads results ................................................
        _gResult = 0;
        _gResultAdd = 0;
        _gResultRemove = 0;
        _gResultPeek = 0;
        if(0 == _gConfiguration._is_dedicated_mode) {
                for(int iThread = 0; iThread < _gNumThreads; ++iThread) {
                        _gResult += _gThreadResultAry[iThread];
                }
        } else {
                int curr_thread=0;
                for(int iThread = 0; iThread < num_add_threads; ++iThread) {
                        _gResultAdd += _gThreadResultAry[curr_thread];
                        ++curr_thread;
                }
                for(int iThread = 0; iThread < num_remove_threads; ++iThread) {
                        _gResultRemove += _gThreadResultAry[curr_thread];
                        ++curr_thread;
                }
                for(int iThread = 0; iThread < num_peek_threads; ++iThread) {
                        _gResultPeek += _gThreadResultAry[curr_thread];
                        ++curr_thread;
                }
        }

        //print benchmark results ..................................................
        for (int iDb=0; iDb<_num_ds; ++iDb) {
                System_err_println("    " + std::string(_gDS[iDb]->name()) + " Num elm: " + Integer::toString(_gDS[iDb]->size()));
        }

        //free resources ...........................................................
        delete [] _gRandNumAry;
        delete [] _gThreadResultAry;
        delete [] _gThreads;

        _gRandNumAry = null;
        _gThreadResultAry = null;
        _gThreads = null;

        //return benchmark results ................................................
        _gResult         /= (long)(_gEndTime - _gStartTime);
        _gResultAdd      /= (long)(_gEndTime - _gStartTime);
        _gResultRemove   /= (long)(_gEndTime - _gStartTime);
        _gResultPeek     /= (long)(_gEndTime - _gStartTime);
}

FCBase<FCIntPtr>* CreateDataStructure(char* final alg_name, LearningEngine* learner) {

        //queue ....................................................................
        if(0 == strcmp(alg_name, "fcqueue")) {
	        return (new SmartQueue<FCIntPtr,false,false>(null, null));
        }
        if(0 == strcmp(alg_name, "smartqueue")) {
	        if ( 0 != _gConfiguration._internal_reward_mode ) 
	                return (new SmartQueue<FCIntPtr,true,true>(_mon, learner));
		else
		        return (new SmartQueue<FCIntPtr,true,false>(_mon, learner));
        }
        if(0 == strcmp(alg_name, "msqueue")) {
                return (new MSQueue<FCIntPtr>());
        }
        if(0 == strcmp(alg_name, "basketsqueue")) {
                return (new BasketsQueue<FCIntPtr>());
        }
        //not working yet
        //if(0 == strcmp(alg_name, "ctqueue")) {
        //        return (new ComTreeQueue<FCIntPtr>());
        //}
        if(0 == strcmp(alg_name, "oyqueue")) {
                return (new OyamaQueue<FCIntPtr>());
        }
        //not working yet
        //if(0 == strcmp(alg_name, "oyqueuecom")) {
        //          return (new OyamaQueueCom<FCIntPtr>());
        //}

        //skiplist .................................................................
        if(0 == strcmp(alg_name, "fcskiplist")) {
	        return (new SmartSkipList<FCIntPtr,false,false>(null, null));
        }
        if(0 == strcmp(alg_name, "smartskiplist")) {
	        if ( 0 != _gConfiguration._internal_reward_mode ) 
	                return (new SmartSkipList<FCIntPtr,true,true>(_mon, learner));
		else
		        return (new SmartSkipList<FCIntPtr,true,false>(_mon, learner));
        }
        if(0 == strcmp(alg_name, "lfskiplist")) {
                return (new LFSkipList<FCIntPtr>());
        }
        if(0 == strcmp(alg_name, "lazyskiplist")) {
                return (new LazySkipList<FCIntPtr>());
        }

        //pairing heaps.............................................................
        if(0 == strcmp(alg_name, "fcpairheap")) {
	        return (new SmartPairHeap<FCIntPtr,false,false>(null, null));
        }
        if(0 == strcmp(alg_name, "smartpairheap")) {
	        if ( 0 != _gConfiguration._internal_reward_mode ) 
	                return (new SmartPairHeap<FCIntPtr,true,true>(_mon, learner));
		else
		        return (new SmartPairHeap<FCIntPtr,true,false>(_mon, learner));
        }

        //stack ....................................................................
        if(0 == strcmp(alg_name, "fcstack")) {
                return (new FCStack<FCIntPtr>());
        }
        //if(0 == strcmp(alg_name, "smartstack")) {
        //      return (new SmartStack<FCIntPtr>(_mon));
        //}
        if(0 == strcmp(alg_name, "lfstack")) {
                return (new LFStack<FCIntPtr>());
        }
        if(0 == strcmp(alg_name, "elstack")) {
                return (new EliminationStack<FCIntPtr>());
        }

        if(0 == strcmp(alg_name, "heartbeat")) {
	        return (new Hb());
        }
        if(0 == strcmp(alg_name, "lazycounter")) {
	        return (new LazyCounter(_gConfiguration._no_of_threads));
        }

        return null;
}


void PrepareActions() {
        final int add_limit = _gConfiguration._add_ops;
        final int remove_limit = add_limit + _gConfiguration._remove_ops;

        for(int iAction=0; iAction < _ACTIONS_ARY_SIZE; ++iAction) {
                final int rand_num = _gRand.nextInt(1024*1024)%100;

                if(rand_num < 0 || rand_num >= 100) {
                        System_err_println("PrepareActions: Error random number" + Integer::toString(rand_num));
                        System::exit(1);
                }

                if(rand_num < add_limit)
                        _gActionAry[iAction] = 1;
                else if(rand_num < remove_limit)
                        _gActionAry[iAction] = 2;
                else
                        _gActionAry[iAction] = 3;
        }
}
void PrepareRandomNumbers(final int size) {
        _gRandNumAry = new int[size];
        for (int iRandNum = 0; iRandNum < size; ++iRandNum) {
                if(0 == _gConfiguration._capacity) {
                        _gRandNumAry[iRandNum] = iRandNum+2;
                } else { 
                        _gRandNumAry[iRandNum] = _gRand.nextInt(_gConfiguration._capacity) + 2;
                        if(_gRandNumAry[iRandNum] <= 0 || _gRandNumAry[iRandNum] >= (_gConfiguration._capacity+2)) {
                                System_err_println("PrepareRandomNumbers: Error random number" + Integer::toString(_gRandNumAry[iRandNum]));
                                System::exit(1);
                        }
                }
        }
}

int NearestPowerOfTwo(final int x) {
        int mask = 1;
        while(mask < x) {
                mask <<= 1;
        }
        return mask;
}
