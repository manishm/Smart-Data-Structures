//FIXME: This file is woefully poor at regression testing. 

//////////////////////////////////////////////////////////////////////////////// 
// File    : regress.cpp                                                        
// Author  : Jonathan Eastep   email: jonathan.eastep@gmail.com                  
// Written : 16 February 2011
//                                        
// Copyright (C) 2011 Jonathan Eastep
//
// Multi-Platform regression test for Smart Data Structures
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

#include <iostream>
#include <string>
#include <iomanip>
#include <pthread.h>
#include <stdlib.h>
#include <sys/time.h>
#include <errno.h>

#include "Heartbeat.h"
#include "FCQueue.h"
#include "SmartQueue.h"
#include "MSQueue.h"
#include "BasketsQueue.h"
#include "ComTreeQueue.h"
#include "OyamaQueue.h"
#include "OyamaQueueCom.h"
#include "FCSkipList.h"
#include "SmartSkipList.h"
#include "LFSkipList.h"
#include "LazySkipList.h"
#include "FCPairingHeap.h"
#include "SmartPairingHeap.h"
#include "FCStack.h"
#include "LFStack.h"
#include "EliminationStack.h"
#include "cpp_framework.h"

using namespace CCP;
using namespace std;

#define NUMTRIALS (20)

FCBase<FCIntPtr>*  queue;
FCBase<FCIntPtr>*  queue2;

typedef FCIntPtr      lli;

//#define SMART
#define QUEUE2

volatile _u64  barrier  = 0;
volatile _u64  barrier2 = 0;

int        _gNumThreads       = 14;
int        _gIsDedicatedMode  = 0;
const int  MAX_THREADS        = 200;
const int  els                = (_gNumThreads-1) * 100;
const int  MAX_ELS            = MAX_THREADS * 100;

char            pad1[CACHE_LINE_SIZE];
volatile char   results1[CACHE_LINE_SIZE*MAX_ELS] = {false};
volatile char   results2[CACHE_LINE_SIZE*MAX_ELS] = {false};
char            pad2[CACHE_LINE_SIZE];


void Sleep(unsigned int nsec)
{
        struct timespec timeout0;
	struct timespec timeout1;
	struct timespec* tmp;
	struct timespec* t0 = &timeout0;
	struct timespec* t1 = &timeout1;

	t0->tv_sec = nsec / 1000000000;
	t0->tv_nsec = (nsec % 1000000000);

	while ((nanosleep(t0, t1) == (-1)) && (errno == EINTR))
	{
                tmp = t0;
		t0 = t1;
		t1 = tmp;
	}
}

void SetupRandSleep()
{
        srand( time(NULL) );
}

void RandSleep()
{
        unsigned int range = 20000;
        unsigned int thetime = rand() % range + 1;
        Sleep(thetime);
}

void * func(void* args)
{
        ptr_t tid = (ptr_t) args;
	//cerr << "tid = " << tid << endl;

	FAADD(&barrier, 1);
	while( barrier != _gNumThreads );

	int total = (els/(_gNumThreads-1));
	//cerr << "total els = " << total << endl;

	for(int i = 0; i < total; i++)
	{
	        RandSleep();
                //volatile so the compiler doesn't optimize it out
                FCIntPtrNode* volatile res;  
		do {
		      res = (FCIntPtrNode*) queue->remove(tid, NULL);
		      //cerr << dec << "outValue1 = " << (res->getvalue()) << endl;	  
		} while (NULL == res);
		//cerr << dec << "outValue1 = " << (res->getvalue()) << endl;
		results1[CACHE_LINE_SIZE*(res->getvalue()-1)] = true;
                //cerr << "freeing " << res << endl;
		delete res;

#ifdef QUEUE2
		RandSleep();
		do {
		        res = (FCIntPtrNode*) queue2->remove(tid, NULL);
			//cerr << dec << "outValue2 = " << (res->getvalue()) << endl;
		} while (NULL == res);
		//cerr << dec << "outValue2 = " << (res->getvalue()) << endl;
		results2[CACHE_LINE_SIZE*(res->getvalue()-1)] = true;
                //cerr << "freeing " << res << endl;
		delete res;
#endif
	}

	//cerr << "func finished" << endl;
	//cerr.flush();

	FAADD(&barrier2, 1);
	while( barrier2 != _gNumThreads );   
}


bool parallel_test(FCBase<FCIntPtr>* queue, FCBase<FCIntPtr>* queue2)
{
        if ( queue == null ) {
	        cerr << "Failed parallel test of null data structure" << endl;
		return false;
	}

#ifdef QUEUE2
	if ( queue2 == null ) {
	        cerr << "Failed parallel test of null data structure" << endl;
		return false;
	}
#endif

	bool rv = true;

	try {

	        SetupRandSleep();

		pthread_attr_t managerthreadattr;
		static pthread_t managerthread[MAX_THREADS];

		pthread_attr_init(&managerthreadattr);
		pthread_attr_setdetachstate(&managerthreadattr, PTHREAD_CREATE_JOINABLE);

		for(int i = 1; i < _gNumThreads; i++)
		        pthread_create(&managerthread[i], &managerthreadattr, func, (void*) i);

		FAADD(&barrier, 1);
		while( barrier != _gNumThreads );


		for(int i = (0+1); i < (els+1); i++)
		{
		        RandSleep();
			queue->add(0, new FCIntPtrNode(i));
#ifdef QUEUE2
			queue2->add(0, new FCIntPtrNode(i));
#endif
		}

		FAADD(&barrier2, 1);
		while( barrier2 != _gNumThreads ); 

		for(int i = 1; i < _gNumThreads; i++)
		        pthread_join(managerthread[i], NULL);

		// test correctness
		for(int i = 0; i < els; i++)
		{  
		        if ( results1[i*CACHE_LINE_SIZE] == false )
			        rv = false;
			if ( results2[i*CACHE_LINE_SIZE] == false )
			        rv = false;
		}

		// reset
		barrier = 0;
		barrier2 = 0;
		for(int i = 0; i < els; i++)
		{  
		        results1[i*CACHE_LINE_SIZE] = false;
			results2[i*CACHE_LINE_SIZE] = false;
		}

	}

	catch (...) {
	        rv = false;
	}

	if ( rv )
	        cerr << "Passed parallel test of " << queue->name() << endl;
	else
	        cerr << "Failed parallel test of " << queue->name() << endl;

	return rv;
}


enum TESTTYPE {
        FIFO,
	PRI,
	LIFO
};

bool serial_test(FCBase<FCIntPtr>* queue, TESTTYPE type)
{
        if ( queue == null ) {
	        cerr << "Failed serial test of null data structure" << endl;
		return false;
	}

	const int  NUMVALS             = 8;
	FCIntPtr   vals[NUMVALS]       = {1000, (FCIntPtr) &vals[1], 665, 1001, 22, 110, (FCIntPtr) vals, 2000};
	FCIntPtr   sortedvals[NUMVALS] = {22, 110, 665, 1000, 1001, 2000, (FCIntPtr) vals, (FCIntPtr) &vals[1]};

	for(int i = 0; i < NUMVALS; i++)
	        queue->add(0, new FCIntPtrNode( vals[i] ));    

	int errs = 0;
	FCIntPtrNode* res; 

	try {
	        switch( type ) {
		case PRI:
		        for(int i = 0; i < NUMVALS; i++) {
			        res = (FCIntPtrNode*) queue->remove(0, NULL);
			        errs += (res->getvalue() == sortedvals[i]) ? 0 : 1;
				delete res;
		        }
			errs += (queue->remove(0, NULL) == FCBase<FCIntPtr>::_NULL_VALUE) ? 0 : 1;    
			break;
		case FIFO:
		        for(int i = 0; i < NUMVALS; i++) {
			        res = (FCIntPtrNode*) queue->remove(0, NULL);
			        errs += (res->getvalue() == vals[i]) ? 0 : 1;
				delete res;
			}
			errs += (queue->remove(0, NULL) == FCBase<FCIntPtr>::_NULL_VALUE) ? 0 : 1;        
			break;
		case LIFO:
		        for(int i = 0; i < NUMVALS; i++) {
			        res = (FCIntPtrNode*) queue->remove(0, NULL);			  
			        errs += (res->getvalue() == vals[NUMVALS-1-i]) ? 0 : 1;
				delete res;
			}
			errs += (queue->remove(0, NULL) == FCBase<FCIntPtr>::_NULL_VALUE) ? 0 : 1;  
			break;
		default:
		        break;
		}
	}
	catch (...) {
	        errs++;
	}

	bool rv = (errs == 0);
	if ( rv )
	        cerr << "Passed serial test of " << queue->name() << endl;
	else
	        cerr << "Failed serial test of " << queue->name() << endl;

	return rv;
}

bool destruct_test(FCBase<lli>** ds1, FCBase<lli>** ds2, int num_ds, Hb* hbmon)
{
        bool rv = true;
        try {

                for(int i = 0; i < num_ds; i++) {
                        delete ds1[i];
                        delete ds2[i];
                }
		delete hbmon;
        }
        catch(...) {
	        cerr << "Failed destruct test" << endl;
	        rv = false;
        }
	return rv;
}


int main(int argc, char* argv[])
{

        Hb*                   hbmon = new Hb();

        const int             NUMDS = 16;

        FCBase<lli>* ds1[NUMDS] = { new FCQueue<lli>(), new SmartQueue<lli>(hbmon), new MSQueue<lli>(), new BasketsQueue<lli>(),
                                    null /*new ComTreeQueue<lli>()*/, new OyamaQueue<lli>(), null /*new OyamaQueueCom<lli>()*/,
                                    new FCSkipList<lli>(), new SmartSkipList<lli>(hbmon), new LFSkipList<lli>(),
                                    new LazySkipList<lli>(), new FCPairHeap<lli>(), new SmartPairHeap<lli>(hbmon),
                                    new FCStack<lli>(), new LFStack<lli>(), new EliminationStack<lli>() };

#ifdef QUEUE2
        FCBase<lli>* ds2[NUMDS] = { new FCQueue<lli>(), new SmartQueue<lli>(hbmon), new MSQueue<lli>(), new BasketsQueue<lli>(),
                                    null /*new ComTreeQueue<lli>()*/, new OyamaQueue<lli>(), null /*new OyamaQueueCom<lli>()*/,
                                    new FCSkipList<lli>(), new SmartSkipList<lli>(hbmon), new LFSkipList<lli>(),
                                    new LazySkipList<lli>(), new FCPairHeap<lli>(), new SmartPairHeap<lli>(hbmon),
                                    new FCStack<lli>(), new LFStack<lli>(), new EliminationStack<lli>() };
#else
        FCBase<lli>* ds2[NUMDS] = {null};
#endif

        int PRISTART = 7;
        int STACKSTART = 13;


        Memory::read_write_barrier();
        
        bool megatotal = true;
        int megafails = 0;
	int megaskips = 0;

        for(int j = 0; j < NUMDS; j++) {

                bool rv, stotal = true, ptotal = true;

                queue = ds1[j];
                queue2 = ds2[j];

                for(int i = 0; i < NUMTRIALS; i++) {
                        TESTTYPE type = (j < PRISTART) ? FIFO : ((j < STACKSTART) ? PRI : LIFO);
                        rv = serial_test(queue, type);
                        megafails += rv ? 0 : 1;
                        megaskips += (null == queue) ? 1 : 0;
                        stotal &= rv;
                }
                if ( stotal )
                        cerr << "Passed all serial tests of " << queue->name() << endl;

                for(int i = 0; i < NUMTRIALS; i++) {
                        rv = parallel_test(queue, queue2);
                        megafails += rv ? 0 : 1;
#ifdef QUEUE2
                        megaskips += ((null == queue) || (null == queue2)) ? 1 : 0;
#else
                        megaskips += (null == queue) ? 1 : 0;
#endif
                        ptotal &= rv;
                }
                if ( ptotal )
                        cerr << "Passed all parallel tests of " << queue->name() << endl;

                megatotal &= stotal;
                megatotal &= ptotal;
        }

        bool rv = destruct_test(ds1, ds2, NUMDS, hbmon);
        megafails += rv ? 0 : 1;
        megatotal &= rv;

        cerr << "Passed object destruction test" << endl;

        if ( megatotal )
                cerr << "Passed all tests" << endl;
        else
                cerr << "Failed " << megafails << " out of " << (NUMDS*NUMTRIALS*2) 
                     << " tests. " << megaskips << " were due to skips. Check output" << endl;

}
