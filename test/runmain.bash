#!/bin/bash

#////////////////////////////////////////////////////////////////////////////// 
#// File    : runmain.bash                                          
#// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com  
#//           Ms.Moran Tzafrir  email: morantza@gmail.com  
#// Written : 16 February 2011, 13 October 2009 
#//                                                         
#// Bash script for running the multi-platform benchamrk main.cpp
#//                                                       
#// Copyright (C) 2011 Jonathan Eastep, 2009 Moran Tzafrir. 
#//                                                                             
#// This program is free software; you can redistribute it and/or modify 
#// it under the terms of the GNU General Public License as published by 
#// the Free Software Foundation; either version 2 of the License, or 
#// (at your option) any later version.  
#//                                    
#// This program is distributed in the hope that it will be useful, but 
#// WITHOUT ANY WARRANTY; without even the implied warranty of       
#// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
#// General Public License for more details.                    
#//                                                        
#// You should have received a copy of the GNU General Public License   
#// along with this program; if not, write to the Free Software Foundation 
#// Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA       
#//////////////////////////////////////////////////////////////////////////////


#sets the ns between data structure ops (aka load or arrival rate) in the benchmark
delays="6400 3200 1600 800 400 200 100"
#delays="400 200 100"
#delays="800"
#use for heartbeats benchmarks
#delays="7000 3500 2333 1750 1400 1167 1000"
#use for lazycounter benchmarks
#delays="6400 3200 1600 800 400 200 100"

#sets the level of concurrency 
#threads="2 4 6 8 10 12 14"
threads="14"
#use for monitor benchmarks
#threads="15"

#which algorithms to benchmark
algorithms="fcqueue fcskiplist fcpairheap smartqueue smartskiplist smartpairheap msqueue basketsqueue basketsqueue oyqueue oyqueuecom lfskiplist lazyskiplist"
#algorithms="smartqueue smartskiplist smartpairheap"
#algorithms="smartskiplist"
#algorithms="smartqueue smartskiplist smartpairheap"
#use for monitor benchmarks
#algorithms="heartbeat"
#algorithms="heartbeat lazycounter"

#the setting to use for the flat combining scancount parameter
#use 0 if you want it to be nthreads
#scancounts="1 3 5 7 9 11 13 15 17 19 21 23 25"
#scancounts="1 5 9 13 17 21 25 29 33 37 41 45 49"
scancounts="0"

#currently deprecated
#syncintervals="1000000 100000 10000 1000 100 10 1"
syncintervals="0"

#use these for 3 producers, 11 consumers
#dedicated="1"
#addops="21"
#removeops="79"
#use these for 2 producers, 12 consumers
#dedicated="1"
#addops="14"
#removeops="86"
#use these for 1 producer, (threads-1) consumers
dedicated="1"
addops="1"
removeops="99"
#use these for monitor benchmarks
#dedicated="1"
#addops="0"
#removeops="99"
#use these for random insert / delete from all threads
#dedicated="0"
#addops="50"
#removeops="50"

#configure this to have a changing amount of postcomputation work (in ns)
#dynamicworkamt="1"
#dynamicworkintervals="1 10 100 1000 10000"
#dynamicworkintervals="1000 10000 100000 1000000 10000000"
#dynamicworkintervals="1000000"
#use this for none
dynamicworkamt="0"
dynamicworkintervals="0"

#machine learning enablement
#scancount tuning off
#scancounttuning="0"
#lockscheduling="0"
#scancount tuning on
scancounttuning="1"
lockscheduling="0"

#configure to simulate slowdown of rl thread 
#rltime / totaltime
#rltosleepidleratios="1.0 .25 .0625 .015625 .00390625 .000976562 .000244141"
#rltosleepidleratios=".03125 .015625"
rltosleepidleratios="1.0"

#how many data structures to instantiate
#normal
numds="1"
#more than 1
#numds="1 2 4 8 16"
#numds="2 4 8 16"

#configure whether to use internal or external reward
#internalreward="1"
#use this for app-specific reward
internalreward="1"

#how many trials of each experiment to perform
reps="0 1 2 3 4 5 6 7 8 9"

#misc parameters
test="TIME"
capacities="20000"
count=0

rm -f $test

for algorithm in $algorithms; do

for num in $numds; do

for syncinterval in $syncintervals; do

for interval in $dynamicworkintervals; do

for delay in $delays; do

for capacity in $capacities; do

for thread in $threads; do

for scancount in $scancounts; do

for rltosleepidleratio in $rltosleepidleratios; do

	count=$(($count + 1))

        line=""
	if [ "$scancount" = "0" ]; then
                 line="$algorithm $num non 0 non 0 non 0 $count $thread $addops $removeops 0.0 $capacity 10 $dedicated 0 $delay $thread $syncinterval $scancounttuning $lockscheduling $dynamicworkamt $interval $rltosleepidleratio $internalreward"
	else
                 line="$algorithm $num non 0 non 0 non 0 $count $thread $addops $removeops 0.0 $capacity 10 $dedicated 0 $delay $scancount $syncinterval $scancounttuning $lockscheduling $dynamicworkamt $interval $rltosleepidleratio $internalreward"
	fi

        for rep in $reps; do
                 echo "$line" 1>&2
		 echo -n "$line" >> $test
		 ./main_intel64 $line >> $test
		 echo "" >> $test
        done;

done; done; done; done; done; done; done; done; done;