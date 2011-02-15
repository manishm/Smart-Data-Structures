#!/bin/bash

#////////////////////////////////////////////////////////////////////////////////             
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
#////////////////////////////////////////////////////////////////////////////////


#sets the ns between data structure ops (aka load or arrival rate) in the benchmark
#delays="6400 3200 1600 800 400 200 100"
delays="800"
#delays="12800 6400 3200 1600 800 400 200 100"

#sets the level of concurrency 
#threads="2 4 6 8 10 12 14"
threads="14"

#which algorithms to benchmark
algorithms="fcqueue fcskiplist fcpairheap smartqueue smartskiplist smartpairheap msqueue basketsqueue ctqueue lfskiplist lazyskiplist lazyskiplist fcstack lfstack elstack"
#algorithms="smartqueue smartskiplist smartpairheap"
#algorithms="smartpairheap"
#algorithms="smartqueue"

#the setting to use for the flat combining scancount parameter
#use 0 if you want it to be nthreads
#scancounts="1 3 5 7 9 11 13 15 17 19 21 23 25"
scancounts="0"

#currently deprecated
#syncintervals="1000000 100000 10000 1000 100 10 1"
syncintervals="0"

test="TIME"
capacities="20000"

#use these for 1 producer, (threads-1) consumers
dedicated="1"
addops="1"
removeops="99"
#use these for random insert / delete from all threads
#dedicated="0"
#addops="50"
#removeops="50"

#configure this to have a changing amount of postcomputation work each second
#dynamicworkamt="1"
#dynamicworkintervals="1 10 100 1000 10000"
#use this for none
dynamicworkamt="0"
dynamicworkintervals="0"

#machine learning enablement
scancounttuning="1"
lockscheduling="0"

#how many trials of each experiment to perform
reps="0 1 2 3 4 5 6 7 8 9"


count=0

rm -f $test

for algorithm in $algorithms; do

for syncinterval in $syncintervals; do

for interval in $dynamicworkintervals; do

for delay in $delays; do

for capacity in $capacities; do

for thread in $threads; do

for scancount in $scancounts; do

	count=$(($count + 1))

        line=""
	if [ "$scancount" = "0" ]; then
                 line="$algorithm 1 non 0 non 0 non 0 $count $thread $addops $removeops 0.0 $capacity 10 $dedicated 0 $delay $thread $syncinterval $scancounttuning $lockscheduling $dynamicworkamt $interval"
	else
                 line="$algorithm 1 non 0 non 0 non 0 $count $thread $addops $removeops 0.0 $capacity 10 $dedicated 0 $delay $scancount $syncinterval $scancounttuning $lockscheduling $dynamicworkamt $interval"
	fi

        for rep in $reps; do
                 echo "$line" 1>&2
		 echo -n "$line" >> $test
		 ./main_intel64 $line >> $test
		 echo "" >> $test
        done;

done; done; done; done; done; done; done;