#ifndef __CONFIGURATION_HPP__
#define __CONFIGURATION_HPP__

//------------------------------------------------------------------------------
// START
// File    : Configuration.h
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com
//           Ms.Moran Tzafrir  email: morantza@gmail.com
// Written : 13 April 2009
// 
// Configuration 
//
// Copyright (C) 2011 Jonathan Eastep, 2009 Moran Tzafrir.
//
// You can use this file only by explicit written approval from Jonathan Eastep
// per Moran's original license.
//------------------------------------------------------------------------------

#include <string>
#include "cpp_framework.h"

class Configuration {
public:
        //..................................................
        char    _alg1_name[1024];
        int     _alg1_num;

        char    _alg2_name[1024];
        int     _alg2_num;

        char    _alg3_name[1024];
        int     _alg3_num;

        char    _alg4_name[1024];
        int     _alg4_num;

        //..................................................
        int     _test_no;
        int     _no_of_threads;

        int     _add_ops;
        int     _remove_ops;

        float   _load_factor;
        int     _capacity;

        int     _throughput_time;
        int     _is_dedicated_mode;

        int     _tm_status;
        int     _read_write_delay;
        int     _fc_passes;
        int     _barrier_interval;
        int     _scancount_tuning;
        int     _lock_scheduling;
        int     _dynamic_work_size;
        int     _dynamic_work_intervals;

        float   _rl_to_sleepidle_ratio;
        int     _internal_reward_mode;

        //..................................................
        bool read() {
                try {
                        //read configuration from input stream
                        int num_read = fscanf(stdin, "%s %d %s %d %s %d %s %d %d %d %d %d %f %d %d %d %d %d %d %d %d %d %d %f", 
                                              _alg1_name, &_alg1_num, _alg2_name, &_alg2_num, _alg3_name, &_alg3_num, _alg4_name, &_alg4_num,
                                              &_test_no, &_no_of_threads, &_add_ops, &_remove_ops, &_load_factor, &_capacity, &_throughput_time, 
                                              &_is_dedicated_mode, &_tm_status, &_read_write_delay, &_fc_passes, &_barrier_interval, 
                                              &_scancount_tuning, &_lock_scheduling, &_dynamic_work_size, &_dynamic_work_intervals,
                                              &_rl_to_sleepidle_ratio, &_internal_reward_mode );

                        return (26 == num_read);
                } catch (...) {
                        return false;
                }
        }

        //..................................................
        bool read(int argc, char **argv) {
                try {
                        int curr_arg=1;

                        strcpy(_alg1_name, argv[curr_arg++]);
                        _alg1_num               = CCP::Integer::parseInt(argv[curr_arg++]);
                        strcpy(_alg2_name, argv[curr_arg++]);
                        _alg2_num               = CCP::Integer::parseInt(argv[curr_arg++]);
                        strcpy(_alg3_name, argv[curr_arg++]);
                        _alg3_num               = CCP::Integer::parseInt(argv[curr_arg++]);
                        strcpy(_alg4_name, argv[curr_arg++]);
                        _alg4_num               = CCP::Integer::parseInt(argv[curr_arg++]);

                        _test_no                = CCP::Integer::parseInt(argv[curr_arg++]);
                        _no_of_threads          = CCP::Integer::parseInt(argv[curr_arg++]);
                        _add_ops                = CCP::Integer::parseInt(argv[curr_arg++]);
                        _remove_ops             = CCP::Integer::parseInt(argv[curr_arg++]);
                        _load_factor            = (float)atof(argv[curr_arg++]);
                        _capacity               = CCP::Integer::parseInt(argv[curr_arg++]);
                        _throughput_time        = CCP::Integer::parseInt(argv[curr_arg++]);
                        _is_dedicated_mode      = CCP::Integer::parseInt(argv[curr_arg++]);
                        _tm_status              = CCP::Integer::parseInt(argv[curr_arg++]);
                        _read_write_delay       = CCP::Integer::parseInt(argv[curr_arg++]);
                        _fc_passes              = CCP::Integer::parseInt(argv[curr_arg++]);
                        _barrier_interval       = CCP::Integer::parseInt(argv[curr_arg++]);
                        _scancount_tuning       = CCP::Integer::parseInt(argv[curr_arg++]);
                        _lock_scheduling        = CCP::Integer::parseInt(argv[curr_arg++]);
                        _dynamic_work_size      = CCP::Integer::parseInt(argv[curr_arg++]);
                        _dynamic_work_intervals = CCP::Integer::parseInt(argv[curr_arg++]);
                        _rl_to_sleepidle_ratio  = (float)atof(argv[curr_arg++]);
                        _internal_reward_mode   = CCP::Integer::parseInt(argv[curr_arg++]);

                        return true;
                } catch (...) {
                        return false;
                }
        }
};

#endif
