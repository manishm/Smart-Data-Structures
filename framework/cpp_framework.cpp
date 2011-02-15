#include "cpp_framework.h"

//------------------------------------------------------------------------------
// START
// File    : cpp_framework.cpp
// Authors : Jonathan Eastep   email: jonathan.eastep@gmail.com
//           Ms.Moran Tzafrir  email: morantza@gmail.com
// Written : 16 February 2011, 13 April 2009
// 
// Copyright (C) 2009 Moran Tzafrir.
// You can use this file only by explicit written approval from Jonathan Eastep
// per Moran's original license
//------------------------------------------------------------------------------

const int CCP::Integer::MIN_VALUE = INT_MIN;
const int CCP::Integer::MAX_VALUE = INT_MAX;
const int CCP::Integer::SIZE = 32;

__thread__	CCP::Thread* CCP::_g_tls_current_thread = null;
