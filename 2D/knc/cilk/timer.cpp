/*
 *  timer.cpp
 *  
 *  Common timer functionality
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#include <ctime>
#include "timer.hpp"

using namespace std;

#ifdef __MIC__
Metrics __global_metrics("Device");
#else
Metrics __global_metrics("Host");
#endif

double Timer::usec()
{
  static int sec = -1;
  struct timeval tv;
  gettimeofday(&tv, 0);
  if(sec < 0) sec = tv.tv_sec;
  return (tv.tv_sec - sec)*1e6 + tv.tv_usec;
}

ostream & operator<<(ostream & os, Metrics const & m)
{
  os << endl << "======== " << m.name << " ========" << endl;

  typedef Metrics::timer_map_t::const_iterator iterator_t;
  for (iterator_t it = m.timers.begin(); it != m.timers.end(); it++) {
    string const & name = it->first;
    Timer const & t = it->second;
    os << name << ": " << (1e-6*t.value()) << " seconds" << endl;
  }
  return os;
}

