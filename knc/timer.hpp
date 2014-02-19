/*
 *  timer.h
 *  
 *  Common timer functionality
 *
 *  Created by John Linford on 4/8/08.
 *  Copyright 2008 Transatlantic Giraffe. All rights reserved.
 *
 */

#ifndef __TIMER_H__
#define __TIMER_H__

#include <iostream>
#include <string>
#include <map>
#include <sys/time.h>

namespace fixedgrid {

class Timer
{
public:

	Timer() : t0(0), elapsed(0) { }

	void start() {
		t0 = usec();
	}

	void stop() {
		elapsed += usec() - t0;
	}

	double value() const {
		return elapsed;
	}

private:

	static double usec();

    double t0;
    double elapsed;
};


class Metrics
{
public:

	typedef std::map<std::string, Timer> timer_map_t;

	Metrics(std::string _name) : name(_name)
	{ }

	Timer & operator[](std::string const & name) {
		return timers[name];
	}

	Timer const & operator[](std::string const & name) const {
		return timers.at(name);
	}

	friend std::ostream & operator<<(std::ostream & os, Metrics const & m);

private:

	std::string name;
	timer_map_t timers;
};

} // namespace fixedgrid

#endif
