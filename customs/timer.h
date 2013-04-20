
#ifndef TIMER_H
#define TIMER_H

#include <sys/time.h>

struct timeval tv_subtract(struct timeval a, struct timeval& b);

double tv_to_usec(struct timeval a);

double tv_to_sec(struct timeval a);

typedef struct timer{
	
	//varuables
	struct timeval start_time;
	unsigned long total_us;
	int started;

	//functions
	double sec();
	double stop();
	
	double now();
	void start();
	
	timer();
	
} Timer;

#ifdef TIME_MEASUREMENT
#	define TIMER_START(timer_name) timer_name.start()
#	define TIMER_STOP(timer_name) timer_name.stop()
#	define TIMER_SEC(timer_name) timer_name.sec()
#else
#	define TIMER_START(timer_name)
#	define TIMER_STOP(timer_name) (0)
#	define TIMER_SEC(timer_name) (0.0)
#endif

#endif


