
#ifndef TIMER_H
#define TIMER_H
#include <sys/time.h>

struct timeval tv_subtract(struct timeval a, struct timeval& b){
	a.tv_usec -= b.tv_usec; 	a.tv_sec -= b.tv_sec;	return a;
}

double tv_to_usec(struct timeval a){
	return a.tv_usec + 1e+6 * a.tv_sec;
}

double tv_to_sec(struct timeval a){
	return a.tv_sec + 1e-6 * a.tv_usec;
}

typedef struct timer{
	
	//varuables
	struct timeval start_time;
	unsigned long total_us;
	int started;

	//functions
	double sec(){
		return (double)1e-6 * total_us; }
		
	double stop(){
		if(started){
			static struct timeval stop_time;
			gettimeofday(&stop_time, NULL);
			total_us += tv_to_usec(tv_subtract(stop_time, start_time));
			started = false;
		}
		return sec();
	}
	
	double now(){
		if(started){
			static struct timeval stop_time;
			gettimeofday(&stop_time, NULL);
			return tv_to_sec(tv_subtract(stop_time, start_time)) + sec();
		}else{
			return sec();			
		}

	}

	void start(){
		started = true;
		gettimeofday(&start_time, NULL);
	}
	
	timer(){ gettimeofday(&start_time, NULL);	total_us = 0;	started = false;}
	
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


