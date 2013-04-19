
#include "../include/timer.h"

struct timeval tv_subtract(struct timeval a, struct timeval& b){
	a.tv_usec -= b.tv_usec; 	a.tv_sec -= b.tv_sec;	return a;
}

double tv_to_usec(struct timeval a){
	return a.tv_usec + 1e+6 * a.tv_sec;
}

double tv_to_sec(struct timeval a){
	return a.tv_sec + 1e-6 * a.tv_usec;
}

double Timer::sec(){
	return (double)1e-6 * total_us; }

double Timer::stop(){
	if(started){
		static struct timeval stop_time;
		gettimeofday(&stop_time, NULL);
		total_us += tv_to_usec(tv_subtract(stop_time, start_time));
		started = false;
	}
	return sec();
}

double Timer::now(){
	if(started){
		static struct timeval stop_time;
		gettimeofday(&stop_time, NULL);
		return tv_to_sec(tv_subtract(stop_time, start_time)) + sec();
	}else{
		return sec();
	}
	
}

void Timer::start(){
	started = true;
	gettimeofday(&start_time, NULL);
}

Timer::timer(){ gettimeofday(&start_time, NULL);	total_us = 0;	started = false;}


