
#ifndef ECMP_H
#define ECMP_H

#include <assert.h>
#include "typedef.h"


const dist_t INF = (10000000);

const unsigned INITIAL_RATE = 1;

// initialize d[i][j] with w[k] or INF
void initialize(dist_t* d, flag_t* p, dist_t* w);
void ecmpWarshallfloyd(dist_t* d, flag_t* p);
void buildNexthopList(node_t i, node_t j, flag_t* p, flag_t* visited);

void allocFlow(node_t h, node_t d, flow_t flow, unsigned rate, flag_t* p, flow_t* f);


#ifdef TIME_MEASUREMENT
Timer init, ecmpwf, nextlist, trace, calcl;

void showTimeStatisticsEcmp();
#endif


void calculateFlowFromWeight(double& L, dist_t w[E]);








#endif


