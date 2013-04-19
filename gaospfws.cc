


#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#ifndef N
#   include "sample/w8_10.h"
#endif
#include <math.h>
#include <time.h>
#include <omp.h>




#ifndef A
#	define A 25
#endif
#ifndef B
#	define B 300
#endif
#ifndef C
#	define C 25
#endif

#   define POP (A+B+C)        // = B + C
#	define P (B+C)
#	define ALPHA ((double)A/POP)
#	define BETA ((double)A/POP)

#ifndef M
#   define M 0.001
#endif
#ifndef K
#   define K 0.8
#endif

#ifndef WMAX
#   define WMAX 20
#endif


#ifndef GENMAX
#   define GENMAX 1000000
//#   define GENMAX 1000
#endif
#ifdef DEBUG
#   undef GENMAX
#   define GENMAX 1000
#endif


#ifndef TMLIM
#	define TMLIM 3600
#endif
//#define GOAL 0.04


//-----------other macro

#define PRINT
//#define TIME_MEASUREMENT


//#define DEBUG 1

//-----------list datatype

typedef unsigned long flag_t;
typedef unsigned char edge_t;
typedef unsigned char node_t;
typedef unsigned dist_t;
typedef double flow_t;
typedef signed int gene_t;

const unsigned msb = (sizeof(flag_t)*8-1);


/* on constant memory */
flow_t* const c = (flow_t*)malloc(sizeof(flow_t)*N*N);        //Capacity
flow_t* const D = (flow_t*)malloc(sizeof(flow_t)*N*N);        //Demand
node_t* const e = (node_t*)malloc(sizeof(node_t)*E*2);        //Edgelist

#define RANDOM random
#define SRANDOM srandom

#include "include/print.h"
#include "include/timer.h"
#include "include/bitop.h"

Timer tm;

Timer eval, gene;

#include "include/gputools.h"
#include "include/ecmp.h"
#include "include/ga.h"






int main(){
	
    
    /* initialize constants */
    memset(c, 0, sizeof(flow_t)*N*N);
    memset(D, 0, sizeof(flag_t)*N*N);
    memset(e, 0, sizeof(flow_t)*E*2);
    
    /* set constant parameters from hedder file (c, D, e, and w) */
    for (node_t i = 0; i < N; i++) {
        for (node_t j = 0; j < N; j++) {
            c[i*N+j] = (flow_t)capacity[i][j];
            D[i*N+j] = (flow_t)demand[i][j];
        }
    }
    
    for (edge_t k = 0; k < E; k++) {
        e[k] = (node_t)edge[k][0];   e[k+E] = (node_t)edge[k][1];
    }
    
	
	
    /* start */
	tm.start();
	
    SRANDOM(time(NULL));
	//for(int i=0 ; i<5 ; i++){ printf("rand = %u\n", (unsigned)RANDOM()%N); }
	
    checkGaParameters();
	
	
    //variables
    double L[A+P] = {};
    dist_t w[A+P][E] = {};
	
	double L_best = 0.0;
	unsigned long G_best = 0;
	double sec_best = 0.0;
	
	
    //init
	geneInitialize(w, L);
	
	
	unsigned long generation = 0;
	
//#pragma omp parallel
	for (generation = 0; generation < GENMAX; generation++) {
		
		
		//evaluation for B and C
		eval.start();
#pragma omp parallel for num_threads(MP)
		for (gene_t g = A; g < A+P; g++) { calculateFlowFromWeight(L[g], w[g]); }
		eval.stop();
		
		gene.start();
		geneEvolution(w, L);
		gene.stop();
		
		
		double sec = tm.now();

		if(L_best != L[0]){
			L_best = L[0];  G_best = generation;  sec_best = sec;
#			ifdef PRINT
			//printf("%.3f\t%10f\n", sec, L[0]);
			printf("%7lu\t%10f\n", (generation+1), L[0]);
#			endif
		}
		
#	ifdef TMLIM
		if(sec > TMLIM){ printf("limit time elapsed (%d)\n", TMLIM); break; }
#	endif
		
#	ifdef GOAL
		if(L[0] <= GOAL){ printf("enough optimized for goal (%f)\n", GOAL); break;}
#	endif
		
		
	}
	
	
#ifdef TIME_MEASUREMENT
	printf("\n");
	printf("----statistics----\n");
	
	showTimeStatisticsEcmp();
	showTimeStatisticsEvolution();
	
	double total = eval.sec() + gene.sec();
	printf("  -total\n");
	printf("eval    : \t%.2f\t(%.2f percent)\n", eval.sec(), eval.sec()/total*100);
	printf("gene    : \t%.2f\t(%.2f percent)\n", gene.sec(), selec.sec()/total*100);


	
#endif

	printf("\n");
	printf("BEST_W:\t");	printGene(w[0]);
	printf("\n");
	printf("BEST_L:     \t%7f\n", L[0]);
	printf("BTIME(sec): \t%.3f\n", sec_best);
	printf("TIME(sec):  \t%.3f\n", tm.now());
	printf("GBEST:      \t%lu\n", G_best);
	printf("GENERATION: \t%lu\n", generation);
	
	
}











