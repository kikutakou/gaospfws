
#ifndef ECMP_H
#define ECMP_H

#include <assert.h>

const dist_t INF = (10000000);

const unsigned INITIAL_RATE = 1;


// initialize d[i][j] with w[k] or INF
void initialize(dist_t* d, flag_t* p, dist_t* w){
    
	
    for (node_t i = 0; i < N; i++) {
        for (node_t j = 0; j < N; j++) {
            if (i != j){ d[i*N+j] = INF; }   //large number if no link
        }
    }
    
    for (edge_t k = 0; k < E; k++){
        node_t ij = e[k]*N + e[k + E];
        d[ij] = w[k];
        p[ij] = (flag_t)1 << (sizeof(flag_t)*8-1);         //MSB if direct link
    }
    
}




void ecmpWarshallfloyd(dist_t* d, flag_t* p){
    
    for (node_t k = 0; k < N; k++) {       //for each transit node k
		
        for(node_t i = 0; i < N; i++) {   //for each src and dest pair i,j
            if(i == k){ continue; }
            for (node_t j = 0; j < N; j++) {
                if(j == i || j == k){ continue; }
                
                dist_t candidate = d[i*N+k] + d[k*N+j];	//candidate route
                
                if( candidate < d[i*N+j] ) {
                    d[i*N+j] = candidate;
                    p[i*N+j] = ((flag_t)1 << k);		//update
                    //printf("new %d - %d via %d [ ", i, j, k); printFrags(p[i*N+j]); printf("]\n");
					
                }else if( candidate == d[i*N+j] ) {
                    p[i*N+j] |= ((flag_t)1 << k);		//append
                    //printf("add %d - %d via %d [ ", i, j, k); printFrags(p[i*N+j]); printf("]\n");
                    
                }
            }
        }
    }
    
}


////__device__
void buildNexthopList(node_t i, node_t j, flag_t* p, flag_t* visited){
    
    if( visited[i] >> j & 1 ){ return; }		//return if visited
    //printf("buildNexthopList %hhu %hhu\n", i, j);

    flag_t tmp = p[i*N+j] & ~((flag_t)1 << msb);          //copy except msb

    node_t k = 0;
    while(tmp > 0){
        while((tmp & 1) == 0){ k++; tmp >>= 1; }      //get next [k]th bit
		assert(k < N);
        tmp &= ~(flag_t)1;		//frag down on correspond bit on tmp
		
        buildNexthopList(i, k, p, visited);
		
		p[i*N+j] &= ~((flag_t)1 << k);		//frag down correspond bit on p[i][j]
		p[i*N+j] |= p[i*N+k];				//merge nearer p[i][j]
    }
	
    if(p[i*N+j] >> msb & 1){				//msb for direct link
		
        p[i*N+j] &= ~((flag_t)1 << msb);         //(frag down msb)
        p[i*N+j] |= ((flag_t)1 << j);			//temporal memory
    }
	
    visited[i] |= ((flag_t)1 << j);		//frag up
    return;
}




////__device__
void allocFlow(node_t h, node_t d, flow_t flow, unsigned rate, flag_t* p, flow_t* f){

	flag_t tmp = p[h*N+d];          //copy
	
    rate *= bitcount(tmp);		//no of branch    
    //printf("rate %u branch %u\n", rate, branch);
    
    node_t k = 0;
    while(tmp > 0){
        while((tmp & 1) == 0){ k++; tmp >>= 1; }      //get next [k]th bit
		assert(k < N);
        tmp &= ~(flag_t)1;		//frag down on correspond bit on tmp
        
        //printf("%u - %u via %u : %f / %u\n", h, d, k, flow, rate);
        f[h*N+k] += flow / rate;
        
        if (k == d) { continue; }
		
        allocFlow(k, d, flow, rate, p, f);
		
    }
    
    
}










#ifdef TIME_MEASUREMENT
Timer init, ecmpwf, nextlist, trace, calcl;

void showTimeStatisticsEcmp(){

	double total = init.sec() + ecmpwf.sec() + nextlist.sec() + trace.sec() + calcl.sec();

	printf("  -evaluation\n");	
	printf("init    : \t%.2f\t(%.2f percent)\n", init.sec(), init.sec()/total*100);
	printf("ecmpwf  : \t%.2f\t(%.2f percent)\n", ecmpwf.sec(), ecmpwf.sec()/total*100);
	printf("next    : \t%.2f\t(%.2f percent)\n", nextlist.sec(), nextlist.sec()/total*100);
	printf("trace   : \t%.2f\t(%.2f percent)\n", trace.sec(), trace.sec()/total*100);
	printf("calcl   : \t%.2f\t(%.2f percent)\n", calcl.sec(), calcl.sec()/total*100);
	printf("total   : \t%.2f\n", total );
	
	
}

#endif




//---------------------- main function ----------------------



void calculateFlowFromWeight(double& L, dist_t w[E]){      // w -> L

	dist_t d[N*N];      //Distance
	flag_t p[N*N];      //Preceding
	flow_t f[N*N];      //Flow
    
    ////----------INIT d----------
	TIMER_START(init);
    initialize(d, p, w);			// w -> d,p
    TIMER_STOP(init);
	
#   if DEBUG == 1
    printf("Distance Table : d \n");      printDist_t(d);
    printf(" tm = %f msec \n", tm.now());
#   endif
    
    
    ////----------GET PRED----------
	TIMER_START(ecmpwf);
    ecmpWarshallfloyd(d, p);		// d -> d,p
	TIMER_STOP(ecmpwf);
	
#   if DEBUG == 1
    printf("Pred Table : p \n");      printFrag_t(p);
    printf(" tm = %f msec \n", tm.now());
#   endif
    
    
    
    
    ////----------GET EDGE_NEXT----------
	TIMER_START(nextlist);
    
	flag_t visited[N] = {};				//flag for visited node[i] of [j]th bit 
	
    for (node_t i = 0; i < N; i++) {
        for (node_t j = 0; j < N; j++) {
            if (i == j) { continue; }
            buildNexthopList(i, j, p, visited);		// p -> p
        }
    }
	TIMER_STOP(nextlist);
    
#   if DEBUG == 1
    printf("Next Table : p \n");      printFrag_t(p);
	printf(" tm = %f msec \n", tm.now());
#   endif
    
    
    
    
    
    ////----------TRACE AND ALLOC----------
	TIMER_START(trace);
    for (node_t i = 0; i < N; i++) {
        for (node_t j = 0; j < N; j++) {
            
            if (i == j) { continue; }
            if (!D[i*N+j]) { continue; } //only D
			
            allocFlow(i, j, D[i*N+j], INITIAL_RATE, p, f);  //p -> f
            
        }
    }
  	TIMER_STOP(trace);
    
#   if DEBUG == 1
    printf("Flow Table : f \n");      printFlow_t(f);
    printf(" tm = %f msec \n", tm.now());
#   endif
    
    
    ////----------CALC L----------
	TIMER_START(calcl);
	
    L = 0.0;
    for (edge_t k = 0; k < E; k++){
        node_t ij = e[k]*N+e[k+E];
        f[ij] /= c[ij];
        if(L < f[ij]){ L = f[ij]; }   //update L
    }
	TIMER_STOP(calcl);
    
#   if DEBUG == 1
    printf("Congestion Table : f \n");      printFlow_t(f);
    printf(" tm = %f msec \n", tm.now());
#   endif
    
	


#	if DEBUG > 2
	printf("--- gene : %2u ( L = %10f )", g, L[g]);  printGene(w[g]);  //printWeight(w[g]);
#	endif

	
	
    return;
}









#endif


