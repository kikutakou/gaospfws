
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>

#ifndef N
#include "w15_10.h"
#endif



#include "cudaemu.h"
#include "cudacommon.h"

#define REPEAT 1

#define CONST


typedef unsigned long long flag_t;
typedef unsigned char edge_t;
typedef unsigned char node_t;
typedef unsigned dist_t;
typedef float flow_t;


#ifdef CONST
__constant__ node_t e[2*E];           //Edge
__constant__ flow_t C[N*N];           //Capacity
__constant__ flow_t D[N*N];           //Demand
__constant__ dist_t w[N*N];           //Weight
#endif

node_t* e_host;             //Edge
flow_t* C_host;             //Capacity
flow_t* D_host;             //Demand
dist_t* w_host;             //Weight


#ifndef __NDEBUG__



inline void printBitStream(flag_t b){
    bool first = 1;
    for(edge_t k = 0; k < E; k++){
        if( ( b >> k ) & 1 ){ if(first){ printf("%d", k); first = 0; }else{ printf(", %d", k); } }
    }
}

inline void printTable(flag_t* t){
    for (node_t i = 0; i < N; i++) {
        for (node_t j = 0; j < N; j++) {
            printf(" ["); printBitStream( t[i*N+j] ); printf("] ");
        }
        printf("\n");
    }
}


#endif




__global__ void gpu_ecmpWarshallfloydEdge(edge_t k, node_t x, node_t y, dist_t* d, flag_t* p){
    
    node_t i = threadIdx.x + blockIdx.x * blockDim.x;
    node_t j = threadIdx.y + blockIdx.y * blockDim.y;
    
    if(i<N && j<N){

        dist_t sub = d[i*N+x] + w[k] + d[y*N+j];
        if( sub < d[i*N+j] ) {
            d[i*N+j] = sub;
            p[i*N+j] = ((flag_t)1 << k);
            //printf("upd %d - %d via %d [ ", i, j, k); printBitStream(p[i*N+j]); printf("]\n");
        }else if( sub == d[i*N+j] ) {
            p[i*N+j] |= ((flag_t)1 << k);
            //printf("add %d - %d via %d [ ", i, j, k); printBitStream(p[i*N+j]); printf("]\n");
        }
    }
}






////__device__
flag_t getNeighborEdges(flag_t sub, const node_t s, const flag_t* p, const node_t* e){
    flag_t tmp = 0;
    for(edge_t k = 0; k < E; k++ ){ if( ( sub >> k ) & 1 ){       //only works where flag_t up
        const node_t near = e[k*2];
        //printf("near = %d\n", near);
        if(near == s){
            //printf("this is source node\n");
            tmp |= ((flag_t)1 << k);
        }else{
            //printf("follow again from %d to %d\n", s, near);
            tmp |= getNeighborEdges(p[s*N+near], s, p, e);
        }
    } }
    return tmp;
}


////__device__
void traceAndAlloc(const node_t h, const node_t d, const flag_t* n, int rate, const flow_t traffic, flow_t* f, node_t* e){
    
    if (h != d) {
        
        char branch = 0;
        for(edge_t k = 0; k < E; k++ ){
            if( ( n[h*N+d] >> k ) & 1 ){ branch++; }
        }
        rate *= branch;
        
        for(edge_t k = 0; k < E; k++ ){
            if( ( n[h*N+d] >> k ) & 1 ){
                node_t next = e[k*2+1];
                f[k] += traffic / rate;
                traceAndAlloc(next, d, n, rate, traffic, f, e);
            }
        }
        
    }
    
}




int main(){
    
    
    printf("sizeof = node_t:%lu, flow_t:%lu, dist_t:%lu, flag_t:%lu\n", sizeof(node_t), sizeof(flow_t), sizeof(dist_t), sizeof(flag_t));

    //property
    cudaDeviceProp  prop;    int whichDevice;
    HANDLE_ERROR( cudaGetDevice( &whichDevice ) );
    HANDLE_ERROR( cudaGetDeviceProperties( &prop, whichDevice ) );
    int mp = prop.multiProcessorCount;
    
    
    //time
    cudaEvent_t start, wf;     float ms=0.0f;
    HANDLE_ERROR( cudaEventCreate( &start ) );
    HANDLE_ERROR( cudaEventCreate( &wf ) );
    HANDLE_ERROR( cudaEventRecord( start, 0 ) );	//start
    
    
    ////----------HOST MEMORY ALLOC----------
    
    /* const variable */
    HANDLE_ERROR( cudaHostAlloc( (void**)&e_host, sizeof(node_t)*2*E, cudaHostAllocDefault ) );        //Edge
    HANDLE_ERROR( cudaHostAlloc( (void**)&C_host, sizeof(flow_t)*N*N, cudaHostAllocDefault ) );        //Capacity
    HANDLE_ERROR( cudaHostAlloc( (void**)&D_host, sizeof(flow_t)*N*N, cudaHostAllocDefault ) );        //Demand
    
    /* given variable */
    HANDLE_ERROR( cudaHostAlloc( (void**)&w_host, sizeof(dist_t)*N*N, cudaHostAllocDefault ) );         //Weight
    
    /* temporary variable */
    dist_t* d_host; HANDLE_ERROR( cudaHostAlloc( (void**)&d_host, sizeof(dist_t)*N*N, cudaHostAllocDefault ) );      //Distance
    flag_t* p_host; HANDLE_ERROR( cudaHostAlloc( (void**)&p_host, sizeof(flag_t)*N*N, cudaHostAllocDefault ) );      //Predtable
    flag_t* n_host; HANDLE_ERROR( cudaHostAlloc( (void**)&n_host, sizeof(flag_t)*N*N, cudaHostAllocDefault ) );      //Nexttable
    flow_t* f_host; HANDLE_ERROR( cudaHostAlloc( (void**)&f_host, sizeof(flow_t)*E, cudaHostAllocDefault ) );        //Flow
    
    //init on host
    for (node_t k = 0; k < E; k++) { e_host[k*2] = (node_t)edge[k][0];   e_host[k*2+1] = (node_t)edge[k][1];  }
    for (node_t i = 0; i < N; i++) { for (node_t j = 0; j < N; j++) { C_host[i*N+j] = (flow_t)capacity[i][j]; } }
    for (node_t i = 0; i < N; i++) { for (node_t j = 0; j < N; j++) { D_host[i*N+j] = (flow_t)distance[i][j]; } }
    for (node_t i = 0; i < N; i++) { for (node_t j = 0; j < N; j++) { w_host[i*N+j] = weightEdge[i*N+j]; } }
    for (node_t i = 0; i < N; i++) { for (node_t j = 0; j < N; j++) { d_host[i*N+j] = i==j ? 0 : 10000000; } }
    memset(p_host, 0, sizeof(flag_t)*N*N);
    memset(n_host, 0, sizeof(flag_t)*N*N);
    memset(f_host, 0, sizeof(flow_t)*E);
    
    //print
    //for (int i = 0; i < N; i++) { for (int j = 0; j < N; j++) { printf(" %d", d[i][j]); } printf("\n"); }
    //for (int i = 0; i < N; i++) { for (int j = 0; j < N; j++) { printf(" %u", p[i][j]); } printf("\n"); }
    //for (int i = 0; i < N; i++) { for (int j = 0; j < N; j++) { printf(" %d", ne[i][j]); } printf("\n"); }
    
    
    
    
    ////----------DEVICE MEMORY ALLOC and TRANSFER----------
#ifndef CONST
    node_t* e; HANDLE_ERROR( cudaMalloc( (void**)&e, sizeof(node_t)*2*E ) );        //Edge      (1x2xE)     120
    flow_t* C; HANDLE_ERROR( cudaMalloc( (void**)&C, sizeof(flow_t)*N*N ) );        //Capacity  (4xNxN)     900
    flow_t* D; HANDLE_ERROR( cudaMalloc( (void**)&D, sizeof(flow_t)*N*N ) );        //Demand    (4xNxN)     900
    dist_t* w; HANDLE_ERROR( cudaMalloc( (void**)&w, sizeof(dist_t)*N*N ) );        //Weight    (4xNxN)     900
#endif
    dist_t* d; HANDLE_ERROR( cudaMalloc( (void**)&d, sizeof(dist_t)*N*N ) );        //Distance  (4xNxN)     900
    flag_t* p; HANDLE_ERROR( cudaMalloc( (void**)&p, sizeof(flag_t)*N*N ) );        //Predtable (8xNxN)     1800
    flag_t* n; HANDLE_ERROR( cudaMalloc( (void**)&n, sizeof(flag_t)*N*N ) );        //Nexttable (8xNxN)     1800
    flow_t* f; HANDLE_ERROR( cudaMalloc( (void**)&f, sizeof(flow_t)*E ) );          //Flow      (4xE)       240
    
    
    
#ifndef CONST
    HANDLE_ERROR( cudaMemcpy( e, e_host, sizeof(node_t)*2*E, cudaMemcpyHostToDevice ) );
    HANDLE_ERROR( cudaMemcpy( C, C_host, sizeof(flow_t)*N*N, cudaMemcpyHostToDevice ) );
    HANDLE_ERROR( cudaMemcpy( D, D_host, sizeof(flow_t)*N*N, cudaMemcpyHostToDevice ) );
    HANDLE_ERROR( cudaMemcpy( w, w_host, sizeof(dist_t)*N*N, cudaMemcpyHostToDevice ) );
#else
    HANDLE_ERROR( cudaMemcpyToSymbol( e, e_host, sizeof(node_t)*2*E ) );
    HANDLE_ERROR( cudaMemcpyToSymbol( C, C_host, sizeof(flow_t)*N*N ) );
    HANDLE_ERROR( cudaMemcpyToSymbol( D, D_host, sizeof(flow_t)*N*N ) );
    HANDLE_ERROR( cudaMemcpyToSymbol( w, w_host, sizeof(dist_t)*N*N ) );
#endif
    HANDLE_ERROR( cudaMemcpy( d, d_host, sizeof(dist_t)*N*N, cudaMemcpyHostToDevice ) );
    HANDLE_ERROR( cudaMemset( p, 0, sizeof(flag_t)*N*N ) );
    HANDLE_ERROR( cudaMemset( n, 0, sizeof(flag_t)*N*N ) );
    HANDLE_ERROR( cudaMemset( f, 0, sizeof(flow_t)*E ) );
    
    
    ////----------GET EDGE_PRED----------

    dim3 blocks(mp,1);
    dim3 threads(N/blocks.x+1, N);
    printf("block/thread = x(%d,%d) y(%d,%d)\n", blocks.x, threads.x, blocks.y, threads.y );
    
    for (edge_t k = 0; k < E; k++) {       //for each edge_t
        node_t x = e_host[k*2],  y = e_host[k*2+1];
#       ifdef __CUDACC__
        gpu_ecmpWarshallfloydEdge<<<blocks, threads>>>(k, x, y, d, p);
#       else
        gpuemulate(blocks,threads) gpu_ecmpWarshallfloydEdge(k, x, y, d_host, p_host);
#       endif
        
    }

    
    
#   if 1
    HANDLE_ERROR( cudaEventRecord( wf, 0 ) );
    HANDLE_ERROR( cudaEventSynchronize( wf ) );
    HANDLE_ERROR( cudaEventElapsedTime( &ms, start, wf ) );
    printf( "time: %f ms\n", ms );
#   ifdef __CUDACC__
    HANDLE_ERROR( cudaMemcpy( p_host, p, sizeof(flag_t)*N*N, cudaMemcpyDeviceToHost ) );
#   endif
    //printf("Pred Table : \n");      printTable(p_host);
#   endif
    

    
    return 0;
    
    
    
    
    
    
    
    
    
    
    
    
    ////----------GET EDGE_NEXT----------
    for (node_t i = 0; i < N; i++) {
        for (node_t j = 0; j < N; j++) {
            if (i == j) { continue; }
            n_host[i*N+j] = getNeighborEdges(p_host[i*N+j], i, p_host, e);
        }
    }
    
#   if 0
    printf("Next Table : \n");      printTable(n_host);
    HANDLE_ERROR( cudaEventRecord( wf, 0 ) );
    HANDLE_ERROR( cudaEventElapsedTime( &ms, start, wf ) );
    printf( "time: %f ms\n", ms );
#   endif
    
    
    ////----------TRACE AND ALLOC----------
    for(edge_t k = 0; k < E; k++ ){ f[k] = 0; }
    
    for (node_t i = 0; i < N; i++) {
        for (node_t j = 0; j < N; j++) {
            
            if (i == j) { continue; }
            traceAndAlloc(i, j, n_host, 1, D_host[i*N+j], f_host, e_host);
            
        }
    }
    
    
    ////----------GET L----------
    
    double l = 0.0;
    
    for(int k = 0; k < E; k++ ){
        
        int x = e_host[k*2],  y = e_host[k*2+1];
        double tmp = f_host[k] / capacity[x][y];
        
        if(l < tmp){ l = tmp; }
    }
    
    
    printf("L = %10f\n", l);
    
    
    
#   if 0
    HANDLE_ERROR( cudaEventRecord( stop, 0 ) );	//stop
    HANDLE_ERROR( cudaEventSynchronize( stop ) );	//stop
    HANDLE_ERROR( cudaEventElapsedTime( &ms, start, stop ) );
    printf( "time: %f ms\n", ms );
#   endif
    
    
    // device mem free
#ifndef CONST
    HANDLE_ERROR( cudaFree( e ) );
    HANDLE_ERROR( cudaFree( C ) );
    HANDLE_ERROR( cudaFree( D ) );
    HANDLE_ERROR( cudaFree( w ) );
#endif
    HANDLE_ERROR( cudaFree( d ) );
    HANDLE_ERROR( cudaFree( p ) );
    HANDLE_ERROR( cudaFree( n ) );
    HANDLE_ERROR( cudaFree( f ) );
    
    // host mem free
    HANDLE_ERROR( cudaFreeHost( e_host ) );
    HANDLE_ERROR( cudaFreeHost( C_host ) );
    HANDLE_ERROR( cudaFreeHost( D_host ) );
    HANDLE_ERROR( cudaFreeHost( w_host ) );
    HANDLE_ERROR( cudaFreeHost( d_host ) );
    HANDLE_ERROR( cudaFreeHost( p_host ) );
    HANDLE_ERROR( cudaFreeHost( n_host ) );
    HANDLE_ERROR( cudaFreeHost( f_host ) );
    
    
    
    //time
    HANDLE_ERROR( cudaEventDestroy( start ) );
    HANDLE_ERROR( cudaEventDestroy( wf ) );
    
    
    
    
    
}











