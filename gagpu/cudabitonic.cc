#include "cudaemu.h"
#include "cudacommon.h"


#ifndef N
#   define N 512
#endif
#define REPEAT 1000

void printArray(int a[N]){
    for (int i = 0; i < N; i++) { printf("%d ", a[i]); }
    printf("\n");
}


void swap(int& i, int& j){
    int tmp = i;    i = j;      j = tmp;
}





__global__ void bitonic(int a[N]){
    
    for (int k = 2; k <= N; k <<= 1) {      //k : maximum size
        
        for (int j = k >> 1; j > 0 ; j = j>>1 ) {   //j : compare with
            
            //as long as j in power of 2 :   i^j = {i+j (if i>j),   i-j (if i<j)}
            
            for( int i = 0; i < N ; i++ ) { // (A)
                
                int o = i^j;
                
                if ( o > i ) {
                    
                    int ai = a[i];
                    int ao = a[o];
                    
                    if( ((i&k)>0)^(ai > ao) ){ a[i] = ao;  a[o] = ai; }
                }
            }
        }
    }
}



__global__ void bitonic_thread(int a[N]){
    
    for (int k = 2; k <= N; k <<= 1) {      //k : maximum size
        
        for (int j = k >> 1; j > 0 ; j = j>>1 ) {   //j : compare with
            
            
            int i = threadIdx.x;
            
            int o = i^j;
            
            if ( o > i ) {
                
                int ai = a[i];
                int ao = a[o];
                
                if( ((i&k)>0)^(ai > ao) ){ a[i] = ao;  a[o] = ai; }
            }

            __syncthreads();

            
        }
    }
}


__global__ void bitonic_shared(int a[N]){
    
    for (int k = 2; k <= N; k <<= 1) {      //k : maximum size
        
        for (int j = k >> 1; j > 0 ; j = j>>1 ) {   //j : compare with
            
            
            int i = threadIdx.x;
            
            int o = i^j;
            
            if ( o > i ) {
                
                int ai = a[i];
                int ao = a[o];
                
                if( ((i&k)>0)^(ai > ao) ){ a[i] = ao;  a[o] = ai; }
            }
            
            __syncthreads();
            
            
        }
    }
}



int main(int argc, char** argv){

    
    //parallel
    //int p; if (argc > 1) { p = atoi(argv[1]); } if(!p){ p = N; }
    
    
    if( ( (N-1) & N ) ){ abort(); }
    
    
    
    int a[N];
    
    //generate random array(with no overlap)
    for (int i = 0; i < N; i++) { a[i] = i+1; }
    srand(time(NULL));
    for (int i = 0; i < 3*N; i++) { swap(a[rand()%N], a[rand()%N]); }
    //printArray(a);
    
    
    
    
    int* a_dev;
    HANDLE_ERROR( cudaMalloc( (void**)&a_dev, sizeof(int)*N ) );
    HANDLE_ERROR( cudaMemcpy( a_dev, a, sizeof(int)*N, cudaMemcpyHostToDevice ) );
    
    int* b_dev;
    HANDLE_ERROR( cudaMalloc( (void**)&b_dev, sizeof(int)*N ) );
    
    
    //time
    cudaEvent_t start, stop;     float ms=0.0f;
    HANDLE_ERROR( cudaEventCreate( &start ) );
    HANDLE_ERROR( cudaEventCreate( &stop ) );
    HANDLE_ERROR( cudaEventRecord( start, 0 ) );	//start
    
    
    for (int r = 0; r < REPEAT; r++) {
        
        
#ifdef __CUDACC__
        //HANDLE_ERROR( cudaMemcpy( b_dev, a_dev, sizeof(int)*N, cudaMemcpyDeviceToDevice ) );
        //bitonic<<<1, 1>>>(b_dev);
        bitonic_thread<<<1, N>>>(b_dev);
        //bitonic_shared<<<1, N>>>(b_dev);
#else
        //HANDLE_ERROR( cudaMemcpy( b_dev, a_dev, sizeof(int)*N, cudaMemcpyDeviceToDevice ) );
        gpuemulate(1,1) bitonic(b_dev);
#endif
        
    }
        
    
    HANDLE_ERROR( cudaEventRecord( stop, 0 ) );	//stop
    HANDLE_ERROR( cudaEventSynchronize( stop ) );	//stop
    HANDLE_ERROR( cudaEventElapsedTime( &ms, start, stop ) );
    
    
    
    
#ifndef METHOD
    HANDLE_ERROR( cudaMemcpy( a, b_dev, sizeof(int)*N, cudaMemcpyDeviceToHost ) );
#endif
    
    
    
    HANDLE_ERROR( cudaFree( b_dev ) );
    HANDLE_ERROR( cudaFree( a_dev ) );
    
    //printArray(a);
    printf( "results :    N = %d arrays,  R = %d times  time = %f ms\n", N, REPEAT, ms );
    
    
	return 0;
};






















