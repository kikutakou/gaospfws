#include <iostream>
#include <iterator>

#ifndef N
#define N 128
#endif


void printArray(int a[N]){
    for (int i = 0; i < N; i++) { printf("%d ", a[i]); }
    printf("\n");
}
void swap(int& i, int& j){
    int tmp = i;    i = j;      j = tmp;
}


void bitonic(int a[N]){
    
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




int main ( int argc, char** argv)
{
    
    
    if( ( (N-1) & N ) ){ abort(); }
    
    int a[N];
    
    //generate random array(with no overlap)
    for (int i = 0; i < N; i++) { a[i] = i+1; }
    srand(time(NULL));
    for (int i = 0; i < 3*N; i++) { swap(a[rand()%N], a[rand()%N]); }
    
    
    
    printArray(a);
    
    bitonic(a);
    
    
    printArray(a);
    
    
	return 0;
};


