

#include <stdio.h>
#include <stdlib.h>

typedef unsigned int uint32_t;

uint32_t xor128(void) {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    static uint32_t w = 88675123;
    uint32_t t;
    
    t = x ^ (x << 11);
    x = y; y = z; z = w;
    return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
}


int main(){
    if(sizeof(uint32_t) != 4){ abort(); }
    
    while (true) {
        printf("xor128 %u\n", xor128() % 20);
    }
    
}


