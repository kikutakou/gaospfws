#ifndef GATOOLS_H
#define GATOOLS_H

typedef unsigned int uint32_t;

uint32_t xor128(void) {
    static uint32_t x = 123456789;
    static uint32_t y = 362436069;
    static uint32_t z = 521288629;
    //static uint32_t w = 88675123;
    static uint32_t w = (unsigned long)time(NULL);
    uint32_t t;
    
    t = x ^ (x << 11);
    x = y; y = z; z = w;
    return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
}

//#define RANDOM xor128



#endif


