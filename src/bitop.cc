
#include "../include/bitop.h"

unsigned bitcount(unsigned char n) {
	if(sizeof(unsigned char) != 8/8) { abort(); }
    n = (n & 0x55) + ((n >> 1) & 0x55);
    n = (n & 0x33) + ((n >> 2) & 0x33);
    n = (n & 0x0f) + ((n >> 4) & 0x0f);
	return (unsigned)n;
}



unsigned bitcount(unsigned short n) {
	if(sizeof(unsigned short) != 16/8) { abort(); }
    n = (n & 0x5555) + ((n >> 1) & 0x5555);
    n = (n & 0x3333) + ((n >> 2) & 0x3333);
    n = (n & 0x0f0f) + ((n >> 4) & 0x0f0f);
    n = (n & 0x00ff) + ((n >> 8) & 0x00ff);
	return (unsigned)n;
}

unsigned bitcount(unsigned n) {
	if (sizeof(unsigned) != 32/8) { abort(); }
    n = (n & 0x55555555) + ((n >> 1) & 0x55555555);
    n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
    n = (n & 0x0f0f0f0f) + ((n >> 4) & 0x0f0f0f0f);
    n = (n & 0x00ff00ff) + ((n >> 8) & 0x00ff00ff);
    n = (n & 0x0000ffff) + ((n >> 16) & 0x0000ffff);
	return (unsigned)n;
}

unsigned bitcount(unsigned long n) {
	if(sizeof(unsigned long) != 64/8) { abort(); }
	n = (n & 0x5555555555555555) + ((n >> 1) & 0x5555555555555555);
    n = (n & 0x3333333333333333) + ((n >> 2) & 0x3333333333333333);
    n = (n & 0x0f0f0f0f0f0f0f0f) + ((n >> 4) & 0x0f0f0f0f0f0f0f0f);
    n = (n & 0x00ff00ff00ff00ff) + ((n >> 8) & 0x00ff00ff00ff00ff);
    n = (n & 0x0000ffff0000ffff) + ((n >> 16) & 0x0000ffff0000ffff);
    n = ((n & 0x00000000ffffffff) + ((n >> 32) & 0x00000000ffffffff));
	return (unsigned)n;
}


