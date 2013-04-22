#ifndef GA_H
#define GA_H


#ifndef MP
#define MP 2
#endif
#ifndef A
#define A 25
#endif
#ifndef B
#define B 300
#endif
#ifndef C
#define C 25
#endif

#define P (B+C)

#ifndef M
#define M 0.001
#endif
#ifndef K
#define K 0.8
#endif
#ifndef WMAX
#define WMAX 20
#endif
#ifndef GENMAX
#define GENMAX 1000
#endif
#ifndef TMLIM
#define TMLIM 3600
#endif

#define INF 10000000
#define INITIAL_RATE 1

typedef unsigned long long flag_t;
typedef unsigned short weit_t;
typedef unsigned dist_t;

typedef struct {
  double L;			//objective value
//	int r;				//ranking
} data_t;

typedef struct {
	flag_t visited[N];
	dist_t d[N*N];
	flag_t p[N*N];
	double f[N*N];
} locl_t;

#endif
