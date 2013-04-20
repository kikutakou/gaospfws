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
#define POP (A+B+C)
#define ALPHA ((double)A/POP)
#define BETA ((double)A/POP)

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
typedef unsigned char edge_t;
typedef unsigned char node_t;
typedef unsigned dist_t;
typedef double flow_t;
typedef signed int gene_t;

//#define RANDOM random
//#define SRANDOM srandom

typedef struct {
  const int (*cap_p)[N];
  const int (*dem_p)[N];
  const int (*edg_p)[2];
  double L;
  dist_t w[E];
  
  flag_t visited[N];
  dist_t d[N*N];
  flag_t p[N*N];
  flow_t f[N*N];
} cell_t;

#endif
