#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "gaospfws.h"

const int (*cap_p)[N];		//capacity
const int (*dem_p)[N];		//demand
const int (*edg_p)[2];		//edge

static data_t g[ A + P ];
static data_t *gp[ A + P ];					//pointer array for sorting
static weit_t cb[ E * B ];					//data array for mating
static locl_t lo;

const double th = ( ( 1.0 - M ) * K * RAND_MAX );	//thresholds for cutoff K in mating
const int mar = (int)( M * RAND_MAX );						//margin for Mutation in mating

void dp_init( data_t *g, locl_t *lp );
void ecmpwarshallfroyd( data_t *g, locl_t *lp );
inline unsigned int bit_count(flag_t n);
void alloc_flow(data_t *g, locl_t *lp, double flow, unsigned int rate, int i, int j);
void buildnexthoplist( data_t *g, locl_t *lp, int i, int j );
void c_flow( data_t *g, locl_t *lp );
void mating( weit_t *c, const weit_t *el, const weit_t *nel );
void geneSorting( data_t **gp );
void geneEvolution( data_t **gp, weit_t *c );

void dp_init( data_t *g, locl_t *lp ) {
  int i,j;
  int ij;
  flag_t msb;
	
  for( i = 0; i < N; i++ ) {
    lp->visited[i] = 0;
    for( j = 0; j < N; j++ ) {
      lp->d[i*N+j] = 0;
      lp->p[i*N+j] = 0;
      lp->f[i*N+j] = 0;
    }
  }
	
  for( i = 0; i < N; i++ ) {
    for( j = 0; j < N; j++ ) {
      if (i != j){ lp->d[i*N+j] = INF; }   //large number if no link
    }
  }
	
  msb = (flag_t)1 << (sizeof(flag_t)*8 - 1);	//most significant bit
	
  for ( i = 0; i < E; i++){
    ij = edg_p[i][0]*N + edg_p[i][1];
    lp->d[ij] = g->w[i];
    lp->p[ij] = msb; //MSB if direct link
  }
	
  return;
}

void ecmpwarshallfroyd( data_t *g, locl_t *lp ) {
  int i,j,k;
  dist_t *d = lp->d;
  flag_t *p = lp->p;
  dist_t c;
	
  for (k = 0; k < N; k++) {		//for each transit node k
    for(i = 0; i < N; i++) {	//for each src and dest pair i,j
      if( i != k ) {
        for (j = 0; j < N; j++) {
          if( j != k && j != i ) {
            c = d[i*N+k] + d[k*N+j];	//candidate route
            if( c < d[i*N+j] ) {
              d[i*N+j] = c;									//update
              p[i*N+j] = ((flag_t)1 << k);
            } else if( c == d[i*N+j] ) {
              p[i*N+j] |= ((flag_t)1 << k);	//add
            }
          }
        }
      }
    }
  }
}

void buildnexthoplist( data_t *g, locl_t *lp, int i, int j ) {
  flag_t *visited = lp->visited;
  dist_t *d = lp->d;
  flag_t *p = lp->p;
  flag_t tmp;
  int k;
	
  if( ( i != j ) && ( ( ( visited[i] >> j ) & 1 ) == 0 ) ) {	//return if visited
    tmp = p[i*N+j] & ~( (flag_t)1 << (sizeof(flag_t)*8-1) );
		
    k = 0;
    while( tmp > 0 ){
      while( (tmp & 1) == 0 ) { k++; tmp >>= 1; }		//get next [k]th bit
      tmp &= ~(flag_t)1;							//frag down corresponding bit
      buildnexthoplist(g, lp, i, k);
      p[i*N+j] &= ~((flag_t)1 << k);	//frag down corresponding bit on p[i][j]
      p[i*N+j] |= p[i*N+k];						//merge the p[i][j] of closer node
    }
		
    if( ( p[i*N+j] >> ( sizeof(flag_t)*8 - 1 ) ) & 1 ) {			//if direct link
      p[i*N+j] &= ~((flag_t)1 << ( sizeof(flag_t)*8 - 1 ) ); //frag down msb
      p[i*N+j] |= ((flag_t)1 << j);			     //temporal memory
    }
		
    visited[i] |= ((flag_t)1 << j);		//frag up
  }
}

inline unsigned int bit_count(flag_t n) {
  n = (n & 0x5555555555555555) + ((n >> 1) & 0x5555555555555555);
  n = (n & 0x3333333333333333) + ((n >> 2) & 0x3333333333333333);
  n = (n & 0x0f0f0f0f0f0f0f0f) + ((n >> 4) & 0x0f0f0f0f0f0f0f0f);
  n = (n & 0x00ff00ff00ff00ff) + ((n >> 8) & 0x00ff00ff00ff00ff);
  n = (n & 0x0000ffff0000ffff) + ((n >> 16) & 0x0000ffff0000ffff);
  n = (n & 0x00000000ffffffff) + ((n >> 32) & 0x00000000ffffffff);
  return (unsigned int)n;
}

void alloc_flow(data_t *g, locl_t *lp, double flow, unsigned int rate, int i, int j) {
	
  int k = 0;
	
  flag_t tmp = lp->p[i*N+j];
	
  rate *= bit_count(tmp);
  
  while(tmp > 0){
    while((tmp & 1) == 0) { k++; tmp >>= 1; }
    tmp &= ~(flag_t)1;
		
    lp->f[i*N+k] += flow / rate;
    if (k != j) {
      alloc_flow(g, lp, flow, rate, k, j);
    }
  }
}

void c_flow( data_t *g, locl_t *lp ) {
  int i;
  int j;
	
  dp_init( g, lp );
  ecmpwarshallfroyd( g, lp );
	
  for( i = 0; i < N; i++ ) {
    for( j = 0; j < N; j++ ) {
      buildnexthoplist( g, lp, i, j );
    }
  }
	
  for ( i = 0; i < N; i++) {
    for ( j = 0; j < N; j++) {
      if( ( i != j ) && dem_p[i][j] ) {
        alloc_flow(g, lp, dem_p[i][j], INITIAL_RATE, i, j);
      }
    }
  }
	
  g->L = 0.0;
  for ( i = 0; i < E; i++){
    j = edg_p[i][0]*N + edg_p[i][1];
    lp->f[j] /= cap_p[ edg_p[i][0] ][ edg_p[i][1] ];
    if( g->L < lp->f[j] ){ g->L = lp->f[j]; }
  }
}

inline void mating( weit_t *c, const weit_t *el, const weit_t *nel ) {
  int i;
  double ran;
	
  for (i = 0; i < E; i++) {
    ran = random();
    ran -= mar;
    
    if( ran < 0.0 ){
      c[i] = ( random() % WMAX ) + 1;
    } else if( ran < th ){
      c[i] = el[i];
    } else {
      c[i] = nel[i];
    }
  }
}

void geneSorting( data_t **gp ){
  int i,j;
	double minL;
  int idx;
  data_t *gp_tmp;
	
	for( i = 0; i < A+P; i++ ) {		//sort
    minL = gp[i]->L; idx = i;
    for( j = i + 1; j < A+P; j++ ) {
      if( minL > gp[j]->L ) {
        minL = gp[j]->L;
        idx = j;
      }
    }
    gp_tmp = gp[i];
    gp[i] = gp[idx];
    gp[idx] = gp_tmp;
  }
	for( i = 0; i < A+P; i++ ) {
		gp[i]->r=i;
	}
}

void geneEvolution( data_t **gp, weit_t *cb ) {
  int i,j;
  data_t *elite, *nonelite;
  weit_t *c_p;
	
  for( i = A; i < A+B; i++ ) {				//elite => keep the same	,
    elite = gp[random() % A ];				//elite
    nonelite = gp[A + random() % P ];	//nonelite
    mating( &cb[ (i-A) * E ], elite->w, nonelite->w );	//mating => ( nonelite + elite )
  }
	
  for( i = A; i < A+B; i++ ) {	//update it back
    for( j = 0; j < E; j++ ) {
      gp[i]->w[j] = cb[ (i-A) * E + j ];			//copy
    }
  }

  for( i = A+B; i < A+P; i++ ) {	  //replace it with random one;
    for( j = 0; j < E; j++ ) {
      gp[i]->w[j] = ( random() % WMAX ) + 1;
      gp[i]->L = 100.0;
    }
  }
}

void dummy(data_t **gp, weit_t *c){
	int i, j;
	int a[P];
	int min, idx, tmp;
	for( i = 0; i < P; i++ ) {
		a[i] = random();
	}
	for( i = 0; i < P; i++ ) {		//sort
    min = a[i]; idx = i;
    for( j = i + 1; j < P; j++ ) {
      if( min > a[j] ) {
        min = a[j];
        idx = j;
      }
    }
    tmp = a[i];
    a[i] = a[idx];
    a[idx] = tmp;
  }
	for( i = A; i < A+P; i++ ) { gp[i]->L = 100.0; }
}


int main(){
	printf("NODE:\t%u\tPATH:\t%u\tSEED:\t%u\tGENMAX:\t%u\n", N, H, S, GENMAX);
	printf("A:\t%3d\tB:\t%3d\tC:\t%3d\tMUT:\t%.3f\tCUT:\t%.3f\n", A, B, C, M, K);
#ifdef _OPENMP
	printf("OPENMP:\t%d\n", MP);
#endif
	printf("\nsec/gen \tL\n");
	
  int i,j;
  double L_best = 0.0;
  double sec_best = 0.0;
  int G_best = 0;
	static struct timeval st, cu;
	double tm, tm_best;
	
  srandom(0);
	
	gettimeofday(&st, NULL);
	
	cap_p = capacity;
	dem_p = demand;
	edg_p = edge;
	
  for( i = 0 ; i < P+A; i++ ) {		//init process
    g[i].L = 100.0;
    for ( j = 0; j < E; j++) {
      g[i].w[j] = ( random() % WMAX ) + 1;		//initialize
    }
    gp[i] = &(g[i]); // initialize pointer setting
  }
	
#pragma omp parallel num_threads(MP) private(i)
  for( i = 0; i < GENMAX; i++ ){
#pragma omp for private(lo)
    for( j = A; j < A+P; j++) {
      c_flow( gp[j], &(lo) );							//calculation flow
    }
#pragma omp single
		{
			geneSorting(gp);			//gene sort
			geneEvolution(gp, cb);		//gene evolution
			gettimeofday(&cu, NULL);
			tm = 1e-6*(cu.tv_usec-st.tv_usec)+(cu.tv_sec-st.tv_sec);
			if(L_best != gp[0]->L){
				L_best = gp[0]->L;  G_best = i;	tm_best=tm;
				printf("%7u\t%10f\n", (i+1), gp[0]->L);
			}
		}
  }
	
	printf("\n");
	printf("BEST_L:     \t%7f\n", gp[0]->L);
	printf("BTIME(sec): \t%.3f\n", tm_best);
	printf("TIME(sec):  \t%.3f\n", tm);
	printf("GBEST:      \t%d\n", G_best);
	printf("GENERATION: \t%d\n", GENMAX);
  return 0;
}
