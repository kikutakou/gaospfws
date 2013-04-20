#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#include "gaospfws.h"

static data_t g[ A + P ];
static data_t *gp[ A + P ];					//pointer array
static cell_t ce[ P ];					//pointer array
static dist_t cb[ E * B ];

const double th = ( ( 1.0 - M ) * K * RAND_MAX );	//thresholds for cutoff K in mating
const int mar = (int)( M * RAND_MAX );						//margin for Mutation in mating

void dp_init( data_t *g, cell_t *cep );
void ecmpwarshallfroyd( data_t *g, cell_t *cep );
inline unsigned int bit_count(flag_t n);
void alloc_flow(data_t *g, cell_t *cep, double flow, unsigned int rate, int i, int j);
void buildnexthoplist( data_t *g, cell_t *cep, int i, int j );
void c_flow( data_t *g, cell_t *cep );
void mating( dist_t *c, dist_t *el, dist_t *nel );
void geneEvolution( data_t **gp, dist_t *c );

void dp_init( data_t *g, cell_t *cep ) {
  int i,j;
  int ij;
  flag_t msb;

  for( i = 0; i < N; i++ ) {
    cep->visited[i] = 0;
    for( j = 0; j < N; j++ ) {
      cep->d[i*N+j] = 0;
      cep->p[i*N+j] = 0;
      cep->f[i*N+j] = 0;
    }
  }

  for( i = 0; i < N; i++ ) {
    for( j = 0; j < N; j++ ) {
      if (i != j){ cep->d[i*N+j] = INF; }   //large number if no link
    }
  }

  msb = (flag_t)1 << (sizeof(flag_t)*8 - 1);	//most significant bit

  for ( i = 0; i < E; i++){		
    ij = g->edg_p[i][0]*N + g->edg_p[i][1];
    cep->d[ij] = g->w[i];
    cep->p[ij] = msb; //MSB if direct link
  }

  return;
}

void ecmpwarshallfroyd( data_t *g, cell_t *cep ) {
  int i,j,k;
  dist_t *d = cep->d;
  flag_t *p = cep->p;
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

void buildnexthoplist( data_t *g, cell_t *cep, int i, int j ) {
  flag_t *visited = cep->visited;
  dist_t *d = cep->d;
  flag_t *p = cep->p;
  flag_t tmp;
  int k;

  if( ( i != j ) && ( ( ( visited[i] >> j ) & 1 ) == 0 ) ) {	//return if visited
    tmp = p[i*N+j] & ~( (flag_t)1 << (sizeof(flag_t)*8-1) );

    k = 0;
    while( tmp > 0 ){
      while( (tmp & 1) == 0 ) { k++; tmp >>= 1; }		//get next [k]th bit
      tmp &= ~(flag_t)1;							//frag down corresponding bit
      buildnexthoplist(g, cep, i, k);
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

void alloc_flow(data_t *g, cell_t *cep, double flow, unsigned int rate, int i, int j) {

  int k = 0;
	
  flag_t tmp = cep->p[i*N+j];

  rate *= bit_count(tmp);
  
  while(tmp > 0){
    while((tmp & 1) == 0) { k++; tmp >>= 1; }
    tmp &= ~(flag_t)1;

    cep->f[i*N+k] += flow / rate;
    if (k != j) {
      alloc_flow(g, cep, flow, rate, k, j);
    }
  }
}

void c_flow( data_t *g, cell_t *cep ) {
  int i;
  int j;

  dp_init( g, cep );
  ecmpwarshallfroyd( g, cep );

  for( i = 0; i < N; i++ ) {
    for( j = 0; j < N; j++ ) {
      buildnexthoplist( g, cep, i, j );
    }
  }

  for ( i = 0; i < N; i++) {
    for ( j = 0; j < N; j++) {
      if( ( i != j ) && g->dem_p[i][j] ) {
        alloc_flow(g, cep, g->dem_p[i][j], INITIAL_RATE, i, j);
      }
    }
  }

  g->L = 0.0;
  for ( i = 0; i < E; i++){
    j = g->edg_p[i][0]*N + g->edg_p[i][1];
    cep->f[j] /= g->cap_p[ g->edg_p[i][0] ][ g->edg_p[i][1] ];
    if( g->L < cep->f[j] ){ g->L = cep->f[j]; }
  }
}

void mating( dist_t *c, dist_t *el, dist_t *nel ) {
  int i;
  double ran;

  for (i = 0; i < E; i++) {
    ran = rand();
    ran -= mar;
    
    if( ran < 0.0 ){
      c[i] = ( rand() % WMAX ) + 1;
    } else if( ran < th ){
      c[i] = el[i];
    } else {
      c[i] = nel[i];
    }
  }
}

void geneEvolution( data_t **gp, dist_t *cb ) {
  int i,j;
  data_t *gp_tmp, *gp_tmp2;
  int idx;
  dist_t *c_p;
  double minL;

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

  // elite => keep the same
  // mating => ( nonelite + elite )

  c_p = cb;

	
  for( i = A; i < A+B; i++ ) {
    gp_tmp = gp[rand() % A ];				//elite
    gp_tmp2 = gp[A + rand() % P ];	//nonelite
    mating( c_p, gp_tmp->w, gp_tmp2->w );
    c_p += E;
  }

  c_p = cb;

  for( i = A; i < A+B; i++ ) {	// update it back
    for( j = 0; j < E; j++ ) {
      gp[i]->w[j] = c_p[j];			//copy
    }
    c_p += E;
  }

  // replace it with rand one;

  for( i = A+B; i < A+P; i++ ) {
    for( j = 0; j < E; j++ ) {
      gp[i]->w[j] = ( rand() % WMAX ) + 1;
      gp[i]->L = 100.0;
    }
  }
}

int main() {
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

  srand(0);


	gettimeofday(&st, NULL);
	


	
  for( i = 0 ; i < P+A; i++ ) {
    g[i].cap_p = capacity;
    g[i].dem_p = demand;
    g[i].edg_p = edge;
    g[i].L = 100.0;
    for ( j = 0; j < E; j++) { 
      g[i].w[j] = ( rand() % WMAX ) + 1;		//initialize
    }
    gp[i] = &(g[i]); // initialize pointer setting
  }  

  for( i = 0; i < GENMAX; i++ ) {
#pragma omp parallel for num_threads(MP)
    for( j = A; j < A+P; j++) { 
      c_flow( gp[j], &(ce[j-A]) );							//calculation flow
    }
    geneEvolution(gp, cb);		//gene evolution

		gettimeofday(&cu, NULL);
		tm = 1e-6*(cu.tv_usec-st.tv_usec)+(cu.tv_sec-st.tv_sec);
    if(L_best != gp[0]->L){
      L_best = gp[0]->L;  G_best = i;	tm_best=tm;
      printf("%7u\t%10f\n", (i+1), gp[0]->L);
    }
  }
	
	printf("\n");
	printf("BEST_L:     \t%7f\n", gp[0]->L);
	printf("BTIME(sec): \t%.3f\n", tm_best);
	printf("TIME(sec):  \t%.3f\n", tm);
	printf("GBEST:      \t%d\n", G_best);
	printf("GENERATION: \t%d\n", i);
  return 0;
}
