
#ifndef GA_H
#define GA_H



void randomGene(dist_t w[E]){
    for (node_t k = 0; k < E; k++) { w[k] = RANDOM() % WMAX + 1; }
    //for (node_t k = 0; k < E; k++) { w[k] = (dist_t)weightEdge[k]; }
}

void mating(dist_t* child, dist_t* elite, dist_t* nonelite){
    for (edge_t k = 0; k < E; k++){
		
        double ran = (double)RANDOM() / RAND_MAX;
		ran -= M;
		
        if(ran < 0.0){
            child[k] = RANDOM() % WMAX + 1;			//printf("p ");
			
        }else if(ran < ((double)(1.0-M)*K)){
            child[k] = elite[k];					//printf("e ");
			
        }else{
            child[k] = nonelite[k];					//printf("n ");
        }
        
    }
    
}


void geneInitialize(dist_t w[A+P][E], double L[A+P]){
	
    for (gene_t g = 0; g < P+A; g++){ randomGene(w[g]); }
    for (gene_t g = 0; g < A; g++){ L[g] = 100.0; }
	
}




void geneRanking(gene_t rank[A+P], double L[A+P]){
    
	//initialize rank
    for (gene_t g = 0; g < P+A; g++){ rank[g] = g; }
	
	//ranking for all genes (bouble sort)
    for (gene_t i = 0; i < A+P; i++) {
        
        double minL = L[i];
        gene_t m = i;

        //find minimum
        for (gene_t j = i+1; j < A+P; j++) {
            if (L[j] < minL) { minL = L[j];  m = j; }
        }
        
		//exchange
        L[m] = L[i];			L[i] = minL;
        gene_t tmp = rank[m];	rank[m] = rank[i];  rank[i] = tmp;
        //printf(" rank : "); for (gene_t g = 0; g < A+P; g++) { printf("%u ", rank[g]); } printf("\n");
        
    }
    
}





void eliteSelection(dist_t w[P+A][E], gene_t rank[A+P]){
    
	static dist_t w_elite[A+B][E];

    //for elite	A		w[rank] -> w_elite
    for (gene_t g = 0; g < A; g++){
        memcpy(w_elite[g], w[rank[g]], sizeof(dist_t)*E);
        //printf("new gene #[%3u] : from Elite #[%3u]\n", g, rank[g]);
    }
    
    //nonelite B
    for (gene_t g = A; g < A+B; g++){

        gene_t elite = rank[RANDOM() % A];
        gene_t nonelite = rank[A + RANDOM() % P];
		
        //printf("new gene #[%3u] : from Father #[%3u] and Mother #[%3u]\n", g, father, mother);
        mating(w_elite[g], w[elite], w[nonelite]);
    }

    //copy elite to w
	memcpy(w, w_elite, sizeof(dist_t)*(A+B)*E);

	//replace C
	for (gene_t g = A+B; g < A+P; g++){
        //printf("new gene #[%3u] : RANDOM\n", g);
        randomGene(w[g]);
    }
	
}



#ifdef TIME_MEASUREMENT
Timer sort, selec, repro;


void showTimeStatisticsEvolution(){
	
	double total = sort.sec() + selec.sec();

	printf("  -regeneration\n");
	printf("sort    : \t%.2f\t(%.2f percent)\n", sort.sec(), sort.sec()/total*100);
	printf("selec   : \t%.2f\t(%.2f percent)\n", selec.sec(), selec.sec()/total*100);
	printf("total   : \t%.2f\n", total );
	
}

#endif






//---------------------- main function ----------------------


void checkGaParameters(){
	
	if(M > 1.0){ fprintf(stderr, "ERROR  : parameter wrong  M(%.1f) > 1.0\n", M); abort(); }
	if(K > 1.0){ fprintf(stderr, "ERROR  : parameter wrong  K(%.1f) > 1.0\n", K); abort(); }
	
#ifdef PRINT
	printf("NODE:\t%u\tPATH:\t%u\tSEED:\t%u\tGENMAX:\t%u\n", N, H, S, GENMAX);
	printf("A:\t%3d\tB:\t%3d\tC:\t%3d\tMUT:\t%.3f\tCUT:\t%.3f\n", A, B, C, M, K);
	printf("POP:\t%3d\tA(%.2f)\tB(%.2f)\tC(%.2f)\t\n", POP, (double)A/POP, (double)B/POP, (double)C/POP);
#ifdef _OPENMP
	printf("OPENMP:\t%d\n", MP);
#endif
	
	printf("\n");
	printf("sec/gen \tL\n");
#endif
	
#if DEBUG > 2
	printf("\nA(%d) : 0 - %d,  B(%d) : %d - %d,  C(%d) : %d - %d\n", A, A-1, B, A, A+B-1, C, A+B, A+B+C-1);
#endif
	
}






void geneEvolution(dist_t w[A+P][E], double L[A+P]){
	
    static gene_t rank[A+P];

	//printf("------------------------------------\n");
	//for (gene_t g = 0; g < A+P; g++) { printf("--- gene : %2u ( L = %10f )", g, L[g]);  printGene(w[g]); }
	
	//rank
	TIMER_START(sort);
	geneRanking(rank, L);     //sort for A+P
	TIMER_STOP(sort);
	
	//printf(" rank : "); for (gene_t g = 0; g < A+P; g++) { printf("%u ", rank[g]); } printf("\n");
	
	//selection
	TIMER_START(selec);
	eliteSelection(w, rank);	//renew w by rank
	TIMER_STOP(selec);
	
	//for (gene_t g = 0; g < A+P; g++) { printf("=== gene : %2u ( L = %10f )", g, L[g]);  printGene(w[g]); }
	
	
	
#   if DEBUG >= 2
	
	//printf("GENERATION[%2u]( L = %10f )", r, L[0]);  printGene(w[0]);  //printWeight(w[g]);
	//printStatistics(w, A);
#   endif
	
	
	
}







#endif


