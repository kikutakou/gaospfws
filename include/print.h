#ifndef PRINT_H
#define PRINT_H


void printFrags(flag_t b){
    bool first = 1;
    if( (b >> msb ) & 1 ){ if(first){ printf("-"); first = 0; }else{ printf(", -"); } }
    for(edge_t k = 0; k < N; k++){ if( ( b >> k ) & 1 ){ if(first){ printf("%d", k); first = 0; }else{ printf(", %d", k); } } }
}

void printFrag_t(flag_t* t){
    
    for (node_t i = 0; i < N; i++) {
        printf("   ");
        for (node_t j = 0; j < N; j++) {
            printf(" ["); printFrags( t[i*N+j] ); printf("] ");
        }
        printf("\n");
    }
}

void printFlow_t(flow_t* f){
    for (node_t i = 0; i < N; i++) {
        if(!i){printf("f = [");}else{printf("     ");}
        for (node_t j = 0; j < N; j++) {
            if(!j){ printf(" [ %.3f", f[i*N+j]);}else{ printf(", %.3f", f[i*N+j]); }
        }
        if(i<N-1){printf("], \n");}else{printf("] ]\n");}
    }
}

void printDist_t(dist_t* w){
    for (node_t i = 0; i < N; i++) {
        if(!i){printf("    [");}else{printf("     ");}
        for (node_t j = 0; j < N; j++) {
            if(!j){ if(w[i*N+j]>99999){printf(" [   INF");}else if(w[i*N+j]){printf(" [ %5d", w[i*N+j]);}else{printf(" [   nil");} }
            else{ if(w[i*N+j]>99999){printf(",   INF");}else if(w[i*N+j]){printf(", %5d", w[i*N+j]);}else{printf(",   nil");} }
        }
        if(i<N-1){printf("], \n");}else{printf("] ]\n");}
    }
}
void printWeight(dist_t* w){
    dist_t wn[N*N] = {};  for (edge_t k = 0; k < E; k++) { wn[e[k]*N+e[k+E]] = w[k]; } printDist_t(wn);
}

void printGene(dist_t* w){
    for (edge_t k = 0; k < E; k++) { if(!k){ printf(" w = [ %3d", w[k]); }else{ printf(", %3d", w[k]); } } printf("]\n");
}

void printEdge(){
    for (edge_t k = 0; k < E; k++) { printf(" edge(%u) : %2u - %2u\n", k, e[k], e[k+E]); }
}


#endif


