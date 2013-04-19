#define N 8
#define E 32
#define H 10
#define S 0
const int weight[8][8] = {
 { 0, 0, 1, 10, 1, 0, 2, 2 }, 
 { 0, 0, 8, 999, 0, 1, 1, 0 }, 
 { 7, 6, 0, 9, 1, 7, 0, 1 }, 
 { 9, 1, 9, 0, 0, 3, 3, 0 }, 
 { 999, 0, 999, 0, 0, 0, 0, 0 }, 
 { 0, 9, 7, 9, 0, 0, 0, 0 }, 
 { 7, 1, 0, 11, 0, 0, 0, 1 }, 
 { 7, 0, 6, 0, 0, 0, 9, 0 }
};
const int weightEdge[32] = { 1, 10, 1, 2, 2, 8, 999, 1, 1, 7, 6, 9, 1, 7, 1, 9, 1, 9, 3, 3, 999, 999, 9, 7, 9, 7, 1, 11, 1, 7, 6, 9 };
const int capacity[8][8] = { 
 { 0, 0, 90, 87, 13, 0, 77, 64 }, 
 { 0, 0, 62, 70, 0, 94, 57, 0 }, 
 { 90, 62, 0, 59, 53, 59, 0, 45 }, 
 { 87, 70, 59, 0, 0, 90, 26, 0 }, 
 { 13, 0, 53, 0, 0, 0, 0, 0 }, 
 { 0, 94, 59, 90, 0, 0, 0, 0 }, 
 { 77, 57, 0, 26, 0, 0, 0, 44 }, 
 { 64, 0, 45, 0, 0, 0, 44, 0 }
};
const int demand[8][8] = { 
 { 0, 0, 1, 1, 0, 0, 0, 0 }, 
 { 0, 0, 0, 0, 0, 0, 0, 0 }, 
 { 0, 1, 0, 0, 1, 0, 0, 0 }, 
 { 0, 0, 0, 0, 1, 0, 1, 0 }, 
 { 0, 0, 0, 0, 0, 0, 0, 0 }, 
 { 0, 0, 0, 0, 0, 0, 0, 1 }, 
 { 0, 0, 0, 1, 0, 0, 0, 0 }, 
 { 0, 0, 0, 1, 0, 1, 0, 0 }
};
const int edge[32][2] = { { 0, 2 }, { 0, 3 }, { 0, 4 }, { 0, 6 }, { 0, 7 }, { 1, 2 }, { 1, 3 }, { 1, 5 }, { 1, 6 }, { 2, 0 }, { 2, 1 }, { 2, 3 }, { 2, 4 }, { 2, 5 }, { 2, 7 }, { 3, 0 }, { 3, 1 }, { 3, 2 }, { 3, 5 }, { 3, 6 }, { 4, 0 }, { 4, 2 }, { 5, 1 }, { 5, 2 }, { 5, 3 }, { 6, 0 }, { 6, 1 }, { 6, 3 }, { 6, 7 }, { 7, 0 }, { 7, 2 }, { 7, 6 } };
