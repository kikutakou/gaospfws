#define N 7
#define E 28
#define H 10
#define S 0
const int weight[7][7] = {
 { 0, 0, 1, 3, 3, 1, 1 }, 
 { 0, 0, 1, 993, 995, 993, 2 }, 
 { 1000, 4, 0, 995, 999, 0, 1 }, 
 { 4, 1, 5, 0, 0, 0, 0 }, 
 { 997, 1, 1, 0, 0, 0, 0 }, 
 { 6, 1, 0, 0, 0, 0, 1 }, 
 { 998, 3, 1, 0, 0, 994, 0 }
};
const int weightEdge[28] = { 1, 3, 3, 1, 1, 1, 993, 995, 993, 2, 1000, 4, 995, 999, 1, 4, 1, 5, 997, 1, 1, 6, 1, 1, 998, 3, 1, 994 };
const int capacity[7][7] = { 
 { 0, 0, 11, 60, 54, 67, 73 }, 
 { 0, 0, 93, 75, 53, 87, 46 }, 
 { 11, 93, 0, 24, 37, 0, 73 }, 
 { 60, 75, 24, 0, 0, 0, 0 }, 
 { 54, 53, 37, 0, 0, 0, 0 }, 
 { 67, 87, 0, 0, 0, 0, 55 }, 
 { 73, 46, 73, 0, 0, 55, 0 }
};
const int demand[7][7] = { 
 { 0, 1, 0, 0, 0, 0, 0 }, 
 { 0, 0, 1, 0, 0, 0, 0 }, 
 { 1, 0, 0, 0, 2, 0, 0 }, 
 { 0, 0, 0, 0, 0, 0, 1 }, 
 { 0, 0, 0, 0, 0, 0, 0 }, 
 { 0, 0, 1, 0, 1, 0, 0 }, 
 { 1, 0, 0, 0, 0, 1, 0 }
};
const int edge[28][2] = { { 0, 2 }, { 0, 3 }, { 0, 4 }, { 0, 5 }, { 0, 6 }, { 1, 2 }, { 1, 3 }, { 1, 4 }, { 1, 5 }, { 1, 6 }, { 2, 0 }, { 2, 1 }, { 2, 3 }, { 2, 4 }, { 2, 6 }, { 3, 0 }, { 3, 1 }, { 3, 2 }, { 4, 0 }, { 4, 1 }, { 4, 2 }, { 5, 0 }, { 5, 1 }, { 5, 6 }, { 6, 0 }, { 6, 1 }, { 6, 2 }, { 6, 5 } };