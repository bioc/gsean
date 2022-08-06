/* NetPreProc.c
   C code for the R package NetPreProc
*/

#include <stdlib.h>

struct elem  {int value;
		      struct elem * next;
		     };
typedef struct elem ELEMENT;

static ELEMENT * tail = NULL;

/* Initialiation of a void list 
Input:
tail : pointer to the tail of the list. N.B. tail must be adummy allocated element of type ELEMENT
Output: pointer to the head of the list. The head points dorectly to the dummy tail
*/
ELEMENT * init_list(void) {
    if (tail == NULL) {
    tail = (ELEMENT *) malloc(sizeof(ELEMENT)); 
    tail->value=0;
    tail->next=tail;
  }
  ELEMENT * head = tail; 
  return(head);
}

/* Insertion of an element at the head of the list (the first place)
Input:
h : pointer to the head of the list. 
val: integer value to be inserted
Output:
Pointer to the head of the list (that now points to newly inserted element)
*/
ELEMENT * head_insert(ELEMENT * h, int val) {
  ELEMENT * x = (ELEMENT *) malloc(sizeof(ELEMENT)); 
  x->value=val;
  x->next = h;
  h = x;
  return(h);
}


/* Free the memory allocated forthe list
Input:
h : pointer to the head of the list. 
*/
void free_list(ELEMENT * h) {
  ELEMENT * p, *q;
  p = h;
  while(p!=tail) {
    q = p;
	p = p->next;
	free(q);  
  }    
}

/* It finds an element in the list
Input:
h : pointer to the head of the list. 
val : element to be found
Output
If val is prsent in the list it return 1, otherwise 0
*/
int find_list(ELEMENT * h, int val) {
  ELEMENT * p;
  for (p=h; p!=tail; p=p->next)
    if (p->value == val)
	  return(1);
  return(0);
}


/* It coints the elemnts common to 2 lists
Input:
h1 : pointer to the head of the first list. 
h2 : pointer to the head of the second list. 
Output
the number elements common to the to lists
*/
int n_common(ELEMENT * h1, ELEMENT * h2) {
  ELEMENT * p;
  int n = 0;
  for (p=h1; p!=tail; p=p->next) 
	 n += find_list(h2, p->value);
  return(n);
}


void  chua_like_norm(double * W,   int * n){
  register int i, j;
  int m = *n;
  int t1, t2, t3;
  ELEMENT * neighbour [m];  /* list of neighbours for each element */ 
 
  int n_neigh [m]; /* number of neighbours for each element */
  
  /* creation of list of neighbourhoods */
  for (i=0; i <m; i++) {
    neighbour[i] = init_list();
	n_neigh[i] = 0;
	W[i + i * m] = 0.0001;
	for (j=0; j < m; j++) 
	  if (W[i*m + j] > 0) {
	    neighbour[i] = head_insert(neighbour[i], j);
		n_neigh[i]++;
	  }
  }
  /* computing Chua normalization */
  for(i=0; i<(m-1); i++)
	for(j=i+1; j<m; j++) {
	  t1 = n_common(neighbour[i], neighbour[j]);
	  t2 = n_neigh[i] - t1;
	  t3 = n_neigh[j] - t1;	  
	  W[i*m + j] = ((double)(2*t1) / ((double)(t2 + 2*t1 +1)) ) * ( (double)(2*t1) / ((double)(t3 + 2*t1 +1)) );
	  W[j*m + i] = W[i*m + j];
	}
  for (i=0; i < m; i++) {
    W[i + i * m] = 0; 
	free_list(neighbour[i]);
  } 
}




/* 
  norm_lapl_graph: Normalized graph Laplacian.
  Given an adjacency matrix of a graph, it computes the corresponding Normalized graph Laplacian
  INPUT:
  W: pointer to a symmetric matrix
  diag: pointer to a vector representing the diagonal of the matrix D^-1/2, where
  D_ii = \sum_j W_ij and D_ij = 0 if i!=j
  n: dimension of the square matrix W
  OUTPUT
  W is the normalized matrix
*/
void  norm_lapl_graph(double * W, double * diag, int * n){
  register int i, j, k, m;
  double x;
  m = *n;
  for (i=0; i <m; i++) {
    x = diag[i];
    for (j=0; j <m; j++){
	   k = j*m + i;	   
	   W[k] = W[k] * x;	
	}
  }
  
  for (i=0; i <m; i++) {
    for (j=0; j <m; j++){ 
	   x = diag[j];
	   k = j*m + i;	   
	   W[k] = W[k] * x;
	}
  }
}
  

/* 
  norm_2: Probabilistic normalization.
  Given an adjacency matrix of a graph, it computes the corresponding Normalized graph:
  D^-1 * W, where D_ii = \sum_j W_ij and D_ij = 0 if i!=j
  INPUT:
  W: pointer to the adjacency matrix
  diag: pointer to a vector representing the diagonal of the matrix D^-1, where
  D_ii = \sum_j W_ij and D_ij = 0 if i!=j
  n: dimension of the square matrix W
  OUTPUT
  W is the normalized matrix
*/
void  norm_2(double * W, double * diag, int * n){
  register int i, j, k, m;
  double x;
  m = *n;
  for (i=0; i <m; i++) {
    x = diag[i];
    for (j=0; j <m; j++){
	   k = j*m + i;	   
	   W[k] = W[k] * x;	
	}
  }
}


