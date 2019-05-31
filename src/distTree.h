/***** distTree.h *********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Tue Sep 18 14:58:05 2012
 **************************************************/
#ifndef DISTTREE
#define DISTTREE

#include "interface.h"
/* define distance matrix */
typedef struct matrix{
  int n;
  double **d;
  double **d2;
  double *r;
  double **s;
  int *ind;
  char **label;
  int minI, minJ;
  int numTaxa;
} Matrix;

/* define tree node */
typedef struct node{
  struct node *left;
  struct node *right;
  char *label;
  double dist;
  double height;    /* tree height for upgma */
  int id;
  unsigned int  isRoot : 1;
} Node;

Node **newTree(int n);
void freeTree(Node **tree, int n);
void printMatrix(Args *args, Matrix *mat);
void freeMatrix(Matrix *mat);
Matrix *newMatrix(int n);
Node *newNode(int i, int root);
void constrTree(Args *args, Matrix *mat, Node **tree);
void recalcDist(Args *args, Matrix *mat);
#endif
