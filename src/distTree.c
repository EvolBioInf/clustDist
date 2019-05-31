/***** distTree.c *********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Tue Sep 18 14:58:02 2012
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "interface.h"
#include "distTree.h"
#include "eprintf.h"

Node **newTree(int n){
  Node **tree;
  int i, numNodes;

  numNodes = 2 * n - 1; /* allocate enough nodes for rooted or unrooted tree */
  tree = (Node **)emalloc(numNodes*sizeof(Node *));
  for(i=0;i<numNodes-1;i++)
    tree[i] = newNode(i,0);
  tree[i] = newNode(i,1);  /* add root */
  return tree;
}

void freeTree(Node **tree, int n){
  int i, numNodes;

  numNodes = 2*n-1;
  for(i=0;i<numNodes;i++){
    if(tree[i]->label != NULL)
      free(tree[i]->label);
    free(tree[i]);
  }
  free(tree);
}

void printMatrix(Args *args, Matrix *mat){
  int i, j, n, max, len;
  char *label;
  static int round = 0;

  printf("****** Round %d ******\n",++round);
  n = mat->n;
  max = 0;
  for(i=0;i<n;i++){
    len = strlen(mat->label[mat->ind[i]]);
    if(len > max)
      max = len;
  }
  for(i=0;i<n;i++){
    label = mat->label[mat->ind[i]];
    printf("%s",label);
    for(j=strlen(label);j<max+2;j++)
      printf(" ");
    for(j=0;j<n;j++){
      if(i<j)
  	printf("%6.2f ",mat->d[i][j]);
      else if(i>j){
  	if(args->u || n <= 3)
  	  printf("%6.2f ",mat->d[i][j]);
  	else if(n > 3)
  	  printf("%6.2f ",mat->s[i][j]);
      }else
  	printf("%6d ",0);
    }
    if(!args->u && n > 3)
      printf("%5.2f",mat->r[i]);
    printf("\n");
  }
}

void freeMatrix(Matrix *mat){
  int i, n, numNodes;

  n = mat->numTaxa;
  numNodes = 2*n-1;
  for(i=0;i<n;i++){
    free(mat->d[i]);
    free(mat->d2[i]);
    free(mat->s[i]);
  }
  for(i=0;i<numNodes;i++)
    if(mat->label[i] != NULL)
      free(mat->label[i]);

  free(mat->d);
  free(mat->d2);
  free(mat->s);
  free(mat->label);
  free(mat->r);
  free(mat->ind);
  free(mat);
}

Matrix *newMatrix(int n){
  Matrix *mat;
  int i, j, numNodes;

  mat = (Matrix *)emalloc(sizeof(Matrix));
  mat->n = n;
  numNodes = 2*n-1;
  mat->numTaxa = n;
  mat->minI = -1;
  mat->minJ = -1;
  mat->d = (double **)emalloc(n*sizeof(double *));
  mat->d2 = (double **)emalloc(n*sizeof(double *));
  mat->s = (double **)emalloc(n*sizeof(double *));
  mat->label = (char **)emalloc(numNodes*sizeof(char *));
  for(i=0;i<n;i++){
    mat->d[i] = (double *)emalloc(n*sizeof(double));
    mat->d2[i] = (double *)emalloc(n*sizeof(double));
    mat->s[i] = (double *)emalloc(n*sizeof(double));
    for(j=0;j<n;j++)
      mat->d[i][j] = mat->d2[i][j] = mat->s[i][j] = 0.;
  }
  for(i=0;i<numNodes;i++)
    mat->label[i] = NULL;
  mat->r = (double *)emalloc(n*sizeof(double));
  mat->ind = (int *)emalloc(n*sizeof(int));
  for(i=0;i<n;i++)
    mat->ind[i] = i;

  return mat;
}

Node *newNode(int i, int root){
  Node *n;

  n = (Node *)emalloc(sizeof(Node));
  n->left = NULL;
  n->right = NULL;
  n->dist = 0.;
  n->label = NULL;
  n->id = i;
  if(root)
    n->isRoot = 1;
  else
    n->isRoot = 0;
  n->height = 0.;

  return n;
}

void constrTree(Args *args, Matrix *mat, Node **tree){
  int i, j, mi, mj, t, n;
  Node *parent, *child1, *child2, *np;
  double d;

  n = mat->n;
  mi = mat->minI;
  mj = mat->minJ;
  t = 2*mat->numTaxa - n;
  mat->label[t] = (char *)emalloc(512*sizeof(char));
  mat->label[t][0] = '\0';
  strncat(mat->label[t],mat->label[mat->ind[mi]],255);
  strcat(mat->label[t],",");
  strncat(mat->label[t],mat->label[mat->ind[mj]],255);
  parent = tree[t];
  child1 = tree[mat->ind[mi]];
  child2 = tree[mat->ind[mj]];
  parent->left = child1;
  np = child1;
  while(np->right)
    np = np->right;
  np->right = child2;
  for(np=child1;np->right != NULL; np=np->right)
    np->right = child2;
  if(args->u){
    d = mat->d[mi][mj]/2.;
    parent->height = d;
    child1->dist = d - child1->height;
    child2->dist = d - child2->height;
  }else{
    child1->dist = 1./2./(n-2)*((n-2)*mat->d[mi][mj] + mat->r[mi] - mat->r[mj]);
    child2->dist = 1./2./(n-2)*((n-2)*mat->d[mi][mj] + mat->r[mj] - mat->r[mi]);
  }
  /* adjust indexes */
  j = 0;
  for(i=0;i<n;i++){
    if((i != mi) && (i != mj))
      mat->ind[j++] = mat->ind[i];
  }
  mat->ind[j] = t;
}

void recalcDist(Args *args, Matrix *mat){
  int i, ii, j, jj, n;
  double mm;

  n = mat->n;
  for(i=0,ii=0;i<n;i++){
    if(i != mat->minI && i != mat->minJ){
      for(j=0,jj=0;j<n;j++){
	if(j != mat->minJ && j != mat->minI){
	  mat->d2[ii][jj] = mat->d[i][j];
	  jj++;
	}
      }
      ii++;
    }
  }
  mm = mat->d[mat->minI][mat->minJ];
  for(i=0,ii=0;i<n;i++){
    if(i != mat->minI && i != mat->minJ){
      if(args->u){
	mat->d2[n-2][ii] = (mat->d[mat->minI][i] + mat->d[mat->minJ][i])/2.;
	mat->d2[ii][n-2] = mat->d2[n-2][ii];
      }else{
	mat->d2[n-2][ii] = (mat->d[mat->minI][i] + mat->d[mat->minJ][i] - mm)/2.;
	mat->d2[ii][n-2] = mat->d2[n-2][ii];
      }
      ii++;
    }
  }
  for(i=0;i<n-1;i++){
    mat->d[i][i] = 0;
    for(j=0;j<n-1;j++)
      /* mat->d[i][j] = mat->d2[i][j]; */
      mat->d[i][j] = mat->d2[i][j] > 0. ? mat->d2[i][j] : 0.;
  }
}
