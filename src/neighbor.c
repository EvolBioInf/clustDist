/***** neighbor.c *********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Tue Sep 18 13:49:34 2012
 **************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "interface.h"
#include "eprintf.h"
#include "distTree.h"

void calcR(Matrix *mat);
void calcS(Matrix *mat);
void recalcDist(Args *args, Matrix *mat);
void findMin(Args *args, Matrix *mat);

void cluster(Args *args, Matrix *mat, Node **tree){
  int treeLen;
  float d;  
  Node *parent, *child1, *child2, *child3, *np;

  while(mat->n > 3){
    calcR(mat);
    calcS(mat);
    findMin(args,mat);
    if(args->m)
      printMatrix(args,mat);
    constrTree(args,mat,tree);
    recalcDist(args, mat);
    mat->n--;
  }
  if(args->u){
    findMin(args,mat);
    if(args->m)
     printMatrix(args,mat);
    constrTree(args,mat,tree);
    recalcDist(args,mat);
    mat->n--;
  }
  if(args->m)
    printMatrix(args,mat);
  /* finish tree */
  if(args->u)
    treeLen = 2 * mat->numTaxa - 1;
  else
    treeLen = 2 * mat->numTaxa - 2;
  parent = tree[treeLen-1];
  parent->isRoot = 1;
  child1 = tree[mat->ind[0]];
  parent->left = child1;
  child2 = tree[mat->ind[1]];
  np  = child1;
  while(np->right)
    np = np->right;
  np->right = child2;
  if(!args->u){
    np = np->right;
    child3 = tree[mat->ind[2]];
    np->right = child3;
    child1->dist = (mat->d[0][1] + mat->d[0][2] - mat->d[1][2])/2.;
    child2->dist = (mat->d[1][0] + mat->d[1][2] - mat->d[0][2])/2.;
    child3->dist = (mat->d[2][0] + mat->d[2][1] - mat->d[0][1])/2.;
  }else{
    d = mat->d[0][1]/2.;
    parent->height = d;
    child1->dist = d - child1->height;
    child2->dist = d - child2->height;
  }
}

void calcR(Matrix *mat){
  int i, j, n;
  
  n = mat->n;
  for(i=0;i<n;i++){
    mat->r[i] = 0;
    for(j=0;j<n;j++){
      mat->r[i] += mat->d[i][j];
    }
  }
}

void calcS(Matrix *mat){
  int i, j, n;

  n = mat->n;
  for(i=0;i<n-1;i++)
    for(j=i+1;j<n;j++)
      mat->s[i][j] = mat->d[i][j] - (mat->r[i] + mat->r[j])/(double)(n-2);
}

void findMin(Args *args, Matrix *mat){
  int i, j, n;
  double min;
  
  n = mat->n;
  min = DBL_MAX;

  for(i=0;i<n-1;i++)
    for(j=i+1;j<n;j++)
      if(args->u){
	if(min > mat->d[i][j]){
	  min = mat->d[i][j];
	  mat->minI = i;
	  mat->minJ = j;
	}
      }else{
	if(min > mat->s[i][j]){
	  min = mat->s[i][j];
	  mat->minI = i;
	  mat->minJ = j;
	}
      }
}
