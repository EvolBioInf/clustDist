/***** newick.c ***********************************
 * Description: 
 * Author: Bernhard Haubold, haubold@evolbio.mpg.de
 * Date: Tue Sep 18 08:07:43 2012
 **************************************************/
#include <stdio.h>
#include <float.h>
#include "distTree.h"

void printSubtree(Node *node);
void printNode(Node *node);

void printNewickTree(Node *node){
  printSubtree(node);
  printf(";\n");
}

void printSubtree(Node *node){
  Node *child;
  if(node){
    if(node->left){  /* internal node */
      printf("(");
      child = node->left;
      while(child->right){
	printSubtree(child);
	printf(",");
	child = child->right;
      }
      printSubtree(child);
      printf(")");
      if(!node->isRoot)
	printNode(node);
    }else            /* leaf */
      printNode(node);
  }
}

void printNode(Node *node){
  if(node->label)
    printf("%s",node->label);
  /* avoid negative branch lengths */
  if(node->dist >= 0.)
    printf(":%.6f",node->dist);
  else
    printf(":0.0");
}
