/***** clustDist.c **********************************************************
 * Description: Cluster distance matrix into phylogeny
 * Author: Bernhard Haubold
 * Email: haubold@evolbio.mpg.de
 * License: GNU General Public License, https://www.gnu.org/licenses/gpl.html
 * Date: Fri May 31 09:31:00 2019
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <assert.h>
#include "interface.h"
#include "eprintf.h"
#include "distTree.h"

void cluster(Args *args, Matrix *mat, Node **tree);
void printNewickTree(Node *node);

void scanFile(FILE *fp, Args *args){
  int i, j, n;
  Matrix *mat;
  Node **tree, *root;
  char *buf;

  buf = emalloc(1000);
  while(fscanf(fp,"%d",&n) != EOF){
    mat = newMatrix(n);
    tree = newTree(n);
    /* read distances */
    for(i=0;i<n;i++){
      assert(fscanf(fp,"%s",buf));
      tree[i]->label = strdup(buf);
      mat->label[i] = strdup(buf);
      for(j=0;j<n;j++)
        assert(fscanf(fp,"%lf",&(mat->d[i][j])));
    }
    /* average distances */
    for(i=1;i<n;i++) {
      for(j=0;j<i;j++){
        mat->d[i][j] = (mat->d[i][j] + mat->d[j][i])/2.;
        mat->d[j][i] = mat->d[i][j];
      }
    }
    cluster(args, mat, tree);
    if(args->u)
      root = tree[2*n-2];
    else
      root = tree[2*n-3];
    printNewickTree(root);
    freeMatrix(mat);
    freeTree(tree, n);
  }
  free(buf);
}

int main(int argc, char *argv[]){
  FILE *fp;
  Args *args = getArgs(argc, argv);

  setprogname2(argv[0]);
  if(args->v)
    printSplash(args);
  if(args->h || args->err)
    printUsage();
  if(args->nf == 0) {
    fp = stdin;
    scanFile(fp, args);
  } else {
    for(int i = 0; i < args->nf; i++) {
      fp = efopen(args->fi[i], "r");
      scanFile(fp, args);
      fclose(fp);
    }
  }
  freeArgs(args);
  free(progname());
  return 0;
}
