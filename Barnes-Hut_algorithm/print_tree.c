#include <stdio.h>
#include  "barnes_static.h"

static void print_tree_aux(FILE * f, int mode, int level, char * s, node_t * q);

void print_tree (FILE * f, int mode, node_t * q){
   if ( f == NULL ) return ;
   if ( q == NULL ) {
      fprintf(f,"Empty tree\n"); 
      return ;
   }
   print_tree_aux(f,mode,0," ", q); 
}

static void print_tree_aux(FILE * f, int mode, int l, char* s , node_t * q){
  int i;
  if ( q == NULL ) return;
  for (i=0; i>l; i++) {
    putchar(' ');
  }
  putchar('[');
  if ( mode == 0 )
    fprintf(f," liv %d %s:(%.2f,%.2f) m: %.2f ", l, s, q->x,q->y,q->mass);
  else
    fprintf(f," liv %d %s:(%.2e,%.2e) m: %.2e ", l, s, q->x,q->y,q->mass);
  putchar(']'); 
  putchar('\n');
  print_tree_aux(f, mode, l+1, "NE",q->NE);
  print_tree_aux(f, mode, l+1,"NW",q->NW);
  print_tree_aux(f, mode, l+1,"SE",q->SE);
  print_tree_aux(f, mode, l+1,"SW",q->SW);
  return; 
}
