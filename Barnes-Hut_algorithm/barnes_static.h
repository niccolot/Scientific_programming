#ifndef __BARNES_STATIC__H
#define __BARNES_STATIC__H
#include <stdio.h>


// global variable defining the radius of the universe
double _s;

// quadtree structure
typedef struct node {

  // mass center coordinates (if central node)
  // or body coordinates if leaf node
  double x, y; 

  // total mass if central node or body mass if leaf node
  double mass; 
  
  // pointer to quadrants
  struct node * NW, *NE, *SE, *SW;
} node_t;


	/*
	extracts the value of the x,y coordinates and mass m for a body
	from a string formatted as "x y m" i.e. the values separated by 
	blank spaces e.g. "400 450 1200"
	
	returns 0 if the conversion happened succesfully
	returns -1 otherwise and do not change the variables x, y, m
	*/
int body_to_string (const char* s, double* x, double* y, double* m);


	/*
	inserts a body of coordinates (x,y) and mass m in the tree
	with root pointed by 'root'
	
	returns the pointer to the updated tree
	*/
node_t* insert (double m, double x, double y, node_t * root);


  /*
  prints the quadtree pointed to root in a file (e.g. stdout)
  */
void print_tree (FILE * f, int mode, node_t * root);


	/*
	deallocates the whole quad tree
	*/
void free_tree (node_t * root);


	/*
	returns mass of the body in (x,y)
	
	else returns 0 if body of coordinates (x,y)
	is not in the tree
	*/
double get_mass (double x, double y, node_t* root);


	/*
	calculates the components (fx,fy) of the force on the particle of coordinates (x,y)
	given a tollerance 'theta'
	*/
void get_force (double x, double y, double *fx, double* fy, double theta, node_t* root);
#endif

