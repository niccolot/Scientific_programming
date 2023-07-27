#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "barnes_static.h"

#define G 0.0000000000667 // gravitational coupling constant


int get_quadrant(double x, double y, double x0, double y0);
node_t *new_node(double x, double y, double m);
node_t *insert_aux(double m, double x, double y, node_t *root, double x0, double y0, int h);
double get_mass_aux(double x, double y, node_t *root, double x0, double y0, int h);
void get_force_aux(double x, double y, double m, double *fx, double *fy, double theta, node_t* root, int h);
double l2_norm(double x1, double y1, double x2, double y2);


int string_to_body(const char* s, double* x, double* y, double* m){
	/*
	extracts the value of the x,y coordinates and mass m for a body
	from a string formatted as "x y m" i.e. the values separated by 
	blank spaces e.g. "400 450 1200"
	
	returns 0 if the conversion happened succesfully
	returns -1 otherwise and do not change the variables x, y, m
	*/
	
	double x_temp, y_temp, m_temp;
	
	x_temp = *x;
	y_temp = *y;
	m_temp = *m;
	
	int retval = sscanf(s, "%lf %lf %lf", x, y, m);
	
	if(retval == 3) return 0;
	
		else{
			*x = x_temp;
			*y = y_temp;
			*m = m_temp;
		
			return -1;
		}
}


node_t* insert(double m, double x, double y, node_t *root){
	/*
	inserts a body of coordinates (x,y) and mass m in the tree
	with root pointed by 'root'
	
	returns the pointer to the updated tree
	*/
	
	
	int h = 1;
	double x0 = 0; 
	double y0 = 0;
	
	// empty root, first body to be inserted in the quadtree
	if(root == NULL){
		root = new_node(x,y,m);
		return root;
	}  

	return insert_aux(m,x,y,root,x0,y0,h);
}


void free_tree (node_t * root){
	/*
	deallocates the whole quad tree
	*/
	
	if(root != NULL){
		free_tree(root->NE);
		free_tree(root->SE);
		free_tree(root->SW);
		free_tree(root->NW);
		free(root);
	}
}


double get_mass(double x, double y, node_t* root){
	/*
	returns mass of the body in (x,y)
	
	else returns 0 if body of coordinates (x,y)
	is not in the tree
	*/
	
	if(root == NULL) return 0;
	
	if(root->x == x && root->y == y) return root->mass;
	
	double x0 = 0;
	double y0 = 0;
	int h = 1;
	
	
	return get_mass_aux(x,y,root,x0,y0,h);
	
}


void get_force(double x, double y, double *fx, double *fy, double theta, node_t* root){
	/*
	calculates the components (fx,fy) of the force on the particle of coordinates (x,y)
	given a tollerance 'theta'
	*/
	
	
	*fx = 0;
	*fy = 0;
	
	if(root == NULL) return;
	
	double m = get_mass(x,y,root);
	int h = 0;
	
	if(m == 0) return ;
	
	get_force_aux(x,y,m,fx,fy,theta,root,h);
	
	return;
}


////////////////////...UTILITY FUNCTIONS...////////////////////////////////
int get_quadrant(double x, double y, double x0, double y0){
	/*
	returns the right subquadrant where to insert the body
	based on the body coordinates (x,y) and the current 
	quadrant's coordinates (x0,y0)
	
	e.g. if the whole space has radius 1 and is subdivided
	into 4 quadrants, the NE quadrant will have center 
	coordinates (0.5,0.5)
	
	
	return 1 for NE quadrant
	return 2 for SE quadrant
	return 3 for SW quadrant
	return 4 for NW quadrant
	*/


	if(x > x0){
		if(y >= y0){
			return 1;
		}
		else return 2;
	}
	
	else{
		if(y < y0){
			return 3;
		}
		else return 4;
	}
}


node_t *new_node(double x, double y, double m){
	/*
	initializes a node
	*/
	
	node_t *temp;
	temp = malloc(sizeof(node_t));
	if(temp == NULL) return NULL;
	
	temp -> x = x;
	temp -> y = y;
	temp -> mass = m;
	temp -> NW = NULL;
	temp -> NE = NULL;
	temp -> SE = NULL;
	temp -> SW = NULL;
	
	return temp;
}


node_t *insert_aux(double m, double x, double y, node_t *root, double x0, double y0, int h){
	/*
	insert a body recursively keeping track of the current quadrant's gemetrical center
	(not mass center) (x0,y0) and of the depth in the tree, given by 'h' 
	*/
	
	
	// empty leaf, just insert the body
	if(root -> mass == 0){
		root = new_node(x,y,m);
		return root;
	}
	
	
	// leaf with a body in it
	//
	// from the leaf node another 4 quadrants are produced and 
	// the old body is placed in the right quadrant
	//
	// the insertion of the new body then proceeds recursively
	// with the case of central node
	if((root->NW) == NULL && (root->NE) == NULL && (root->SW) == NULL && (root->SE) == NULL){
		
		node_t *temp = new_node(root->x,root->y,root->mass);
		
		int pos = get_quadrant(temp->x,temp->y,x0,y0);
		
		if(pos == 1){
			
			// the old body is assigned to the right quadrant
			// and the other quadrant are inizialized with arbitrary 0 values
			root->NE = temp;
			root->SE = new_node(0,0,0);
			root->SW = new_node(0,0,0);
			root->NW = new_node(0,0,0);
		}
		
		else if(pos == 2){
			root->NE = new_node(0,0,0);
			root->SE = temp;
			root->SW = new_node(0,0,0);
			root->NW = new_node(0,0,0);
		}
		
		else if(pos == 3){
			root->NE = new_node(0,0,0);
			root->SE = new_node(0,0,0);
			root->SW = temp;
			root->NW = new_node(0,0,0);
		}
		
		else if(pos == 4){
			root->NE = new_node(0,0,0);
			root->SE = new_node(0,0,0);
			root->SW = new_node(0,0,0);
			root->NW = temp;
		}
	}
	
	
	// central node (it has at least a child node)
	// continue to search recursively for an empty leaf
	// deepening the search in the right quadrants depending
	// on the coordinates of the body to be inserted
	int pos = get_quadrant(x,y,x0,y0);
		
	if(pos == 1){
			
		//center of the new quadrant
		x0 += _s/pow(2,h);
		y0 += _s/pow(2,h);
			
		// insertion proceeds recursively in the new quadrant
		root->NE = insert_aux(m,x,y,root->NE,x0,y0,h+1);
	}
		
	else if(pos == 2){
		x0 += _s/pow(2,h);
		y0 -= _s/pow(2,h);
		root->SE = insert_aux(m,x,y,root->SE,x0,y0,h+1);
	}
		
	else if(pos == 3){
		x0 -= _s/pow(2,h);
		y0 -= _s/pow(2,h);
		root->SW = insert_aux(m,x,y,root->SW,x0,y0,h+1);
	}
		
	else if(pos == 4){
		x0 -= _s/pow(2,h);
		y0 += _s/pow(2,h);
		root->NW = insert_aux(m,x,y,root->NW,x0,y0,h+1);
	}
		
	// node's mass and mass center update
	root->x += x*m/(root->mass);
	root->x *= (root->mass)/(root->mass+m);
	root->y += y*m/(root->mass);
	root->y *= (root->mass)/(root->mass+m);
	root->mass += m;
	
	return root;
}


double get_mass_aux(double x, double y,node_t *root, double x0, double y0, int h){

	
	if(root == NULL) return 0;
	
	// check if the node is indeed a lead (i.e. a body and not a mass center)
	if(root -> NW == NULL && root -> NE == NULL && root -> SE == NULL && root -> SW == NULL){
	
		if(root->x == x && root->y == y) return root->mass;
		
		return 0;
	}
	
	// binary-like search in the appropriate quadrant
	int pos = get_quadrant(x,y,x0,y0);
	
	if(pos==1){
	
		// center of the next subquadrant
		//
		// universe_radius/pow(2,h) is the radius of a subquadrant 
		// of a node at depth 'h' in the quadtree
		x0 += _s/pow(2,h);  
		y0 += _s/pow(2,h);
		
		// search in the next subquadrant
		return get_mass_aux(x,y,root->NE,x0,y0,h+1);
	}
	
	else if(pos==2){
		x0 += _s/pow(2,h);
		y0 -= _s/pow(2,h);
		return get_mass_aux(x,y,root->SE,x0,y0,h+1);
	}
	
	else if(pos==3){
		x0 -= _s/pow(2,h);
		y0 -= _s/pow(2,h);
		return get_mass_aux(x,y,root->SW,x0,y0,h+1);
	}
	
	else if(pos==4){
		x0 -= _s/pow(2,h);
		y0 += _s/pow(2,h);
		return get_mass_aux(x,y,root->NW,x0,y0,h+1);
	}
	
	return 0;		
}


void get_force_aux(double x, double y, double m, double *fx, double *fy, double theta, node_t* root, int h){
	/*
	calculates recursively the force following the barnes-hut approximation depending on the tollerance 'theta'
	*/
	
	
	// invalid node
	if(root -> mass == 0 || (root -> x == x && root -> y == y)) return;
	
	double d = l2_norm(x,y,root->x,root->y);
	double size = _s/pow(2,h); // size of the current quadrant
	
	
	// the node's mass center is far enough from the  body
	// to be considered a single point
	// or we arrived to a leaf node
	if(size/d < theta || (root -> NW == NULL && root -> NE == NULL && root -> SE == NULL && root -> SW == NULL)){
		
		double Dx = root->x - x;
		double Dy = root->y - y;
		double f = G*m*(root->mass)/pow(d,2); // modulus of the force
		
		// force components
		*fx += f*Dx/d;
		*fy += f*Dy/d;
		
		return;
	}
	
	get_force_aux(x, y, m, fx, fy, theta, root->NE, h+1);
	get_force_aux(x, y, m, fx, fy, theta, root->SE, h+1);
	get_force_aux(x, y, m, fx, fy, theta, root->SW, h+1);
	get_force_aux(x, y, m, fx, fy, theta, root->NW, h+1);

}


double l2_norm(double x1, double y1, double x2, double y2){
	/*
	L2 distance between (x1,y1) and (x2,y2) vectors
	*/
	
	return sqrt(pow(x1-x2, 2) + pow(y1-y2, 2));
}
