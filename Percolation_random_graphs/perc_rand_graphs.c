#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


struct ConnComp{
    struct ConnComp *parent; // node identifying the cluster
    unsigned int size; // size of the connected cluster
};


struct ConnComp *initialize(unsigned long long int n);
struct ConnComp *component_of(struct ConnComp *pcomp);
unsigned int merge_components(struct ConnComp *pcomp1, struct ConnComp *pcomp2);
void generate_list(double *list, int m);
unsigned long long int get_node(unsigned long long int n);


int main(){

    unsigned long long int n = 1000; // nodes in the graph
    // array of connected components, every element is a node but not
    // all nodes corresponds to the parent of a connected component,
    // only when comp->parent = comp is the case
    struct ConnComp *comp;
    unsigned long long int mean_clust_size;
    unsigned int max_clust_size;
    int m = 100; // number of order parameter values used in the simulation
    double c_list[m];
    int measures = 1000; // number of datapoints for each order parameter value
    unsigned long long int site1, site2;
    struct ConnComp *pcomp1, *pcomp2;
    unsigned int new_size;
    srand(time(0));
    FILE *pf_trajectories;

    generate_list(c_list, m);

    pf_trajectories = fopen("n1000.txt", "w");
    
    // measures
    for(int i=0; i<measures; i++){
        
		// trajectory in function of the order parameter
		for(int c=0; c<m; c++){

            comp = initialize(n);
			mean_clust_size = (unsigned long long int) n;
			max_clust_size = (unsigned int) 1;
            
			// evolution of the graph, a graph of average degree c and N nodes has cN/2 links 
			for(int j=0; j < (int) (c_list[c]*n*0.5); j++){
                
				site1 = get_node(n);
                
                do{
                    site2 = get_node(n);
                }while(site1 == site2);
                
				pcomp1 = component_of(comp + site1);
                pcomp2 = component_of(comp + site2);
                
				if(pcomp1 != pcomp2){

					mean_clust_size += 2*(unsigned long long int) pcomp1->size*pcomp2->size;
					new_size = merge_components(pcomp1, pcomp2);
	
					if(new_size > max_clust_size){
						max_clust_size = new_size;
					}
				}            
			}

            // the (square) of the largest cluster size have to be subtracted from mean_clust_size
            // in order to remove the dominating component and actually see the divergence for c=1
            fprintf(pf_trajectories, "%f\t%f\t%f\n", c_list[c], (float)  (mean_clust_size-pow(max_clust_size,2))/n, (float) max_clust_size/n);
            free(comp);
		}

        fprintf(pf_trajectories, "\n");
    }


    fclose(pf_trajectories);

    return 0;
}


struct ConnComp *initialize(unsigned long long int n){
    /*
    at the beginning every separate node is a single
    connected component of size 1
    */

    struct ConnComp *comp;
    comp = malloc(sizeof(struct ConnComp)*n);

    for(int i=0; i<n; i++){
        comp[i].parent = comp + i;
        comp[i].size = (unsigned int) 1;
    }

    return comp;
}


struct ConnComp *component_of(struct ConnComp *pcomp){
    /*
    returns the connected component a certain element belongs to
    */

    while(pcomp->parent != pcomp){
        pcomp = pcomp->parent;
    }
    return pcomp;
}


unsigned int merge_components(struct ConnComp *pcomp1, struct ConnComp *pcomp2){
    /*
    merges the smallest cluster into the largest
    */

    if(pcomp1->size < pcomp2->size){
        pcomp1->parent = pcomp2;
        pcomp2->size += pcomp1->size;
        return pcomp2->size;
    }
    else{
        pcomp2->parent = pcomp1;
        pcomp1->size += pcomp2->size;
        return pcomp1->size;
    }
}


void generate_list(double *list, int m){
    /*
    list of order parameter c (average degree of the graph) values
    */

    for(int i=0; i<m; i++){
        list[i] = 0.02*i;
    }
}


unsigned long long int get_node(unsigned long long int n){
    /*
    random node index for n nodes
    */
    
    unsigned long long int x;
    x = rand();
    x <<= 15;
    x ^= rand();
    x %= n;

    return x;
}
