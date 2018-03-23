/*! \file node_shuffle.c
 *  \brief This files contains the routines in order to create all sorts of 
 *         hints, sturctural information in order to generate sparse matricies
 *         later on.
 *  This file contains the functions in order to generate and manipulate the
 *  so called index matricies. These matricies are generate using the 
 *  overall nodal information. The index matricies give hints about how much 
 *  memory is needed for sparse matrix creation and how many entries these 
 *  sparse matricies will have.
 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"hex_fem_solver.h"

/*! \brief function that mappes each index in from a reduced result vector
 *         to a global result vector.
 *  This function generates a translation a translation table that mappes 
 *  reduced (i.e. a vector without containing values for fixed nodes) to a 
 *  global index (i.e. a vector having indicies for all nodes).
 *  \param object a sturcture structure containing the finite element problem
 */
int* generate_reduced_index_to_index_vector_translation(structure object) {

  int i,j;
  int* table = malloc(sizeof(int)*object.n_nodes);
  node* local_nodes = object.nodes;
  
  j=0;
  for (i=0;i<object.n_nodes;i++) {
    if ( local_nodes[i].fixed == 1 ) {
      table[i] = -1;
    } else {
      table[i] = j;
      j++;
    }
  }
  return(table);
}

/*! \brief function that creates the index and the reduced index matricies
           where reduced means without fixed nodes.
 *  This function generates both the reduced and non reduced index matricies.
 *  The reduced index matrix does not account for fixed nodes. The index
 *  matricies are used to generate the according sparse matricies later on, and
 *  to estimate their memory usage.
 *  \param object a sturcture structure containing the finite element problem
 *  \param reduced_translation_table a translation able obtained using
 *         generate_reduced_index_to_index_vector_translation()
 */
char** generate_index_and_reduced_index_matrix(structure object,
					       int* reduced_translation_table) {

  int i,j,k;

  int reduced_k, reduced_j;
  
  int* n_indicies;
  long long offset; // should speed things a bit up
  
  element ele_buffer;

  int r_size;
  node * local_nodes = object.nodes;
  int n_nodes = object.n_nodes;

  char* I_index;
  char* R_index;

  char** index_matricies = (char**)malloc(sizeof(char*)*2);
  
  I_index = (char*)malloc(sizeof(char)*(long long)object.n_nodes
			  *(long long)object.n_nodes);

  memset(I_index,0,sizeof(char)*(long long)object.n_nodes
	 *(long long)object.n_nodes);

  r_size = reduced_size(local_nodes,n_nodes);
  R_index = (char*)malloc((long long)r_size*(long long)r_size*sizeof(char));
  memset(R_index,0,(long long)r_size*(long long)r_size*sizeof(char));
  

  for(i=0;i<object.n_elements;i++) {
    
    n_indicies = object.elements[i].node_indicies;
    
    for(j=0;j<8;j++) {
      offset =
	(long long)object.n_nodes*(long long)n_indicies[j];
      for(k=0;k<8;k++) {
	I_index[(long long)offset+(long long)n_indicies[k]] =
	  1;
	if (local_nodes[n_indicies[j]].fixed != 1 &&
	    local_nodes[n_indicies[k]].fixed != 1 ) {
	  reduced_j = reduced_translation_table[n_indicies[j]];	  
	  reduced_k = reduced_translation_table[n_indicies[k]];
	  R_index[(long long)r_size*(long long)reduced_j
		  +(long long)reduced_k]=1;
	}
      }
    }
  }
  index_matricies[0] = I_index;
  index_matricies[1] = R_index;
  return(index_matricies);
}

/*! \brief function obtains the number of non zero elements in an index
 *         matrix
 *  This function obtains the number of non zero elements in an index matrix.
 *  The index matrix here is stored in one character per index form and not in 
 *  binary form.
 *  \param index_matrix The index matrix where to count non zero elements.
 *  \param range The range of the index matrix in question.
 */
int number_of_non_zero_elements_in_index_matrix(char* index_matrix,
						int range) {

  int counter = 0;
  int i;

  for(i=0;i<range*range;i++) {
    if (index_matrix[i] == 1) {
      counter++;
    }
  }
  return counter;
}

