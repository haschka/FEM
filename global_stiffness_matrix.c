/*! \file global_stiffness_matrix.c
 * \brief The file that implements global stiffness matrix assembly.
 * This file contains the necessary functions that assembling the 
 * global stiffness matrix in full and reduced form ( without fixed nodes )
 * The global stiffness matrix is assembled by calculating stiffness matricies
 * for each element, summing them into the global stiffness matrix.
 */  

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<unistd.h>

#if defined(__SSSE3__) 
#include<immintrin.h>
#endif

#include"hex_fem_solver.h"

/*! \brief A function that allows to obtain the reduced size of ( 
 *         without fixed nodes ) of the finite element system.
 *  The reduced size function allows one to obtain the reduced size of the 
 *  finite element system, this is the number of nodes - the number of 
 *  fixed nodes.
 *  \param nodes A pointer to all the nodes of the system
 *  \param n_nodes The total number of all nodes in the system
 */
int reduced_size(node* nodes, int n_nodes) {
  int i,j;
  j=0;
  for(i=0;i<n_nodes;i++) {
    if ( nodes[i].fixed != 1 ) {
      j++;
    }
  }
  return j;
}

/*! \brief A function that allows one to create a foce vector from the force
 *         values stored in nodes that are not fixed.
 *  This function maps the forces applied onto non fixed nodes into a vector. 
 *  This is necessary in order to solve the finite element system using the 
 *  reduced stiffness matrix. 
 */
double * create_reduced_force_vector(structure object) {

  int i,j,k;
  int r_size = reduced_size(object.nodes, object.n_nodes);

  double * force = (double*)malloc(r_size*3*sizeof(double));

  node * nodes = object.nodes;
  j = 0;
  for (i=0;i<object.n_nodes;i++) {
    if (nodes[i].fixed != 1) {
      for(k=0;k<3;k++) {
	force[j*3+k]=nodes[i].f[k];
      }
      j++;
    }
  }
  return(force);
}

/*! \brief The function that assembles the local stiffness matrices into 
 *         global ones.
 *  This function generates both a global stiffness matrix and a reduced 
 *  global stiffness matrix in sparse form, by calling the necessary functions
 *  to calculate the local stiffness matricies and summing the results into 
 *  the global stiffness matricies. The function essentially allocates the 
 *  space necessary to store the sparse matricies by calling the functions 
 *  creating the index matrix, filling the sparse matricies then up with the
 *  results from the local stiffness matricies. 
 */   
two_sparse_matricies create_g_and_r_stiffness_matrix_sparse(structure object) {

  int i,j,k,l,m;
  
  element ele_buffer;
  
  double* K_local;

  two_sparse_matricies retval;
  
  sparse_matrix K_global; // sparse matrix type for sparse solvers
  sparse_matrix K_reduced;
  char** index_matricies; // container for both matricies
  char* I_global; // Index array for generating sparse matricies
  char* I_reduced;
  int* reduced_translation_table;

  int reduced_k, reduced_j;
  int r_size;
  
  reduced_translation_table =
    generate_reduced_index_to_index_vector_translation(object);
  
  index_matricies =
    generate_index_and_reduced_index_matrix(object, reduced_translation_table);
  
  I_global = index_matricies[0];
  I_reduced = index_matricies[1];
						         
  K_global = generate_sparse_matrix_from_index_matrix(I_global, object.n_nodes);
  free(I_global);
  order_sparse_matrix(K_global);

  r_size = reduced_size(object.nodes, object.n_nodes);
  
  K_reduced = generate_sparse_matrix_from_index_matrix(I_reduced,r_size);
  free(I_reduced);
  order_sparse_matrix(K_reduced);
  free(index_matricies);
  
  for(i=0;i<object.n_elements;i++) {
    
    ele_buffer = object.elements[i];
    K_local = create_stiffness_matrix_for_element(object,ele_buffer);
    
    for(j=0;j<8;j++) {
      for(k=0;k<8;k++) {
	for(l=0;l<3;l++) {
	  for(m=0;m<3;m++) {
	    add_to_value_in_sparse_matrix(K_global,
					  3*ele_buffer.node_indicies[k]+m,
					  3*ele_buffer.node_indicies[j]+l,
					  K_local[k*3*33+j*3+m*33+l]);
	  }
	}
	if ( object.nodes[ele_buffer.node_indicies[k]].fixed != 1 &&
	     object.nodes[ele_buffer.node_indicies[j]].fixed != 1 ) {
	  for(l=0;l<3;l++) {
	    for(m=0;m<3;m++) {

	      reduced_k =
		reduced_translation_table[ele_buffer.node_indicies[k]];
	      reduced_j =
		reduced_translation_table[ele_buffer.node_indicies[j]];

	      add_to_value_in_sparse_matrix(K_reduced,
					    3*reduced_k+m,
					    3*reduced_j+l,
					    K_local[k*3*33+j*3+m*33+l]);
	    }
	  }
	}
      }
    }
    free(K_local);
  }
  free(reduced_translation_table);
  retval.a = K_global;
  retval.b = K_reduced;

  //  K_global = drop_zeros_in_sparse_matrix(K_global);
  return(retval);
}

/*! \brief A function that calculates the forces applied onto the fixed nodes
 *         from a system already solved for displacements. 
 *  The forces_from_reduced function calculates the forces applied onto the
 *  nodes that are fixed (have no degrees of freedom) in our system
 *  \param object The structure structure of the finite element system
 *  \param gloabl_stiffness_matrix The global non reduced stiffness matrix
 */
void forces_from_reduced(structure object,
			 sparse_matrix* global_stiffness_matrix) {

  double* double_g_s_m;
  sparse_matrix sparse_g_s_m;
  
  sparse_g_s_m = *global_stiffness_matrix;
  sparse_g_s_m.eye = create_sparse_eye(sparse_g_s_m);
  sparse_g_s_m.has_eye = 1;
    
  int i,j,k,l;

  int n_nodes = object.n_nodes;
  
  node * nodes = object.nodes;

  for (j=0;j<n_nodes;j++) {
    if ( nodes[j].fixed == 1) {
      for (i=0;i<n_nodes;i++) {
	for(k=0;k<3;k++) {
	  for(l=0;l<3;l++) { 
	    nodes[j].f[k] +=
	      get_value_from_sorted_sparse_matrix_with_eye(sparse_g_s_m,
							   i*3+l,j*3+k)
	      *nodes[i].u[l];
	  }
	}
      }
    }
  }
}

/*! \brief a function that attributes the nodal displacements from reduced,
 *         without taking account fixed nodes to the form to the global form.
 *  The node_displacements_from_reduced function attributes the displacements
 *  to the nodes in the structure structure describing the finite element
 *  problem from the solution of the K u = F equation. u is effectively called
 *  reduced here, and holds the displacement of all non fixed nodes.
 *  \param object The structure structure containing the fem problem.
 *  \param reduced The displacements obtained from the solution K u = F
 *                 where u hold all the displacements for the non fixed 
 *                 vectors.
 */ 
void node_displacements_from_reduced(structure object, double* reduced) {
  int i, j, k;

  node * nodes = object.nodes;

  j = 0;
  
  for (i=0;i<object.n_nodes;i++) {
    if (nodes[i].fixed != 1) {
      for( k=0;k<3;k++) {
	nodes[i].u[k]=reduced[j*3+k];
      }
      j++;
    } else {
      for ( k=0;k<3;k++) {
      nodes[i].u[k] = 0.f;
      }
    }
  }
}

//
// Conjungated gradient solver
//

#if defined(__AVX2__) 
#include "cg_avx2_without_gather"
#elif defined(__SSSE3__)
#include "cg_sse"
#else 
#include "cg_scalar"
#endif
// #include "cg_avx2"

