/*! \file hex_fem_solver.c
 * \brief The main execution contorl file containing the main() function.
 * The hex_fem_solver.c file contains the main execution control functions,
 * as such main() is implemented in here and do_sparse_cg() which solves
 * the finite element system given by the input file. These functions are
 * very short in oder to allow a clean overview of the program.
 * Most called functions are then to be found in different files. 
 */
#include<stdio.h>
#include<stdlib.h>
#include"hex_fem_solver.h"

/*!
 * \brief This is a simple functions that prints out resulting forces
 *        and displacements to a file.
 * The function print_forces_displacements 
 * prints out the results of the finite element calculation to the file.
 * \param object The structure definition that contains the results. 
 * \param file	The file to print the results to.
 * \param what_to_print 0: print displacements, 1: print forces
 */
void print_forces_displacements(structure object, FILE* file,
				int what_to_print) {
  int i,j;
  node current_node;
  for( i = 0;i<object.n_nodes;i++) {
    current_node = object.nodes[i];
    switch(what_to_print) {
    case 0:
      fprintf(file,"Node #: %10d displ: x0: % 12.6E x1: % 12.6E x2: % 12.6E \n",
	      current_node.id,
	      current_node.u[0], current_node.u[1], current_node.u[2]);
      break;
    case 1:
      fprintf(file,"Node #: %10d force: x0: % 12.6E x1: % 12.6E x2: % 12.6E \n",
	      current_node.id,
	      current_node.f[0], current_node.f[1], current_node.f[2]);
      break;
    default:
      printf("Error in argument of - print_forces_displacements()\n");
    }
  }
  printf("\n");

}

/*! \brief This is the main control function that runs the program.
 *  It is called directly from main(). The function do_sparse_cg
 *  highlights the program flow calling subsequent functions, 
 *  like reading the mesh from the file, preprocessing it solving 
 *  the finite element problem. This function mainly guides the overall
 *  program flow and details are implemented in the called functions.
 *  \param inputfile The general inputfile to the program. 
 */ 
int do_sparse_cg(FILE* inputfile) {

  structure mesh = read_mesh_from_file(inputfile);

  two_sparse_matricies s_mats = create_g_and_r_stiffness_matrix_sparse(mesh);

  sparse_matrix global_K = s_mats.a;
  
  sparse_matrix reduced_K = s_mats.b;
  
  double* reduced_u = obtain_reduced_displacements_cg(mesh, reduced_K);

  printf("finished - reduction \n");
  
  node_displacements_from_reduced(mesh, reduced_u);

  printf("finished - displacements \n");
  
  forces_from_reduced(mesh, &global_K);

  printf("finished - forces \n");
  
  free(reduced_u);
  //  free(reduced_K);
  free_sparse_matrix(global_K);
  if(mesh.outputs.out_forces) {
    print_forces_displacements(mesh,stdout,0);
  }
  if(mesh.outputs.out_displacements){ 
    print_forces_displacements(mesh,stdout,1);
  }
  postprocessing(mesh);
  
  return(0);
}
  
int main(int argc, char** argv) {

  FILE* inputfile = fopen(argv[1],"r");

  int testvalue; 
  
  do_sparse_cg(inputfile);
  
  return(0);
}  
