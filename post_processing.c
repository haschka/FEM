/*! \file post_processing.c
 *  \brief This file includes all the functions that can be evaluated
 *         after the finite element solution. i.e. obtain strains per element
 *         etc.
 */ 
#include<stdlib.h>
#include<stdio.h>
#include<unistd.h>
#include<string.h>
#include<math.h>
#include"hex_fem_solver.h"

/*! \brief function that fills a linear array of displacements for one 
 *         element
 *  The function elementry_displacements_from_nodal creates a single array
 *  for one element of the fem structure that holds all the displacments.
 *  It is basically a reshuffeling of displacements in order to allow easy 
 *  access to the displacements of an element's nodes without having to dive
 *  into the complete fem structure. Note that the array has to be allocated
 *  beforehand here. 
 */
void elementry_displacements_from_nodal(element ele, structure object,
					double* elementry_displacements) {

  int i,j;

  node current_node;
  
  for(i=0;i<8;i++) {
    current_node = object.nodes[ele.node_indicies[i]];
    for (j=0;j<3;j++) {
      elementry_displacements[i*3+j] = current_node.u[j];
    }
  }
}

/*! \brief Calculates the partial invers of the local stiffness matrix that
 *         leads to static condensation.
 *  This is basically the matrix method of the static condensation. In this 
 *  case the resulting matrix multiplied to the stiffness matrix would yeald
 *  static condensation. In practice this matrix is just needed for post 
 *  processing where it is applied to the strain-displacement relation in order
 *  to propagate the effect of static condensation into the strain-displacement
 *  relation when calculating strains from displacements during post processing.
 *  \param m In general the local stiffness matrix of an element in uncondensed
 *           form.
 */
double* condensation_inverse(double* m) {

  int i,j,k;
  double key_value;
  double row_divisor;

  int range = 33;
    
  double* in = (double*)malloc(sizeof(double)*1089);

  memset(in,0,1089*sizeof(double));
  
  for(i=0;i<33;i++) {
    in[i*33+i] = 1.;
  }
  
  for(i=32;i>23;i--) {
    key_value = m[i*range+i];
    for (j=0;j<range;j++) {
      m[i*range+j] /= key_value;
      in[i*range+j] /= key_value;
    }
    for (j=0;j<i;j++) {
      row_divisor = m[j*range+i];
      for (k=0;k<range;k++) {
	m[j*range+k] -= row_divisor*m[i*range+k];
	in[j*range+k] -= row_divisor*in[i*range+k];
      }
    }
    for (j=i+1;j<range;j++) {
      row_divisor = m[j*range+i];
      for (k=0;k<range;k++) {
	m[j*range+k] -= row_divisor*m[i*range+k];
	in[j*range+k] -= row_divisor*in[i*range+k];
      }
    }   
  }
  return(in);
}

/*! \brief Generates a strain-displacement relation that includes the effect
 *         of pseudo strain displacent by condensation.
 *  The B_condensate function applies the partial invers of the local stiffness
 *  matrix onto the strain-displacement and pseudo-strain-displacement (G)
 *  matricies that the pseudo strain-displacements are condensated into
 *  the strain displacement matrix. This is necessary in order to obtain
 *  correct strain-displacements from displacenets evalated using a stiffness
 *  matrix modified by static condensation.
 *  \param strain_displacement The uncondensed strain displacement matrix.
 *  \param pseudo_strain_displacement The pseudo strain displacement matrix
 *                                    caused by incompatible displacement modes.
 *  \param cond_inv The partial invers of the local stiffness matrix for the
 *                  element in question.
 */
double* B_condensate(double* strain_displacement,
		     double* pseudo_strain_displacement,
		     double* cond_inv) {


  int i,j,k;
  
  double *B_condensate = (double*)malloc(sizeof(double)*144);
  memset(B_condensate,0,sizeof(double)*144);
  
  double *B = strain_displacement;
  double *G = pseudo_strain_displacement;

  for(i=0;i<6;i++){
    for(j=0;j<24;j++){
      for(k=0;k<24;k++) {
	B_condensate[i*24+j]+=B[i*24+k]*cond_inv[j*33+k];
      }
      for(k=24;k<33;k++) {
	B_condensate[i*24+j]+=G[i*9+(k-24)]*cond_inv[j*33+k];
      }
    }
  }
  return(B_condensate);
}

/* \brief calculates the strain tensor for every corner of an element.
 * This function creates an array containing the strain tensor for every corner
 * of an element ( where the nodes are, but note that shared nodes might 
 * have a different strain field for each element ). The returned matrix has 8
 * colums and 6 rows. Each column contains the strain tensor for each 
 * of the corresponding 8 nodes. The six rows are the elements of the strain
 * tensor as follows: epsilon_11, epsilon_22, epsilon_33, epsilon_23,
 * epsilon_13, epsilon_12.
 * \param ele The element to calculate the strain tensors for.
 * \param elements_and_nodes the complete solved fem structure structure.
 * \param elementry_displacements The displacements casted into an array for 
 *                                the element in questions. This array can 
 *                                be obtained by calling 
 *                                elementry_displacements_from_nodal().
 */
double* generate_new_strain_field_for_element(element ele,
					      structure elements_and_nodes,
					      double* elementry_displacements) {

  // Stain Field - Matrix:

  //        Node 1  Node 2  Node 3  Node 4  Node 5  Node 6  Node 7  Node 8
  //  e11    sf[0]   sf[1]   sf[2]   sf[3]   sf[4]   sf[5]   sf[6]   sf[7]
  //  e22    sf[8]   sf[9]  sf[10]  sf[11]  sf[12]  sf[13]  sf[14]  sf[15]
  //  e33   sf[16]  sf[17]  sf[18]  sf[19]  sf[20]  sf[21]  sf[22]  sf[23]
  // 2e23   sf[24]  sf[25]  sf[26]  sf[27]  sf[28]  sf[29]  sf[30]  sf[31]
  // 2e13   sf[32]  sf[33]  sf[34]  sf[35]  sf[36]  sf[37]  sf[38]  sf[39]
  // 2e12   sf[40]  sf[41]  sf[42]  sf[43]  sf[44]  sf[45]  sf[46]  sf[47]
  
  int i,j,n,m,l,i_three;
  int index;
  
  double* jacobi;
  double* jac_null;
  double* jac_null_inv;
  double* stress_strain;
  double* strain_displacement[8];
  double* pseudo_strain_displacement[8];
  double* q_points;

  double* condensed_strain_displacement;
  
  // strain_displacement_transpose . stress_strain . strain_displacement
  double* BT_D_B;
  double* BT_D_G;
  double* GT_D_G;
  double* D_G;
  
  double* cond_inv;

  double xi, eta, zeta, jacobi_det, jac_null_det, jac_null_jac_quot;


  int * local_nodes_index;
  node * local_nodes = (node*)malloc(sizeof(node)*8);

  double* k = malloc(sizeof(double)*1089);

  double* strain_field = (double*)malloc(sizeof(double)*48);

  memset(strain_field,0, sizeof(double)*48);
  
  memset(k,0,1089*sizeof(double));
  
  local_nodes_index = ele.node_indicies;

  for(i=0;i<8;i++) {
    local_nodes[i] = elements_and_nodes.nodes[ele.node_indicies[i]];
  }

  jac_null = create_jacobi(local_nodes,0.,0.,0.);

  jac_null_det = det_three_by_three(jac_null);
  
  jac_null_inv =
    create_invers_three_by_three(jac_null);

  free(jac_null);
  
  q_points = generate_quadrature_points();
  stress_strain = create_stress_strain(ele);
  
  for(i=0;i<8;i++) {

    i_three = i*3;
    
    xi = q_points[i_three];
    eta = q_points[i_three+1];
    zeta = q_points[i_three+2];
    
    jacobi = create_jacobi(local_nodes,xi,eta,zeta);
    strain_displacement[i] = create_strain_displacement(local_nodes,jacobi,
						     xi,eta,zeta);

    jacobi_det = det_three_by_three(jacobi);
    jac_null_jac_quot = jac_null_det/jacobi_det;
    
    pseudo_strain_displacement[i] =
      create_pseudo_strain_displacement(local_nodes,jac_null_inv,
					jac_null_jac_quot,
					xi,eta,zeta);

    BT_D_B = create_BT_D_B(strain_displacement[i],stress_strain);

    D_G = create_D_G(stress_strain,pseudo_strain_displacement[i]);

    GT_D_G=create_GT_D_G(pseudo_strain_displacement[i], D_G);

    BT_D_G=create_BT_D_G(strain_displacement[i],D_G);

    free(D_G);

    // Matrix Block
    //
    //         BT_D_B    BT_D_G 
    //  K += 
    //        (BT_D_G)T  GT_D_G
    
    for(n=0;n<24;n++) {
      for(m=0;m<24;m++) {
	k[n*33+m] += BT_D_B[n*24+m]*jacobi_det;
      }
      for(m=24;m<33;m++) {
	k[n*33+m] += BT_D_G[n*9+(m-24)]*jacobi_det;
      }
    }
    for(n=24;n<33;n++) {
      for(m=0;m<24;m++) {
	k[n*33+m] += BT_D_G[m*9+(n-24)]*jacobi_det;
      }
      for(m=24;m<33;m++) {
	k[n*33+m] += GT_D_G[(n-24)*9+(m-24)]*jacobi_det;
      }
    }
    free(jacobi);
    free(BT_D_B);
    free(BT_D_G);
    free(GT_D_G);
  }

  cond_inv = condensation_inverse(k);
  
  for(i=0;i<8;i++) {
    condensed_strain_displacement = B_condensate(strain_displacement[i],
						 pseudo_strain_displacement[i],
						 cond_inv);
    
    for(j=0;j<6;j++) {
      for(l=0;l<24;l++) {
	// sum over all displacements (24) in order to get
	// epsilons in one point
	strain_field[j*8+i] += condensed_strain_displacement[j*24+l]
	  *elementry_displacements[l];
      }
    }
    
    free(strain_displacement[i]);
    free(pseudo_strain_displacement[i]);
  }
  free(cond_inv);
  free(jac_null_inv);
  free(q_points);
  free(stress_strain);
  free(condensed_strain_displacement);
  return(strain_field);
}					    


/*! \brief Function that calculates the volume of a hexaheadral element after
 *         deformation.
 *  This function calculates the volume of a hexheadral element taking 
 *  the displacements from the finite element analysis into account.
 *  Even though in order to obtain the displacements bubble functions
 *  have been used, the deformation (wobbeling cause by the bubble functions) 
 *  is not taken into account any further.
 *  The volume is calculated by using the gaussian quadrature algorithm for 
 *  a hexaheadral element.
 *  \param ele The element to obtain the volume from
 *  \param object The whole solved FEM structure structure.
 */
double volume_of_hexahedral_element_sortie(element ele, structure object) {

  int i,j,i_three;
  double xi, eta, zeta, jacobi_det;
  node * local_nodes = (node*)malloc(sizeof(node)*8);
  node this_node;
  
  double volume = 0;

  double * jacobi;
  double * q_points = generate_quadrature_points();
  
  for(i=0;i<8;i++) {
    local_nodes[i] = object.nodes[ele.node_indicies[i]];
    for(j=0;j<3;j++) {
      local_nodes[i].x[j] += local_nodes[i].u[j];
    }
  }
  
  for(i=0;i<8;i++) {
    
    i_three = i*3;
    xi = q_points[i_three];
    eta = q_points[i_three+1];
    zeta = q_points[i_three+2];

    jacobi = create_jacobi(local_nodes,xi,eta,zeta);
    jacobi_det = det_three_by_three(jacobi);

    volume += jacobi_det;
    free(jacobi);
  }
  
  free(local_nodes);
  free(q_points);
  
  return(volume);
}

/*! \brief Function that calculates the volume of a hexaheadral element before
 *         deformation.
 *  This function calculates the volume of a hexheadral element in its 
 *  original form as specified in the input file.
 *  \param ele The element to obtain the volume from
 *  \param object The whole solved FEM structure structure.
 */
double volume_of_hexahedral_element_entree(element ele, structure object) {

  int i,i_three;
  double xi, eta, zeta, jacobi_det;
  node * local_nodes = (node*)malloc(sizeof(node)*8);

  double volume = 0;

  double * jacobi;
  double * q_points = generate_quadrature_points();
  
  for(i=0;i<8;i++) {
    local_nodes[i] = object.nodes[ele.node_indicies[i]];
  }

  for(i=0;i<8;i++) {
    
    i_three = i*3;
    xi = q_points[i_three];
    eta = q_points[i_three+1];
    zeta = q_points[i_three+2];

    jacobi = create_jacobi(local_nodes,xi,eta,zeta);
    jacobi_det = det_three_by_three(jacobi);

    volume += jacobi_det;
    free(jacobi);
  }
  
  free(local_nodes);
  free(q_points);
  
  return(volume);
}

/*! \brief Function that casts the strain field of an element into a 
 *         stain tensor for an element
 *  This function sums accross the columns to calculate an approximative
 *  strain tensor for an element. 
 *  \param strainfield the strain field as obtained by calling 
 *                     generate_new_strain_field_for_element()
 */      
double* generate_elementry_strain_tensor_from_field(double* strainfield) {
  int i;
  double *elementry_strain_matrix = (double*)malloc(9*sizeof(double));
  memset(elementry_strain_matrix,0,sizeof(double)*9);
  
  for(i=0;i<8;i++) {  

    elementry_strain_matrix[0] += strainfield[i];
    elementry_strain_matrix[4] += strainfield[8+i];
    elementry_strain_matrix[8] += strainfield[16+i];
    elementry_strain_matrix[1] += strainfield[24+i];
    elementry_strain_matrix[5] += strainfield[32+i];
    elementry_strain_matrix[2] += strainfield[40+i];

  }

  elementry_strain_matrix[0] /= 8;
  elementry_strain_matrix[4] /= 8;
  elementry_strain_matrix[8] /= 8;

  elementry_strain_matrix[5] /= 16; // divided by 2*8 
  elementry_strain_matrix[2] /= 16;
  elementry_strain_matrix[1] /= 16;

  elementry_strain_matrix[3] = elementry_strain_matrix[1];
  elementry_strain_matrix[6] = elementry_strain_matrix[2];
  elementry_strain_matrix[7] = elementry_strain_matrix[5];

  return(elementry_strain_matrix);
}

/*! \brief a function that calculates the eigenvalues of a strain tensor
 *  This function calculates the eigenvalues of a strain tensor, hence, 
 *  it yealds the principal strains for an element.
 */
double* create_principal_strain_from_tensor(double* elementry_strain_matrix){

  return(create_eigenvalues_three_by_three_real_sym(elementry_strain_matrix));
}

/*! \brief The main postprocessing entry function.
 *  This function calls the postprocessing functions as defined 
 *  by required outputs in the input file
 *  \param object the solved fem structure structure
 */
void postprocessing(structure object) {

  int i, retval;

  double* local_strain_field;
  double* strain_tensor;
  double* eigenvalues;
  double elementry_displacements[24];

  element current_element;
  
  printf("\n\nPostprocessing: \n");

  for(i=0;i<object.n_elements;i++) {

    current_element = object.elements[i];
    
    printf("Element %i:\n",i);
    if (object.outputs.out_volumes_before) {
      printf("   Volume before deformation: % 12.6E\n",
	     volume_of_hexahedral_element_entree(current_element, object));
    }
    if (object.outputs.out_volumes_after) {
      printf("    Volume after deformation: % 12.6E\n",
	     volume_of_hexahedral_element_sortie(current_element, object));
    }
    
    if (object.outputs.out_strain_tensor ||
	object.outputs.out_principal_strains) {

      elementry_displacements_from_nodal(object.elements[i], object,
					 elementry_displacements);
      local_strain_field =
	generate_new_strain_field_for_element(object.elements[i],
					      object, elementry_displacements);

      strain_tensor =
	generate_elementry_strain_tensor_from_field(local_strain_field);

      free(local_strain_field);

      if (object.outputs.out_strain_tensor) {
	printf("   strain tensor:  % 12.6E % 12.6E % 12.6E \n"
	       "                   % 12.6E % 12.6E % 12.6E \n"
	       "                   % 12.6E % 12.6E % 12.6E \n",
	       strain_tensor[0], strain_tensor[1], strain_tensor[2],
	       strain_tensor[3], strain_tensor[4], strain_tensor[5],
	       strain_tensor[6], strain_tensor[7], strain_tensor[8]);
      }
      
      if (object.outputs.out_principal_strains) {
	eigenvalues = create_principal_strain_from_tensor(strain_tensor);
	printf("   principal strains:  % 12.6E % 12.6E % 12.6E \n",
	       eigenvalues[0], eigenvalues[1], eigenvalues[2]);
      }
      
      free(strain_tensor);
      free(eigenvalues);
    }
  }
}
