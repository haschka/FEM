/*! \file stiffness_matrix_for_element.c
 * \brief This file implements the calculation routines in order to obtain
 *        the stiffness matrix for a single element. 
 * This file contains the functions that are necessary to calculate the 
 * stiffness matrix for a single hexaheadral element. The stiffness matricies
 * of single hexaheadral elements are then summed later on into the global 
 * stiffness matrix, due to linearity.
 */

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include"hex_fem_solver.h"

/*! \brief function that generates gaussian quadrature points for a 
 *         hexaheadral element.
 */
double* generate_quadrature_points() {

  int i;

  int prefix_table[24] = { -1, -1, -1,
			    1, -1, -1,
			    1,  1, -1,
			   -1,  1, -1,
			   -1, -1,  1,
			    1, -1,  1,
			    1,  1,  1,
			   -1,  1,  1 };

  double* table = malloc(sizeof(double)*24);  

  double factor = 1/sqrt(3.f);
  
  for(i=0;i<24;i++){
    table[i] = factor*prefix_table[i];
  }
  return(table);
} 

/*! \brief creates the jacobi matrix for transformations 
 *         partial(x,y,z)/partial(xi,eta,zeta) for a hexaheadral element.
 *  This function calculates the jacobi matrix
 *  partial(x,y,z)/partial(xi,eta,zeta) at point xi, eta, zeta for a given 
 *  hexaheadral element.
 *  \param local_nodes a pointer to an array containing the eight nodes 
 *                     corresponding to this element.
 *  \param xi xi coordinate value.
 *  \param eta eta coordinate value.
 *  \param zeta zeta coordinate value.
 */
double* create_jacobi(node * local_nodes, // elements nodes (0-8)
		      double xi,
		      double eta,
		      double zeta) {

  int i,j; // matrix indicies (j*3+i) J(i,j)
  int k; // derived shape function indicies
  double* jacobi = (double*) malloc(sizeof(double)*9);

  memset(jacobi,0,9*sizeof(double));
  
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      for (k=0;k<8;k++) {
	jacobi[j*3+i] += local_nodes[k].x[i]*shape_func_n_d_g(k,j,xi,eta,zeta);
      }
    }
  }
  return (jacobi);
}

/*! \brief creates the strain stress relation matrix for an element.
 * This function creates the strain stress relationship matrix D according to 
 * the three dimensional hooks law for one element.
 */
double* create_strain_stress(element ele) {

  int i,j ;
  
  double prefix;
  
  double nu = ele.poisson_number;
  double E = ele.E_modulus;
  double* d = (double*) malloc (sizeof(double)*36);
  if (d == NULL) {
    printf("Could not allocate memory for stress strain matrix \n");
    printf("- create_stress_strain() \n");
  }
  
  prefix = 1./E;

  // zero out all elements;
  memset(d,0,36*sizeof(double));
  
  for(i=0;i<3;i++) {
    // of diagnoal elements
    for(j=0;j<3;j++) {
      d[j*6+i] = -1.*nu*prefix;
    }
    // diagonal elements overwrites previous nu setting for upper block
    d[i*6+i] = prefix;
    d[(3+i)*6+3+i] = 2*prefix*(1+nu);
  }
  
  return d;
}

/*! \brief creates the stress strain relation matrix for an element.
 * This function creates the stress strain relationship matrix according to 
 * the three dimensional hooks law for one element.
 */
double* create_stress_strain(element ele) {

  int i,j ;
  
  double prefix;
  
  double nu = ele.poisson_number;
  double E = ele.E_modulus;
  double* d = (double*) malloc (sizeof(double)*36);
  if (d == NULL) {
    printf("Could not allocate memory for stress strain matrix \n");
    printf("- create_stress_strain() \n");
  }
  
  prefix = E/((1+nu)*(1-2*nu));

  // zero out all elements;
  memset(d,0,36*sizeof(double));
  
  for(i=0;i<3;i++) {
    // of diagnoal elements
    for(j=0;j<3;j++) {
      d[j*6+i] = prefix*nu;
    }
    // diagonal elements overwrites previous nu setting for upper block
    d[i*6+i] = prefix*(1 - nu);
    d[(3+i)*6+3+i] = prefix*((1-2*nu)*0.5f);
  }

  return d;
}

/*! \brief A function that creates the pseudo strain displacement matrix (G)
 *         for a hexheadral element with bubble functions.
 * This function allows one to obtain the pseudo strain displacement caused
 * by the bubble functions of the element. In the manual denoted as (D) 
 * matrix. 
 * \param localnodes an array containing the 8 nodes of this element.
 * \param jac_null_inv the inverse of the jacobi matrix and hence
 *                     (partial(x,y,z)/partial(xi,eta,zeta))^-1 at 
 *                     xi = eta = zeta = 0.
 * \param jac_null_jac_quot the quotient of the determinants
 *                          of the jacobi matrix at xi, eta, zeta = 0 
 *                          over the determinat of the jacobi matrix 
 *                          at xi = xi, eta = eta, zeta = zeta
 *                          ( in the quadrature points for instance ).
 * \param xi the xi coordinate value.
 * \param eta the eta coordinate value.
 * \param zeta the zeta coordinate value.
 */
double* create_pseudo_strain_displacement(node* localnodes,
					  double* jac_null_inv,
					  double jac_null_jac_quot,
					  double xi, double eta, double zeta) {
  int i,j,i_three;
  double* G = (double*)malloc(sizeof(double)*54);
  double x[3];

  memset(G,0,54*sizeof(double));
  for(i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      x[j] = bubble_func_n_d_l(i,j,jac_null_inv,xi,eta,zeta)
	*jac_null_jac_quot;
    }
    i_three = i*3;

    // 9 elements per row

    G[i_three] = x[0];          // 9
    G[i_three+10] = x[1];       // 18
    G[i_three+20] = x[2];       // 27
    G[i_three+27] = x[1]; G[i_three+28] = x[0]; // 36
    G[i_three+37] = x[2]; G[i_three+38] = x[1]; // 45
    G[i_three+45] = x[2]; G[i_three+47] = x[0];
  }
  return(G);
}
    
    
/*! \brief A function that creates the strain displacement matrix (B)
 *         for a hexheadral element.
 * This function allows one to obtain the strain displacement relation for a
 * hexaheadral element. In the manual denoted as (B) matrix. 
 * \param localnodes an array containing the 8 nodes of this element.
 * \param jacobi the jacobi matrix and hence 
 *                (partial(x,y,z)/partial(xi,eta,zeta)) at xi, eta, zeta.
 * \param xi the xi coordinate value.
 * \param eta the eta coordinate value.
 * \param zeta the zeta coordinate value.
 */
double* create_strain_displacement(node* localnodes, double* jacobi,
				   double xi, double eta, double zeta) {

  int i,j,i_three;

  double* B = (double*)malloc(sizeof(double)*144);
  double x[3]; // derived shape function buffer

  double* invers =
    create_invers_three_by_three(jacobi);

  memset(B,0,144*sizeof(double));
  
  for (i=0;i<8;i++) {
    for (j=0;j<3;j++) {
      x[j] = shape_func_n_d_l(i,j,invers,xi,eta,zeta);
    }
    i_three = i*3;
    
    // 24 elements per row 
    B[i_three] = x[0];    //24
    B[i_three+25] = x[1]; //48
    B[i_three+50] = x[2]; //72
    B[i_three+72] = x[1]; B[i_three+73]=x[0]; //96
    B[i_three+97] = x[2]; B[i_three+98]=x[1]; //120 
    B[i_three+120] = x[2]; B[i_three+122]=x[0];    
  }
  free(invers);
  return(B);
}

/* \brief a function that returns the matrix resulting from the matrix 
 *        multiplication of the stress_strain relation with the
 *        pseudo strain displacement relation (hence D*G).
 * This function performs the simple matrix multiplication of the stress strain
 * relation with the prseudo strain displacement relation and as such returns a
 * pseudo stress displacement relation DG.
 */
double* create_D_G(double* stress_strain, double* pseudo_strain_displacement) {

  int i,j,k;
  
  double * D_G = (double*)malloc(sizeof(double)*54);
  memset(D_G,0,sizeof(double)*54);
  
  for(i=0;i<9;i++) {
    for (j=0;j<6;j++) {
      for (k=0;k<6;k++) {
	D_G[j*9+i] += stress_strain[j*6+k]*pseudo_strain_displacement[k*9+i];
      }
    }
  }
  return(D_G);
}

/* \brief a function that returns the matrix resulting from the matrix 
 *        multiplication of the pseudo strain displacement relation transposed
 *        with the pseudo stress displacement relation DG.
 * This function perfoms the matrix multiplaction of the 
 * pseudo strain displacement relation transposed with the pseudo stress 
 * displacement relation. As such it returns a pseudo stiffness relation,
 * or the lower right part of the stiffness matrix that is going to be 
 * condensed into the matrix later on. 
 * \param pseudo_strain_displacement the pseudo strain displacement matrix.
 * \param D_G the pseudo stress displacement relation.
*/
double* create_GT_D_G(double* pseudo_strain_displacement, double* D_G) {

  int i,j,k;
  
  double * GT_D_G = (double*)malloc(sizeof(double)*81);
  memset(GT_D_G,0,sizeof(double)*81);
  
  for(i=0;i<9;i++) {
    for(j=0;j<9;j++) {
      for(k=0;k<6;k++) {
	GT_D_G[j*9+i] += pseudo_strain_displacement[k*9+j]*D_G[k*9+i];
      }
    }
  }
  return(GT_D_G);
}

/* \brief a function that returns the matrix resulting from the matrix 
 *        multiplication of the transposed strain displacement 
 *        relation with the pseudo stress displacement relation.
 * This function perfoms the matrix multiplication of the 
 * strain displacement relation transposed with the pseudo stress 
 * displacement relation. As such it returns the off diagonal parts of the 
 * stiffness matrix linking the bubble deformation (pseudo part) with the 
 * standard hexaheadral deformations.
 * \param strain_displacement the strain displacement matrix.
 * \param D_G the pseudo stress displacement relation.
 */
double* create_BT_D_G(double* strain_displacement, double* D_G) {

  int i,j,k;
  
  double * BT_D_G = (double*)malloc(sizeof(double)*216);
  memset(BT_D_G,0,sizeof(double)*216);

  for(i=0;i<9;i++) {
    for(j=0;j<24;j++) {
      for(k=0;k<6;k++) {
	BT_D_G[j*9+i] += strain_displacement[k*24+j]*D_G[k*9+i];
      }
    }
  }
  return(BT_D_G);
}

/* \brief A function that returns the matrix resulting from the matrix 
 *        multiplication of the transposed strain displacement 
 *        relation with the stress displacement relation (BT_D_B).
 * This function calculates the stress displacement relation first by
 * multiplying the stress strain relation with the strain displacement
 * relation, then the function performs the matrix multiplication of the 
 * strain displacement relation transposed with the stress 
 * displacement relation. As such it returns the upper left part of the
 * stiffness matrix.
 * \param strain_displacement the strain displacement matrix (B).
 * \param stress_strain the stress strain matrix(D).
 */
double* create_BT_D_B(double* strain_displacement,
		      double* stress_strain) {
  
  int i,j,k;

  double * BT_D_B = (double*)malloc(sizeof(double)*576);
  
  double * D_B  = (double*)malloc(sizeof(double)*144);

  memset(D_B,0,144*sizeof(double));
  memset(BT_D_B,0,576*sizeof(double));
  
  for (i=0;i<24;i++) {
    for (j=0;j<6;j++) {
      for (k=0;k<6;k++) {
	D_B[j*24+i] += stress_strain[j*6+k]*strain_displacement[k*24+i];
      }
    }
  }

  for(i=0;i<24;i++) {
    for(j=0;j<24;j++) {
      for(k=0;k<6;k++) {
	BT_D_B[j*24+i] += strain_displacement[k*24+j]*D_B[k*24+i];
      }
    }
  }

  free(D_B);
  return(BT_D_B);
}

/*! \brief The condensation function that pulls pseudo stiffnesses caused by 
 *         bubble functions into the (upper left) stiffness matrix.
 * A function that pull the pseudo stiffness cause by bubble function into 
 * the stiffness matrix (makes them act on real nodes). This actually is 
 * a kind of gauss elemination where only the lower right part 
 * (where coefficients from (D) are acting) is made diagonal
 * to one modifying the rest of the matrix. Hence condesating them into the
 * upper left part of the matrix.
 */ 
void condensate(double* k) {

  int i,j,key;
  double key_val;
    
  for(key=32;key>23;key--) {
    key_val = k[key*33+key];
    for(i=0;i<=key;i++) {
      k[key*33+i] /= key_val;
      for(j=0;j<=key;j++) {
	k[j*33+i] -= k[j*33+key]*k[key*33+i];
      }
    }
  }
}

/*! \brief A function that creates the stiffness matrix for an element.
 *  This function is the main call to create the stiffness matrix for an 
 *  element. It calls the various functions in order to build the matrix
 *  of a condensated stiffness matrix. 
 *  \param elements_and_nodes the structure structure of the whole finite 
 *                            element system.
 *  \param ele the element to build the stiffness matrix for. 
 */
double* create_stiffness_matrix_for_element(structure elements_and_nodes,
					    element ele) {
  int i,j,n,m,i_three;
  int index;
  
  double* jacobi;
  double* jac_null;
  double* jac_null_inv;
  double* stress_strain;
  double* strain_displacement;
  double* pseudo_strain_displacement;
  double* q_points;

  // strain_displacement_transpose . stress_strain . strain_displacement
  double* BT_D_B;
  double* BT_D_G;
  double* GT_D_G;
  double* D_G;
  
  double xi, eta, zeta, jacobi_det, jac_null_det, jac_null_jac_quot;
  
  int * local_nodes_index;
  node * local_nodes = (node*) malloc(sizeof(node)*8);
  double* k = malloc(sizeof(double)*1089);

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
    strain_displacement = create_strain_displacement(local_nodes,jacobi,
						     xi,eta,zeta);

    jacobi_det = det_three_by_three(jacobi);
    jac_null_jac_quot = jac_null_det/jacobi_det;
    
    pseudo_strain_displacement =
      create_pseudo_strain_displacement(local_nodes,jac_null_inv,
					jac_null_jac_quot,
					xi,eta,zeta);

    BT_D_B = create_BT_D_B(strain_displacement,stress_strain);

    D_G = create_D_G(stress_strain,pseudo_strain_displacement);

    GT_D_G=create_GT_D_G(pseudo_strain_displacement, D_G);

    BT_D_G=create_BT_D_G(strain_displacement,D_G);

    free(strain_displacement);
    free(pseudo_strain_displacement);
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
  
  condensate(k);
  
  free(jac_null_inv);
  free(q_points);
  free(stress_strain);
  return k;
}
