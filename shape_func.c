/*! \file shape_func.c
 * \brief implements shape and bubble functions for hexaheadral elements
 * This file contains all the functions to evaluate the shape functions 
 * and bubble functions of hexaheadral finite elements.
 */

#include<stdlib.h>
#include<stdio.h>

/*! \brief function to obtain the value of the n-th hexaheadral shape function
 *         at point xi, eta, zeta
 */
double shape_func_n(int n, double xi,double eta,double zeta) {

  int prefix_table[24] = { -1, -1, -1,
			    1, -1, -1,
			    1,  1, -1,
			   -1,  1, -1,
			   -1, -1,  1,
			    1, -1,  1,
			    1,  1,  1,
			   -1,  1,  1 };

  int* prefix = prefix_table+(n*3);
  
  return (1+xi*prefix[0])*(1+eta*prefix[1])*(1+zeta*prefix[2]);
}

/*! \brief function to obtain the value of the n-th hexaheadral 
 *         shape function at points xi, eta, zeta derived to either 
 *         xi, eta, or zeta
 * This function allows one to obtain the value of the n-th shape function
 * at the values xi, eta, zeta derived to either xi, eta or zeta. For furhter
 * details you might refere to the manual that shall come with this program.
 * n_d_g shall help the programmer to indicate that this function derives to 
 * greek coordiantes.
 * \param n n-th shape function
 * \param g 0: derived to xi, 1: derived to eta, 2: derived to zeta.
 * \param xi xi-coordinate value.
 * \param eta eta-coordinate value.
 * \param zeta zeta-coordinate value.
 */
double shape_func_n_d_g(int n, int g,
			double xi, double eta, double zeta) {

  int prefix_table[24] = { -1, -1, -1,
			    1, -1, -1,
			    1,  1, -1,
			   -1,  1, -1,
			   -1, -1,  1,
			    1, -1,  1,
			    1,  1,  1,
			   -1,  1,  1 };

  int* prefix = prefix_table+(n*3);

  switch (g) {
  case 0:
    return (0.125*prefix[0]*(1+eta*prefix[1])*(1+zeta*prefix[2]));
    break;
  case 1:
    return (0.125*prefix[1]*(1+xi*prefix[0])*(1+zeta*prefix[2]));
    break;
  case 2:
    return (0.125*prefix[2]*(1+xi*prefix[0])*(1+eta*prefix[1]));
    break;
  default:
    printf("Error calling shape_func_n_d_g derival parameter is out of \n");
    printf("range (0,1,2) and is specified to be %d which is invalid \n",g);
    exit(-1);
    break;
  }
}

/*! \brief function to obtain the value of the n-th hexaheadral 
 *         shape function at points xi, eta, zeta derived to either 
 *         x, y, or z
 * This function allows one to obtain the value of the n-th shape function
 * at the values xi, eta, zeta derived to either x, y or z. For furhter
 * details you might refere to the manual that shall come with this program.
 * n_d_l shall help the programmer to indicate that this function derives to 
 * latin coordiantes.
 * \param n n-th shape function
 * \param g 0: derived to x, 1: derived to y, 2: derived to z.
 * \param invers requires the invers jacobi matrix
 * \param xi xi-coordinate value.
 * \param eta eta-coordinate value.
 * \param zeta zeta-coordinate value.
 */
double shape_func_n_d_l(int n, int g,
			double* invers, // requires the invers jacobi matrix
			double xi, double eta, double zeta) {

  return(invers[g*3]*shape_func_n_d_g(n,0,xi,eta,zeta)
	 +invers[g*3+1]*shape_func_n_d_g(n,1,xi,eta,zeta)
	 +invers[g*3+2]*shape_func_n_d_g(n,2,xi,eta,zeta) );
  
}

/*! \brief function to obtain the value of the n-th hexaheadral bubble function
 *         at point xi, eta, zeta
 */
double bubble_func(int n, double xi, double eta, double zeta) {

  double t[3];
  
  t[0] = xi;
  t[1] = eta;
  t[2] = zeta;

  return(1.-t[n]*t[n]);
}


/*! \brief function to obtain the value of the n-th hexaheadral 
 *         bubble function at points xi, eta, zeta derived to either 
 *         xi, eta, or zeta
 * This function allows one to obtain the value of the n-th bubble function
 * at the values xi, eta, zeta derived to either xi, eta or zeta. For furhter
 * details you might refere to the manual that shall come with this program.
 * n_d_g shall help the programmer to indicate that this function derives to 
 * greek coordiantes.
 * \param n n-th bubble function
 * \param g 0: derived to xi, 1: derived to eta, 2: derived to zeta.
 * \param xi xi-coordinate value.
 * \param eta eta-coordinate value.
 * \param zeta zeta-coordinate value.
 */
double bubble_func_n_d_g(int n, int g, double xi, double eta, double zeta) {

  double t[3];
  
  t[0] = xi;
  t[1] = eta;
  t[2] = zeta;
  
  if ( n == g ) {
    return((-2.*t[n]));
  } else {
    return(0.);
  }
}

/*! \brief function to obtain the value of the n-th hexaheadral 
 *         bubble function at points xi, eta, zeta derived to either 
 *         x, y, or z
 * This function allows one to obtain the value of the n-th bubble function
 * at the values xi, eta, zeta derived to either x, y or z. For furhter
 * details you might refere to the manual that shall come with this program.
 * n_d_l shall help the programmer to indicate that this function derives to 
 * latin coordiantes.
 * \param n n-th bubble function
 * \param g 0: derived to x, 1: derived to y, 2: derived to z.
 * \param invers requires the invers jacobi matrix
 * \param xi xi-coordinate value.
 * \param eta eta-coordinate value.
 * \param zeta zeta-coordinate value.
 */
double bubble_func_n_d_l(int n, int g,
			 double* invers, // requires the invers jacobi matrix
			 double xi, double eta, double zeta) {

  return(invers[g*3]*bubble_func_n_d_g(n,0,xi,eta,zeta)
	 +invers[g*3+1]*bubble_func_n_d_g(n,1,xi,eta,zeta)
	 +invers[g*3+2]*bubble_func_n_d_g(n,2,xi,eta,zeta) );
}
  
  
  
