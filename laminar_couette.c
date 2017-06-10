/*=========================================================================
*
*     NAME: couetteFlow.c
* 
*     This file computes the velocity profile of laminar Couette flow 
*     
*     To run, enter in three commands in terminal:
*     >> gcc -c couetteFlow.c
*     >> make couetteFlow
*     >> ./couetteFlow couetteFlow.c
*
*     Contains:
*
*       INPUTS: None
*
*       OUTPUTS: Text file of velocity profiles (coming soon)
*
*       RETURNS:
*   
*  
*     Bryan Kaiser
*     16 November 13
*
*--------------------------------------------------------------------------*/

// Header files

#include <stdio.h>
#include <stdlib.h> // needed for random number generation
#include <math.h> // needed for trig, exp, ceil, floor, etc

/*--------------------------------------------------------------------------*/
main()
{

 // Declare variables
 int i,j,Nx,Ny,iw,t,Nt;
 double rho,nu,mu,dt,dx,dy,dp,pi,Ui,Un,Unp1,Utopw,Ubotw,Umonitor;

 // Couette Flow Assumptions
 // Transient
 // 2D duct
 // Newtonian
 // Isothermal
 // Viscous
 // Negligible body forces
 // No external pressure gradients
 // Laminar

 // Timesteppin'
 dt = 0.00035; // 0.0033s [30x3], 0.001s [60x6], 0.00035s [120x12], timestep (s)
 Nt = 500; // Number of timesteps	

 // Grid Parameters
 iw = 1; // Number of wall nodes per side
 Nx = 120; // Number of cells in x [30] [60] [120]
 Ny = 12; // Number of cells in y [3] [6] [12]
 dx = 0.000833333; // 0.0032 (Nx=30), 0.0017 (Nx=60), 0.00083333 (Nx=120), uniform internal nodes Nx (m)
 dy = 0.00038462; // 0.0013 (Ny=3), 0.00071429 (Ny=6), 0.00038462 (Ny=12), uniform internal nodes Ny (m)
 dp = 0.5/(Nx+2*iw); // Pressure difference, L = 0.1m, dP/dx = 5
 double uVeln[Nx][Ny+(2*iw)][Nt]; // Previous timestep grid, velocity
 double p[Nx][Ny+(2*iw)]; // Previous timestep grid, pressure

 // Material Properties
 mu = 0.0000415; // Dynamic viscosity kg/(ms)
 rho = 0.35; // Density kg/m^3
 nu = mu/rho; // Kinematic viscosity m^2/s

 // Initial and Boundary Condition Specifications
 Ui = 0; // initial velocity in x (streamwise direction)
 pi = 0; // gage initial pressure
 Utopw = 0.889; // velocity for top wall moving 1 m/s
 Ubotw = 0; // velocity of bottom wall moving 0 m/s

 // The Initial Velocity & Pressure Field
 for(i=0; i<(Nx); i++){
   for(j=0; j<(Ny+2*iw); j++){
     uVeln[i][j][0] = Ui; // flow field initialization
     p[i][j] = pi+dp*(i-1); // Switch sign for adverse/favorable gradient
   } // for j
 } // for i

 /*// Initial pressure set-up check (print out)
  for(j=0; j<(Ny+2*iw); j++){
   printf("\n");
    for(i=1; i<(10); i++){
      printf("%2.4f  ", p[i][j][0]);
    } // for i
    //printf("\n");
  } // for j */

 
 // Computing the flow field
 for(t=1; t<Nt; t++){ // Iteration loop
  for(i=1; i<(Nx); i++){ // x direction (streamwise)
   for(j=1; j<(Ny+2*iw); j++){ // y direction (wall normal)
   uVeln[i][Ny+iw][t-1] = Utopw;
   uVeln[i][0][t-1] = Ubotw; 
   uVeln[i][j][t] = uVeln[i][j][t-1]+dt*((nu*(uVeln[i][j+1][t-1]-2*uVeln[i][j][t-1]
             +uVeln[i][j-1][t-1])/pow(dy,2))-(1/rho)*(p[i+1][j]-p[i-1][j])/(2*dx)); 
   } // for y
  } // for x
  
  Un = uVeln[60][7][t-1]; // Saved midpoint velocity value
  Unp1 = uVeln[60][7][t]; // Saved midpoint velocity value
  Umonitor = Unp1-Un; // Convergence monitor
  
  // Printing convergence monitor
  //printf("\n");
  //printf("Time: %f \n",dt*t);
  //printf("Convergence: %2.9f \n",Umonitor);
  printf("%2.9f %f \n",Umonitor,dt*t);

 } // For t, up to final time T

//----------------------------------------------------------------------------------

 /*// Printed data
 // Velocity field print out, first 5 x locations
  for(j=0; j<(Ny+2*iw); j++){
   printf("\n");
    for(i=5; i<(10); i++){
      printf("%2.4f  ", uVeln[i][j][Nt-2]);
    } // for i
    //printf("\n");
  } // for j 
 printf("\n \n"); */

return 0;
}

