/***********************************************************************
 * File name: spinodalDecomp.cpp
 * Compile this file on a Linux system (assuming Armadillo is installed) using:
 * g++ spinodalDecomp.cpp -O2 -o spinodalDecomp -larmadillo
 * Please refer to Armadillo documentation for installation and API details
 * The ouput can be visualized on a Python terminal by running the companion
 * plotSpinodalDecomp.py file.
 *  *******************************************************************/
 /**********************************************************************
 * This code implements the 2d Shan-Chen singlecomponent, multiphase LBM 
 * model with BGK collision and explicit (Guo) forcing. 
 * A spinodal decomposition problem is simulated using the Shan-Chen model
 * Copyright (c) Parthib R. Rao; 
 * Address: MS-321, Rice University, Houston,TX 77005; 
 * Email: parthib.rao@rice.edu
 **********************************************************************/	
/***********************************************************************
 * This program is free software; you can redistribute it and/or modify 
 * it under the terms of the GNU General Public License as published by 
 * the Free Software Foundation; either version 2 of the License, or 
 * (at your option) any later version. This program is distributed in 
 * the hope that it will be useful,but WITHOUT ANY WARRANTY; without even 
 * the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License along 
 * with this program; if not, write to the 
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 **********************************************************************/	
#include <iostream>
#include <armadillo>
#include <cmath>
using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
	/**  Lattice related constants **/
	int NV = 9;
	vec w; 
	w << 4./9. << 1./9. << 1./9. << 1./9. << 1./9. << 1./36. << 1./36. << 1./36. << 1./36.;
	double cs_sq = 1.0/3.0; 
	double cs_sq_sq= 1.0/9.0;
	double inv_cs_sq = 1.0/cs_sq; // 1/cs^2

	ivec cx;
	ivec cy;
	ivec opp;		   
	cx << 0 << 1 << 0 << -1 << 0 << 1 << -1 << -1 << 1; 
	cy << 0 << 0 << 1 << 0 << -1 << 1 << 1 << -1 << -1;
	opp << 0 << 3 << 4 << 1 << 2 << 7 << 8 << 5 << 6;

	/** Domain Size **/
	int NY = 32;
	int NX = 32; 

	/** Fluid parameters **/
	double tau = 1.0;  
	double omega = 1.0/tau;
	double rho0 = 200.4; 	//rho1>rho2		      

	/** Simulation constants **/
	double Gc = -130.0;
	double psi0 = 4.0;
	double max = 7.0;
	double min = -7.0;
	int NT = 5000; // Number of iterations
	int tDisp = 100;

	/** Multidimensional array declaration **/
	cube f = zeros<cube>(NY,NX,NV); 		
	cube fTemp = zeros<cube>(NY,NX,NV);
	cube fEq = zeros<cube>(NY,NX,NV);	
	cube S = zeros<cube>(NY,NX,NV); 		

	mat rho = zeros(NY,NX);
	mat psi = zeros(NY,NX); // Interaction potential		
	mat ux = zeros(NY,NX); 		
	mat uy = zeros(NY,NX);
	mat jx = zeros(NY,NX); // x-momentum
	mat jy = zeros(NY,NX); 		
	mat FX = zeros(NY,NX);			
	mat FY = zeros(NY,NX);
	mat cdotu = zeros(NY,NX);	
	vec f_temp(NV); 						
	
	wall_clock timer;
	timer.tic();
	
	/** Initialize a uniformly random distributed density within the domain **/
	mat pert = randu<mat>(NY,NY); // Note randu has range of [0,1]
	pert = pert*(max-min) + min;
	rho= rho0+pert;
	rho.save("rhoInitial.txt",raw_ascii);

	/** Initial condition for both distributions: (t=0) ==> f(i) = w(i) **/ 
	for (int k=0;k<NV;k++){
		f.slice(k)=w(k)*rho;
		fTemp.slice(k)=w(k)*rho; 
		}
	
	/** Start of Main Time Loop **/                   
	for (int iT=0; iT<NT; ++iT)
	{		
		if (iT%tDisp==0)
		{
		//cout<<"For iteration "<<iT<<" the average density is  "<<mean(mean(rho))<<endl;
		}
		/** Compute fluid variables (density, pressure and velocity) **/                     
		rho  = sum(f,2);
		jx.fill(0.); jy.fill(0.);
		for (int k=0;k<NV;k++){
				jx = jx + f.slice(k)*cx(k);
				jy = jy + f.slice(k)*cy(k);			
				}
			psi = psi0*exp(-rho0/rho);		
			ux = (jx + 0.5*FX)/rho; // Force corrected velocity that solves the N-S eqn
			uy = (jy + 0.5*FY)/rho;	
		
	/** Shan-Chen interaction forces (Only fluid-fluid forces, no fluid-solid forces) **/
		FX.fill(0.0); FY.fill(0.0);		
		for (int i=0; i<NX; ++i) {
			for (int j=0; j<NY;++j) {						
				for (int k=0; k<NV;++k)	{					
					int nextI = (i + cx(k) + NX) % NX; 
					int nextJ = (j + cy(k) + NY) % NY;					
					FX(j,i) = FX(j,i) + w(k)*cx(k)*psi(nextJ,nextI);	
					FY(j,i) = FY(j,i) + w(k)*cy(k)*psi(nextJ,nextI);							
						}
					}
				}			
		/** Total fluid-fluid force **/
		FX = -Gc*psi%FX; // Note that % indicate element-wise matrix multiplication and 	
		FY = -Gc*psi%FY; // not the modulo operator

		/** Compute Source term due to interparticle interactions **/
		for (int k=0; k<NV; ++k) {
			// Second-order Explicit Forcing		
			S.slice(k)= w(k)*FX*inv_cs_sq%(cx(k) + inv_cs_sq*(cx(k)*cy(k)*uy+ ux*(cx(k)*cx(k)-cs_sq))) + w(k)*FY*inv_cs_sq%(cy(k) + inv_cs_sq*(cx(k)*cy(k)*ux + uy*(cy(k)*cy(k)-cs_sq)));
			// Optional first-order explicit forcing
			//S.slice(k)= w(k)*cx(k)*FX*inv_cs_sq + w(k)*cy(k)*FY*inv_cs_sq;	
		}
		
		/** Combined forcing and collsion **/
		for(int k=0; k<NV; ++k) {
			cdotu = cx(k)*ux + cy(k)*uy;
			fEq.slice(k) = w(k)*rho + (w(k)*rho % (3.0*cdotu + 4.5*(square(cdotu))- 1.5*(square(ux) + square(uy))));			
			// Collision for fluid nodes. Now f contains post-collision distributions	
			f.slice(k) = (1.0-omega)*f.slice(k) + omega*fEq.slice(k) + (1.0 -(0.5*omega))*S.slice(k);
		} 

		/** Streaming. Periodicity in x-and y-directions is assumed **/
		for (int i=0; i<NX; ++i) {
			for (int j=0; j<NY; ++j) 	{
				for (int k=0; k<NV; ++k)	{			
					int nextI = (i + cx(k) + NX) % NX; 
					int nextJ = (j + cy(k) + NY) % NY;
					fTemp(nextJ,nextI,k) = f(j,i,k);
					}	 
				}
			}
		/** Swap distributions for the next timestep **/
			f = fTemp;	
	} // End of time loop
	double totalTime = timer.toc();
	cout<<endl;
	cout<<"It took "<<totalTime/60.<<" minutes to complete the simulation"<<endl;
	cout<<endl;

	/** Post-Processing **/
	rho.save("rhoFinal.txt",raw_ascii);
	cout << "Maximum Density is "<<rho.max()<<" and Minimum Density is "<< rho.min() <<endl;
	cout<<endl;
	
	double nx1000 = (double)NX/1000.;
	double ny1000 = (double)NY/1000.;
	// Estimate million site updates per second (MSUS)
	double msus = (nx1000*ny1000*(double)NT)/totalTime;
	cout<<"Estimated MSUS for the current hardware is "<<msus<<endl;

return 0;
} // End of main
