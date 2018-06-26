/******************************************************************/
/* Singlephase LBM code for Poiseuille flow
/* with MRT collision model and second-order (explicit) forcing/*
/* Body force driven flow
/* Copyright (c) Parthib R. Rao; 
/* Rice University, Houston TX; parthib.rao@rice.edu
/* The following code is for demonstration of Armadillo's capabilties 
/* for LBM prototyping only.
/******************************************************************/	
#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{	  
/******************************************************************/
				// Import of solid map and processing //
/******************************************************************/	
//mat temp_solid_map;
//temp_solid_map.load("solid_map.txt");
//uvec domain_size;
//domain_size.load("domain_size.txt");

//int NX = domain_size(1);
//int NY = domain_size(0);  

//mat solid_map = reshape(temp_solid_map,NY,NX);  // check reshape 

int NX = 50;
int NY = 50;

mat solid_map = zeros(NY,NX);
for (int i = 0; i < NX; i++)
	{
		solid_map(0,i) = 1;
		solid_map(NY-1, i) = 1;
	}

/******************************************************************/
			//  Lattice related constants //
/******************************************************************/
int NV = 9;		
double cs_sq = 1.0/3.0; 
double cs_sq_sq = 1.0/9.0;
const double PI = 3.141592653589793238463;

vec w; //weights
w << 4./9. << 1./9. << 1./9. << 1./9. << 1./9. << 1./36. << 1./36. << 1./36. << 1./36.;
ivec ex; // column vector of signed integers
ivec ey;
vec cx; // column vector of floats
vec cy;
ivec opp;
	   
ex << 0 << 1 << 0 << -1 << 0 << 1 << -1 << -1 << 1; 
ey << 0 << 0 << 1 << 0 << -1 << 1 << 1 << -1 << -1;
cx << 0.0 << 1.0 << 0.0 << -1.0 << 0.0 << 1.0 << -1.0 << -1.0 << 1.0; 
cy << 0.0 << 0.0 << 1.0 << 0.0 << -1.0 << 1.0 << 1.0 << -1.0 << -1.0;
opp << 0 << 3 << 4 << 1 << 2 << 7 << 8 << 5 << 6;
/******************************************************************/
  
/******************************************************************/
		// Fluid and simulation parameters //
/******************************************************************/
double rho0 = 1.0; 	
			
// MRT parameters // 
vec omega (NV); 
double sc = 0.0; //s1,s4,s6 for conserved moments; has to be one
double se = 1.0; // 1.63 s2 related to bulk viscosity
double sepsilon = 1.0; // 1.14 s3 energy squared;  
double sq = 1.142857143; //1.92 // s5 = s7; heat flux; or w_q
double snu = 1.0; // s8, s9; viscosity
omega << sc << se << sepsilon << sc << sq << sc << sq << snu << snu;
	  
// Simulation related constants //	
double gx = 0.00002666666667; // Force in x and y direction
double gy = 0.0;
int time_save = 100;

int NT = 10000; // Number of iterations

/**********************************************************************/
		// VARIABLE DECLARATION //
/**********************************************************************/
cube f = zeros<cube>(NY,NX,NV); 		
cube fPostColl = zeros<cube>(NY,NX,NV); 	
cube Force = zeros<cube>(NY,NX,NV); 		
vec f_temp(NV); 												
vec fEq(NV);
mat rho = zeros(NY,NX); 
mat pressure = zeros(NY,NX);	
mat uxNS = zeros(NY,NX); 	mat uyNS = zeros(NY,NX);
mat ux = zeros(NY,NX); 		mat uy = zeros(NY,NX);
mat uNS = zeros(NY,NX);
mat FX = zeros(NY,NX);		mat FY = zeros(NY,NX); 	
vec moments(NV); 	vec EquiMoments(NV); vec SourceMoments(NV);

mat velocity_final = zeros(NY+1,NX+1);
vec Domain = zeros(NY);

for (int i = 0; i < NX; i++)
{
		Domain(i) = i;
}

/**********************************************************************/
	
/**********************************************************************/
//			 INITIALIZATION                               //	
/**********************************************************************/
rho.fill(rho0);	
// Initialize distributions assuming zero initial velocity
for (int k=0; k<NV; ++k) {
	for (int i=0; i<NX; ++i) {
		for (int j=0; j<NY; ++j) {		
				f(j,i,k) = rho(j,i)*w(k);
				fPostColl(j,i,k) = f(j,i,k);
			}
		}
	}

/**********************************************************************/
		// Start of Main Time Loop //                      
/**********************************************************************/				
for (int iT=0; iT<NT; ++iT)
{
	
/**********************************************************************/
					// FORCES //	 
/**********************************************************************/
double sum_x, sum_y;

for (int i=0; i<NX; ++i) {
	for (int j=0; j<NY;++j) {
		if (solid_map(j,i) == 0)	{		
			FX(j,i) = 0.0;  	
			FY(j,i) = 0.0;
									
			// Total force density on the fluid 
			FX(j,i) = gx*rho(j,i); 
			FY(j,i) = gy*rho(j,i); 
		}		
	}
}

/**********************************************************************/	
		// COMBINED MRT COLLISION + FORCING //
/**********************************************************************/
double a,b,c,d,e,h;
a= 1.0/9.0; b=1.0/36.0; c=1.0/18.0;
d=1.0/6.0; e=1.0/4.0; h=1.0/12.0;

double uxSq = 0.0; double uySq = 0.0;
double rhoL = 0.0; // L stands for 'local'
double uxL = 0.0; double uyL = 0.0;	
double uSq;

for (int i=0; i<NX; ++i) {
	for (int j=0; j<NY; ++j) {
		
		// Solid nodes ONLY. Collision (fullway Bounceback) on solid (obstacle) nodes
		if (solid_map(j,i) == 1) { 
			for (int k=0; k<NV; ++k) {					
				f_temp(k) = 0.0;
				f_temp(k) = f(j,i,k);
			}
			for (int k=0; k<NV; ++k) {					
				fPostColl(j,i,opp(k)) = f_temp(k);										
				} 	
			}		
		else if (solid_map(j,i) == 0) { // Fluid nodes
		
			uxSq = uxNS(j,i)*uxNS(j,i); 
			uySq = uyNS(j,i)*uyNS(j,i);
			rhoL = rho(j,i);
			uxL = uxNS(j,i); 
			uyL = uyNS(j,i);
			
		
		// Step 2 of MRT collision: Transform to moment space
		moments(0) = f(j,i,0)+f(j,i,1)+f(j,i,2)+f(j,i,3)+f(j,i,4)+f(j,i,5)+f(j,i,6)+f(j,i,7)+f(j,i,8);
		moments(1) = -4.0*f(j,i,0)-f(j,i,1)-f(j,i,2)-f(j,i,3)-f(j,i,4)+2.0*(f(j,i,5)+f(j,i,6)+f(j,i,7)+f(j,i,8));
		moments(2) = 4.0*f(j,i,0)-2.0*(f(j,i,1)+f(j,i,2)+f(j,i,3)+f(j,i,4))+f(j,i,5)+f(j,i,6)+f(j,i,7)+f(j,i,8);
		moments(3) = f(j,i,1)-f(j,i,3)+f(j,i,5)-f(j,i,6)-f(j,i,7)+f(j,i,8);
		moments(4) = -2.0*f(j,i,1)+2.0*f(j,i,3)+f(j,i,5)-f(j,i,6)-f(j,i,7)+f(j,i,8);
		moments(5) = f(j,i,2)-f(j,i,4)+f(j,i,5)+f(j,i,6)-f(j,i,7)-f(j,i,8);
		moments(6) = -2.0*f(j,i,2)+2.0*f(j,i,4)+f(j,i,5)+f(j,i,6)-f(j,i,7)-f(j,i,8);
		moments(7) = f(j,i,1)-f(j,i,2)+f(j,i,3)-f(j,i,4);
		moments(8) = f(j,i,5)-f(j,i,6)+f(j,i,7)-f(j,i,8);	 

		// Step 2b: compute equilibrium moments // Per Guo, Luo etc
		EquiMoments(0) = rhoL; 			
		EquiMoments(1) = rhoL*(-2.0+3.0*rhoL*(uxSq+uySq));			
		EquiMoments(2) = rhoL*(1.0-3.0*rhoL*(uxSq+uySq));
		EquiMoments(3) = rhoL*uxL;
		EquiMoments(4) = -rhoL*uxL; 
		EquiMoments(5) = rhoL*uyL;
		EquiMoments(6) = -rhoL*uyL;
		EquiMoments(7) = rhoL*rhoL*(uxSq - uySq);
		EquiMoments(8) = rhoL*rhoL*uxL*uyL;
		
		// Source (Forcing) term moments
		SourceMoments(0) = 0.0;	
		SourceMoments(1) = 6.0*(uxL*FX(j,i) + uyL*FY(j,i));		
		SourceMoments(2) = -6.0*(uxL*FX(j,i) + uyL*FY(j,i));	
		SourceMoments(3) = FX(j,i);		
		SourceMoments(4) = -FX(j,i);		
		SourceMoments(5) = FY(j,i);		
		SourceMoments(6) = -FY(j,i);		
		SourceMoments(7) = 2.0*(uxL*FX(j,i)-uyL*FY(j,i));		
		SourceMoments(8) = uxL*FY(j,i) + uyL*FX(j,i);
			
		// Step 3: Collide: m* = m - omega(m-meq): Post-collisional moments
		for(int k=0; k<NV; ++k) {
			moments(k) = moments(k) - omega(k)*(moments(k)-EquiMoments(k)) + (1.0-0.5*omega(k))*SourceMoments(k);
			//moments(k) = moments(k) + (1.0-0.5*omega(k))*SourceMoments(k);
			}
			
	// Step 4 Transform from moment space to population space. 
		fPostColl(j,i,0) = a*(moments(0)-moments(1)+moments(2));		
		fPostColl(j,i,1) = a*moments(0)-b*moments(1)-c*moments(2)+d*moments(3)-d*moments(4)+e*moments(7);
		fPostColl(j,i,2) = a*moments(0)-b*moments(1)-c*moments(2)+d*moments(5)-d*moments(6)-e*moments(7);
		fPostColl(j,i,3) = a*moments(0)-b*moments(1)-c*moments(2)-d*moments(3)+d*moments(4)+e*moments(7);
		fPostColl(j,i,4) = a*moments(0)-b*moments(1)-c*moments(2)-d*moments(5)+d*moments(6)-e*moments(7);
		fPostColl(j,i,5) = a*moments(0)+c*moments(1)+b*moments(2)+d*moments(3)+h*moments(4)+d*moments(5)+h*moments(6)+e*moments(8);
		fPostColl(j,i,6) = a*moments(0)+c*moments(1)+b*moments(2)-d*moments(3)-h*moments(4)+d*moments(5)+h*moments(6)-e*moments(8);
		fPostColl(j,i,7) = a*moments(0)+c*moments(1)+b*moments(2)-d*moments(3)-h*moments(4)-d*moments(5)-h*moments(6)+e*moments(8);
		fPostColl(j,i,8) = a*moments(0)+c*moments(1)+b*moments(2)+d*moments(3)+h*moments(4)-d*moments(5)-h*moments(6)-e*moments(8);

			}// end of if_solid			
		}
	}
/**********************************************************************/
				 // STREAMING //
/**********************************************************************/
int nextI,nextJ;
	for (int i=0; i<NX; ++i) {
		for (int j=0; j<NY; ++j) 	{
				for (int k=0; k<NV; ++k)	{				
					nextI = (i + ex(k) + NX) % NX; 
					nextJ = (j + ey(k) + NY) % NY;				
					f(nextJ,nextI,k) = fPostColl(j,i,k);
				}			 
		}
	}
// Periodic boundary conditions along X -direction
// i = 0: Left boundary nodes
	for (int j=0; j<=NY-1; ++j)
		{
			f(j,0,1) = fPostColl(j,NX-1,1);
			f(j,0,5) = fPostColl(j,NX-1,5);
			f(j,0,8) = fPostColl(j,NX-1,8);
		}
	
// i = NX-1: Right boundary nodes
	for (int j=0; j<=NY-1;++j)
		{
			f(j,NX-1,3) = fPostColl(j,0,3);
			f(j,NX-1,7) = fPostColl(j,0,7);
			f(j,NX-1,6) = fPostColl(j,0,6);
		}
/**********************************************************************/
				// FLUID VARIABLES //                      
/**********************************************************************/
double sum_f_cx, sum_f_cy;

for (int i=0; i<NX; ++i) {
	for (int j=0; j<NY; ++j) {
		
		rho(j,i) = 0.0; 
		sum_f_cx = 0.0; sum_f_cy = 0.0;
		uNS(j,i)=0.0; 
		uxNS(j,i)=0.0; uyNS(j,i)=0.0;
		ux(j,i)=0.0; uy(j,i)=0.0;
			
		if (solid_map(j,i) == 0) { // Compute density, etc only for fluid nodes		
			for (int k=0; k<NV; ++k)	{			
				rho(j,i) = rho(j,i) + f(j,i,k);					 
				sum_f_cx = sum_f_cx + cx(k)*f(j,i,k);
				sum_f_cy = sum_f_cy + cy(k)*f(j,i,k);
				}
					
			// Pressure				
			pressure(j,i) = cs_sq*rho(j,i);	

			// Force corrected velocity that is used for computing EDF; NS=Navier-Stokes
			uxNS(j,i) = (sum_f_cx + 0.5*FX(j,i))/rho(j,i);
			uyNS(j,i) = (sum_f_cy + 0.5*FY(j,i))/rho(j,i);
			uNS(j,i) = sqrt(uxNS(j,i)*uxNS(j,i) + uyNS(j,i)*uyNS(j,i));
			}			
		}
	}
} // End of time loop

/**********************************************************************/
				 // Post-Processing //		
/**********************************************************************/
	//uNS.print();
	//uNS.save("velocity.txt",raw_ascii);
	velocity_final = join_horiz(Domain,uNS);
	velocity_final.save("velocityMRT.txt",raw_ascii);
return 0;
}
