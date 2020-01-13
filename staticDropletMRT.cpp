/******************************************************************/
/* Multicomponent, Psuedopotential LBM code with MRT collsion 
/* Static Bubble Test without Gravity /*
/* Copyright (c) Parthib R. Rao; Rice University, Houston TX; parthib.rao@rice.edu
/******************************************************************/	
#include <iostream>
#include <armadillo>
#include <cmath>

using namespace std;
using namespace arma;

int main(int argc, char** argv)
{
/******************************************************************/
			//  Lattice related constants //
/******************************************************************/	
int NV = 9;
vec w; //vector of doubles
w << 4./9. << 1./9. << 1./9. << 1./9. << 1./9. << 1./36. << 1./36. << 1./36. << 1./36.;

ivec cx; // vector of signed integers
ivec cy;
ivec opp;
	   
cx << 0 << 1 << 0 << -1 << 0 << 1 << -1 << -1 << 1; 
cy << 0 << 0 << 1 << 0 << -1 << 1 << 1 << -1 << -1;
opp << 0 << 3 << 4 << 1 << 2 << 7 << 8 << 5 << 6;

double cs_sq = 1.0/3.0; 
double cs_sq_sq= 1.0/9.0;
const double PI  =3.141592653589793238463;
/******************************************************************/

// Geometry parameters
int NY = 50;
int NX = 50;  
int initRadius = 10;
 
/******************************************************************/
		// Fluid and simulation parameters //
/******************************************************************/		
double tau1 = 1.0; 		double omega1 = 1.0/tau1;
double tau2 = 1.0;  	double omega2 = 1.0/tau2;

double rho1 = 3.0; 	//rho1>rho2		      
double rho2 = 1.0; 
double rhoMini = 0.00005;			  

// Pseudopotenial model parameters //
double G11 = 0.0;
double G22 = 0.0;
double G12 = 0.00;

double G1s = 0.0;
double G2s = 0.0;

// MRT parameters // 
vec omega (NV); 
double sc = 1.0; //conserved moments; has to be one
double se= 1.0; // s2 or w_e; related to bulk viscosity
double sepsilon = 1.0; // energy squared; s3 
double sq = 1.0; // heat flux; s5 = s7 or w_q

double snu = 1.0; // viscosity		

omega << sc << se << sepsilon << sc << sq << sc << sq << snu << snu;

// Gravity
float gx = 0.0; 
float gy = 0.0; 

int NT = 50; // Number of iterations

/**********************************************************************/
		// VARIABLE DECLARATION //
/**********************************************************************/
mat rhoCombined (NY,NX);

cube fComp1 = zeros<cube>(NY,NX,NV); 			cube fComp2 = zeros<cube>(NY,NX,NV);
cube fPostCollComp1 = zeros<cube>(NY,NX,NV); 	cube fPostCollComp2 = zeros<cube>(NY,NX,NV);
cube ForceComp1 = zeros <cube> (NY,NX,NV); 		cube ForceComp2 = zeros<cube>(NY,NX,NV);
vec fComp1_temp(NV); 							vec fComp2_temp(NV); 
vec fEqComp1(NV); 								vec fEqComp2(NV);

mat rhoComp1 = zeros(NY,NX);		mat rhoComp2 = zeros(NY,NX);
mat psiComp1 = zeros(NY,NX); 		mat psiComp2 = zeros (NY,NX);
mat pressure = zeros(NY,NX);
vec momentsComp1(NV); 				vec EquiMomentsComp1(NV);
vec momentsComp2(NV); 				vec EquiMomentsComp2(NV);
vec SourceMomentsComp1(NV);			vec SourceMomentsComp2(NV);

mat uxComp1 = zeros(NY,NX); 		mat uyComp1 = zeros(NY,NX);
mat uxComp2 = zeros(NY,NX); 		mat uyComp2 = zeros(NY,NX);
mat uxMix = zeros(NY,NX); 			mat uyMix = zeros(NY,NX);	

mat FXComp1 = zeros(NY,NX);			mat FYComp1 = zeros(NY,NX); 	
mat FXComp2 = zeros(NY,NX); 		mat FYComp2 = zeros(NY,NX);
	 
/**********************************************************************/
	
/**********************************************************************/
//			 INITIALIZATION                               //	
/**********************************************************************/	

// Initialize density on FLUID NODES
for (int i=0; i<NX; ++i) {
	for (int j=0; j<NY; ++j) {
		
		if (((j-NY/2)*(j-NY/2) + ((i-NX/2)*(i-NX/2))) < initRadius*initRadius)	{
			rhoComp1(j,i) = rho1; 	
			rhoComp2(j,i) = rhoMini;
			psiComp1(j,i) = rhoComp1(j,i); 	
			psiComp2(j,i) = rhoComp2(j,i); 
			}
		else 	{
			rhoComp1(j,i) = rhoMini; 
			rhoComp2(j,i) = rho2;
			psiComp1(j,i) = rhoComp1(j,i); 
			psiComp2(j,i) = rhoComp2(j,i); 
			}			
		}
	}
rhoComp1.save("rhoComp1Initial.txt",raw_ascii);
rhoComp2.save("rhoComp2Initial.txt",raw_ascii);

// Initialize distributions assuming zero initial velocity
for (int k=0; k<NV; ++k) {
	for (int i=0; i<NX; ++i) {
		for (int j=0; j<NY; ++j) {		
				fComp1(j,i,k) = rhoComp1(j,i)*w(k);
				fComp2(j,i,k) = rhoComp2(j,i)*w(k); 
				fPostCollComp1(j,i,k) = fComp1(j,i,k);
				fPostCollComp2(j,i,k) = fComp2(j,i,k);
			}
		}
	}
momentsComp1.fill(0.0); momentsComp2.fill(0.0);
EquiMomentsComp1.fill(0.0); EquiMomentsComp2.fill(0.0);
SourceMomentsComp1.fill(0.0); SourceMomentsComp2.fill(0.0);

/**********************************************************************/
		// Start of Main Time Loop //                      
/**********************************************************************/				
for (int iT=0; iT<NT; ++iT)
{ 
		
/**********************************************************************/
		// Density, Pressure and Velocity //                      
/**********************************************************************/
int noOfFluidNodes = 0;
double sum_f1cx, sum_f1cy;
double sum_f2cx, sum_f2cy;
double rhoMix;

for (int i=0; i<NX; ++i) {
	for (int j=0; j<NY; ++j) {

		rhoComp1(j,i) = 0.0; 
		rhoComp2(j,i) = 0.0;
		sum_f1cx = 0.0; sum_f1cy = 0.0;
		sum_f2cx = 0.0; sum_f2cy = 0.0;
	
			for (int k=0; k<NV; ++k)	{
		
				rhoComp1(j,i) = rhoComp1(j,i) + fComp1(j,i,k);
				rhoComp2(j,i) = rhoComp2(j,i) + fComp2(j,i,k);
					 
				sum_f1cx = sum_f1cx + cx(k)*fComp1(j,i,k);
				sum_f1cy = sum_f1cy + cy(k)*fComp1(j,i,k);

				sum_f2cx = sum_f2cx + cx(k)*fComp2(j,i,k);
				sum_f2cy = sum_f2cy + cy(k)*fComp2(j,i,k);
				}

			rhoMix = rhoComp1(j,i) + rhoComp2(j,i);

			// Pseudopotential: psi(x) = rho(x);
			psiComp1(j,i) = rhoComp1(j,i);
			psiComp2(j,i) = rhoComp2(j,i);			
						
			// Pressure				
			pressure(j,i) = cs_sq*(rhoComp1(j,i)+rhoComp2(j,i)) + (0.16667*G12*psiComp1(j,i)*psiComp2(j,i));	
			
			// Component specific bare velocities					
			uxComp1(j,i) = sum_f1cx/rhoComp1(j,i);
			uyComp1(j,i) = sum_f1cy/rhoComp1(j,i);		
			uxComp2(j,i) = sum_f2cx/rhoComp2(j,i); 
			uyComp2(j,i) = sum_f2cy/rhoComp2(j,i);
				
			// Note the difference due to MRT; Mixture velocity that is used for computing EDF
	uxMix(j,i) = ((sum_f1cx + 0.5*FXComp1(j,i))*omega(0) + (sum_f2cx + 0.5*FXComp2(j,i))*omega(0))/(rhoComp1(j,i)*omega(0) + rhoComp2(j,i)*omega(0));
	uyMix(j,i) = ((sum_f1cy + 0.5*FYComp1(j,i))*omega(0) + (sum_f2cy + 0.5*FYComp2(j,i))*omega(0))/(rhoComp1(j,i)*omega(0) + rhoComp2(j,i)*omega(0));

		}
	}

/**********************************************************************/
					// FORCES //	 
/**********************************************************************/
double FXComp1_temp, FXComp2_temp; 
double FYComp1_temp, FYComp2_temp;
double sum_x, sum_y;

for (int i=0; i<NX; ++i) {
	for (int j=0; j<NY;++j) {
		
		FXComp1(j,i) = 0.0;  	FXComp2(j,i) = 0.0;
		FYComp1(j,i) = 0.0;		FYComp2(j,i) = 0.0;		
		FXComp1_temp = 0.0; 	FXComp2_temp = 0.0;
		FYComp1_temp = 0.0;		FYComp2_temp = 0.0;						
		sum_x = 0.0;  
		sum_y = 0.0;
					
		for (int k=0; k<NV;++k)	{
				
			// Find the neighbor of (j,i) in the k-th direction						
			
			int nextI = (i + cx(k) + NX) % NX; 
			int nextJ = (j + cy(k) + NY) % NY;				
					
			// Do the following if neighbor (nextj, nexti) is a fluid node
			FXComp1_temp = FXComp1_temp + w(k)*cx(k)*psiComp1(nextJ,nextI);	
			FYComp1_temp = FYComp1_temp + w(k)*cy(k)*psiComp1(nextJ,nextI);
			
			FXComp2_temp = FXComp2_temp + w(k)*cx(k)*psiComp2(nextJ,nextI);	
			FYComp2_temp = FYComp2_temp + w(k)*cy(k)*psiComp2(nextJ,nextI);				
		} // end for loop for k
							
	// Fluid-Fluid Force (inter + intra component contributions)
	FXComp1(j,i) =  -G11*psiComp1(j,i)*FXComp1_temp - G12*psiComp1(j,i)*FXComp2_temp;
	FYComp1(j,i) =  -G11*psiComp1(j,i)*FYComp1_temp - G12*psiComp1(j,i)*FYComp2_temp;
				
	FXComp2(j,i) =  -G22*psiComp2(j,i)*FXComp2_temp - G12*psiComp2(j,i)*FXComp1_temp;
	FYComp2(j,i) =  -G22*psiComp2(j,i)*FYComp2_temp - G12*psiComp2(j,i)*FYComp1_temp;			

	// Total force on individual components (the last term is due to gravity)
	FXComp1(j,i) = FXComp1(j,i) +  gx*rhoComp1(j,i); //gx*(rhoComp1(j,i)-rhoAve); 
	FYComp1(j,i) = FYComp1(j,i) +  gy*rhoComp1(j,i); // -rhoAve);
	
	FXComp2(j,i) = FXComp2(j,i) + gx*rhoComp2(j,i); //-rhoAve);
	FYComp2(j,i) = FYComp2(j,i) + gy*rhoComp2(j,i); //-rhoAve);
	}
}

/**********************************************************************/	
		// COMBINED MRT COLLISION + FORCING //
/**********************************************************************/
double a,b,c,d,e,h;
a= 1.0/9.0; b=1.0/36.0; c=1.0/18.0;
d=1.0/6.0; e=1.0/4.0; h=1.0/12.0;

double uxSqComp1; double uySqComp1;
double uxSqComp2; double uySqComp2;

// L stands for 'local'
double rhoLComp1; double rhoLComp2; 
double uxLComp1;  double uyLComp1;
double uxLComp2;  double uyLComp2;
double uxMixL; 		double uyMixL;	

for (int i=0; i<NX; ++i) {
	for (int j=0; j<NY; ++j) {
	
		// define local and precomputed variables for purely efficiency purposes
		uxLComp1 = uxComp1(j,i); // Bare velcities, only for meq calc
		uyLComp1 = uyComp1(j,i);
		uxLComp2 = uxComp2(j,i); // Bare velcities, only for meq calc
		uyLComp2 = uyComp2(j,i);
		uxMixL = uxMix(j,i);
		uyMixL = uyMix(j,i);
		
		uxSqComp1 = uxComp1(j,i)*uxComp1(j,i); 	
		uySqComp1 = uyComp1(j,i)*uyComp1(j,i);	
		
		uySqComp2 = uyComp2(j,i)*uyComp2(j,i);
		uxSqComp2 = uxComp2(j,i)*uxComp2(j,i);
		
		rhoLComp1 = rhoComp1(j,i); 
		rhoLComp2 = rhoComp2(j,i);
		
	// Step 2 of MRT collision: Transform to moment space
	momentsComp1(0) = fComp1(j,i,0)+fComp1(j,i,1)+fComp1(j,i,2)+fComp1(j,i,3)+fComp1(j,i,4)+fComp1(j,i,5)+fComp1(j,i,6)+fComp1(j,i,7)+fComp1(j,i,8);
	momentsComp1(1) = -4.0*fComp1(j,i,0)-fComp1(j,i,1)-fComp1(j,i,2)-fComp1(j,i,3)-fComp1(j,i,4)+2.0*(fComp1(j,i,5)+fComp1(j,i,6)+fComp1(j,i,7)+fComp1(j,i,8));
	momentsComp1(2) = 4.0*fComp1(j,i,0)-2.0*(fComp1(j,i,1)+fComp1(j,i,2)+fComp1(j,i,3)+fComp1(j,i,4))+fComp1(j,i,5)+fComp1(j,i,6)+fComp1(j,i,7)+fComp1(j,i,8);
	momentsComp1(3) = fComp1(j,i,1)-fComp1(j,i,3)+fComp1(j,i,5)-fComp1(j,i,6)-fComp1(j,i,7)+fComp1(j,i,8);
	momentsComp1(4) = -2.0*fComp1(j,i,1)+2.0*fComp1(j,i,3)+fComp1(j,i,5)-fComp1(j,i,6)-fComp1(j,i,7)+fComp1(j,i,8);
	momentsComp1(5) = fComp1(j,i,2)-fComp1(j,i,4)+fComp1(j,i,5)+fComp1(j,i,6)-fComp1(j,i,7)-fComp1(j,i,8);
	momentsComp1(6) = -2.0*fComp1(j,i,2)+2.0*fComp1(j,i,4)+fComp1(j,i,5)+fComp1(j,i,6)-fComp1(j,i,7)-fComp1(j,i,8);
	momentsComp1(7) = fComp1(j,i,1)-fComp1(j,i,2)+fComp1(j,i,3)-fComp1(j,i,4);
	momentsComp1(8) = fComp1(j,i,5)-fComp1(j,i,6)+fComp1(j,i,7)-fComp1(j,i,8);

	momentsComp2(0) = fComp2(j,i,0)+fComp2(j,i,1)+fComp2(j,i,2)+fComp2(j,i,3)+fComp2(j,i,4)+fComp2(j,i,5)+fComp2(j,i,6)+fComp2(j,i,7)+fComp2(j,i,8);
	momentsComp2(1) = -4.0*fComp2(j,i,0)-fComp2(j,i,1)-fComp2(j,i,2)-fComp2(j,i,3)-fComp2(j,i,4)+2.0*(fComp2(j,i,5)+fComp2(j,i,6)+fComp2(j,i,7)+fComp2(j,i,8));
	momentsComp2(2) = 4.0*fComp2(j,i,0)-2.0*(fComp2(j,i,1)+fComp2(j,i,2)+fComp2(j,i,3)+fComp2(j,i,4))+fComp2(j,i,5)+fComp2(j,i,6)+fComp2(j,i,7)+fComp2(j,i,8);
	momentsComp2(3) = fComp2(j,i,1)-fComp2(j,i,3)+fComp2(j,i,5)-fComp2(j,i,6)-fComp2(j,i,7)+fComp2(j,i,8);
	momentsComp2(4) = -2.0*fComp2(j,i,1)+2.0*fComp2(j,i,3)+fComp2(j,i,5)-fComp2(j,i,6)-fComp2(j,i,7)+fComp2(j,i,8);
	momentsComp2(5) = fComp2(j,i,2)-fComp2(j,i,4)+fComp2(j,i,5)+fComp2(j,i,6)-fComp2(j,i,7)-fComp2(j,i,8);
	momentsComp2(6) = -2.0*fComp2(j,i,2)+2.0*fComp2(j,i,4)+fComp2(j,i,5)+fComp2(j,i,6)-fComp2(j,i,7)-fComp2(j,i,8);
	momentsComp2(7) = fComp2(j,i,1)-fComp2(j,i,2)+fComp2(j,i,3)-fComp2(j,i,4);
	momentsComp2(8) = fComp2(j,i,5)-fComp2(j,i,6)+fComp2(j,i,7)-fComp2(j,i,8);	 

		// Step 2b: compute equilibrium moments // Per Guo, Luo etc
		EquiMomentsComp1(0) = rhoLComp1;
		EquiMomentsComp2(0) = rhoLComp2; 			
 			
		EquiMomentsComp1(1) = rhoLComp1*(-2.0+3.0*(uxSqComp1+uySqComp1));			
		EquiMomentsComp2(1) = rhoLComp2*(-2.0+3.0*(uxSqComp2+uySqComp2));
		
		EquiMomentsComp1(2) = rhoLComp1*(1.0-3.0*(uxSqComp1+uySqComp1));
		EquiMomentsComp2(2) = rhoLComp2*(1.0-3.0*(uxSqComp2+uySqComp2));
		
		EquiMomentsComp1(3) = rhoLComp1*uxLComp1;
		EquiMomentsComp2(3) = rhoLComp2*uxLComp2;
		
		EquiMomentsComp1(4) = -rhoLComp1*uxLComp1;
		EquiMomentsComp2(4) = -rhoLComp2*uxLComp2;
		 
		EquiMomentsComp1(5) = rhoLComp1*uyLComp1;
		EquiMomentsComp2(6) = -rhoLComp2*uyLComp2;
		
		EquiMomentsComp1(7) = rhoLComp1*(uxSqComp1 - uySqComp1);
		EquiMomentsComp2(7) = rhoLComp2*(uxSqComp2 - uySqComp2);
		
		EquiMomentsComp1(8) = rhoLComp1*uxLComp1*uyLComp1;
		EquiMomentsComp2(8) = rhoLComp2*uxLComp2*uyLComp2;
		
		// Source (Forcing) term moments
		SourceMomentsComp1(0) = 0.0;
		SourceMomentsComp2(0) = 0.0;
		
		SourceMomentsComp1(1) = 6.0*(uxMixL*FXComp1(j,i) + uyMixL*FYComp1(j,i));
		SourceMomentsComp2(1) = 6.0*(uxMixL*FXComp2(j,i) + uyMixL*FYComp2(j,i));
		
		SourceMomentsComp1(2) = -6.0*(uxMixL*FXComp1(j,i) + uyMixL*FYComp1(j,i));
		SourceMomentsComp2(2) = -6.0*(uxMixL*FXComp2(j,i) + uyMixL*FYComp2(j,i));
		
		SourceMomentsComp1(3) = FXComp1(j,i);
		SourceMomentsComp2(3) = FXComp2(j,i);
		
		SourceMomentsComp1(4) = -FXComp1(j,i);
		SourceMomentsComp2(4) = -FXComp2(j,i);
		
		SourceMomentsComp1(5) = FYComp1(j,i);
		SourceMomentsComp2(5) = FYComp2(j,i);
		
		SourceMomentsComp1(6) = -FYComp1(j,i);
		SourceMomentsComp2(6) = -FYComp2(j,i);
		
		SourceMomentsComp1(7) = 2.0*(uxMixL*FXComp1(j,i)-uyMixL*FYComp1(j,i));
		SourceMomentsComp2(7) = 2.0*(uxMixL*FXComp2(j,i)-uyMixL*FYComp2(j,i));
		
		SourceMomentsComp1(8) = uxMixL*FYComp1(j,i) + uyMixL*FXComp1(j,i);
		SourceMomentsComp2(8) = uxMixL*FYComp2(j,i) + uyMixL*FXComp2(j,i);
			
		 // Step 3: Collide: m* = m - omega(m-meq): Post-collisional moments
	//for(int k=0; k<NV; ++k) {
		//momentsComp1(k) = momentsComp1(k) - omega(k)*(momentsComp1(k)-EquiMomentsComp1(k));// + (1.0-0.5*omega(k))*SourceMomentsComp1(k);
		//momentsComp2(k) = momentsComp2(k) - omega(k)*(momentsComp2(k)-EquiMomentsComp2(k));// + (1.0-0.5*omega(k))*SourceMomentsComp2(k);

		//}
			
		// Step 4 Transform from moment space to population space. 
fPostCollComp1(j,i,0) = a*(momentsComp1(0)-momentsComp1(1)+momentsComp1(2));
fPostCollComp2(j,i,0) = a*(momentsComp2(0)-momentsComp2(1)+momentsComp2(2));
		
fPostCollComp1(j,i,1) = a*momentsComp1(0)-b*momentsComp1(1)-c*momentsComp1(2)+d*momentsComp1(3)-d*momentsComp1(4)+e*momentsComp1(7);// + ForceComp1(j,i,1);
fPostCollComp2(j,i,1) = a*momentsComp2(0)-b*momentsComp2(1)-c*momentsComp2(2)+d*momentsComp2(3)-d*momentsComp2(4)+e*momentsComp2(7);// + ForceComp2(j,i,1);

fPostCollComp1(j,i,2) = a*momentsComp1(0)-b*momentsComp1(1)-c*momentsComp1(2)+d*momentsComp1(5)-d*momentsComp1(6)-e*momentsComp1(7);// + ForceComp1(j,i,2);
fPostCollComp2(j,i,2) = a*momentsComp2(0)-b*momentsComp2(1)-c*momentsComp2(2)+d*momentsComp2(5)-d*momentsComp2(6)-e*momentsComp2(7);// + ForceComp2(j,i,2);

fPostCollComp1(j,i,3) = a*momentsComp1(0)-b*momentsComp1(1)-c*momentsComp1(2)-d*momentsComp1(3)+d*momentsComp1(4)+e*momentsComp1(7);// + ForceComp1(j,i,3);
fPostCollComp2(j,i,3) = a*momentsComp2(0)-b*momentsComp2(1)-c*momentsComp2(2)-d*momentsComp2(3)+d*momentsComp2(4)+e*momentsComp2(7);// + ForceComp2(j,i,3);

fPostCollComp1(j,i,4) = a*momentsComp1(0)-b*momentsComp1(1)-c*momentsComp1(2)-d*momentsComp1(5)+d*momentsComp1(6)-e*momentsComp1(7);// + ForceComp1(j,i,4);
fPostCollComp2(j,i,4) = a*momentsComp2(0)-b*momentsComp2(1)-c*momentsComp2(2)-d*momentsComp2(5)+d*momentsComp2(6)-e*momentsComp2(7);// + ForceComp2(j,i,4);

fPostCollComp1(j,i,5) = a*momentsComp1(0)+c*momentsComp1(1)+b*momentsComp1(2)+d*momentsComp1(3)+h*momentsComp1(4)+d*momentsComp1(5)+h*momentsComp1(6)+e*momentsComp1(7);// + ForceComp1(j,i,5);
fPostCollComp2(j,i,5) = a*momentsComp2(0)+c*momentsComp2(1)+b*momentsComp2(2)+d*momentsComp2(3)+h*momentsComp2(4)+d*momentsComp2(5)+h*momentsComp2(6)+e*momentsComp2(7);// + ForceComp2(j,i,5);

fPostCollComp1(j,i,6) = a*momentsComp1(0)+c*momentsComp1(1)+b*momentsComp1(2)-d*momentsComp1(3)-h*momentsComp1(4)+d*momentsComp1(5)+h*momentsComp1(6)-e*momentsComp1(7);// + ForceComp1(j,i,6);
fPostCollComp2(j,i,6) = a*momentsComp1(0)+c*momentsComp2(1)+b*momentsComp2(2)-d*momentsComp2(3)-h*momentsComp2(4)+d*momentsComp2(5)+h*momentsComp2(6)-e*momentsComp2(7);// + ForceComp2(j,i,6);

fPostCollComp1(j,i,7) = a*momentsComp1(0)+c*momentsComp1(1)+b*momentsComp1(2)-d*momentsComp1(3)-h*momentsComp1(4)-d*momentsComp1(5)-h*momentsComp1(6)+e*momentsComp1(7);// + ForceComp1(j,i,7);
fPostCollComp2(j,i,7) = a*momentsComp2(0)+c*momentsComp2(1)+b*momentsComp2(2)-d*momentsComp2(3)-h*momentsComp2(4)-d*momentsComp2(5)-h*momentsComp2(6)+e*momentsComp2(7);// + ForceComp2(j,i,7);

fPostCollComp1(j,i,8) = a*momentsComp1(0)+c*momentsComp1(1)+b*momentsComp1(2)+d*momentsComp1(3)+h*momentsComp1(4)-d*momentsComp1(5)-h*momentsComp1(6)-e*momentsComp1(7);// + ForceComp1(j,i,8);
fPostCollComp2(j,i,8) = a*momentsComp2(0)+c*momentsComp2(1)+b*momentsComp2(2)+d*momentsComp2(3)+h*momentsComp2(4)-d*momentsComp2(5)-h*momentsComp2(6)-e*momentsComp2(7);// + ForceComp2(j,i,8);
			
		}
	}

/**********************************************************************/
				 // STREAMING //
	 //Periodic boundaries in both x and Y directions
/**********************************************************************/
int nextJ, nextI; 
for (int i=0; i<NX; ++i) {
	for (int j=0; j<NY; ++j) 	{
		for (int k=0; k<NV; ++k)	{
			
			nextI = (i + cx(k) + NX) % NX; 
			nextJ = (j + cy(k) + NY) % NY;
			fComp1(nextJ,nextI,k) = fPostCollComp1(j,i,k);
			fComp2(nextJ,nextI,k) = fPostCollComp2(j,i,k);	
		}	 
	}
}

/**********************************************************************/
						// Save data to storage //
/**********************************************************************/
//int ignore;
//ignore = system("mkdir -p output"); // create a folder if it doesnt exist
//imageWriter.writeScaledPpm(rhoComp1,"rhoComp1);
//if (iT%time_save == 0){
	//rhoComp1.save("rhoComp1");
//}


} // End of time loop

/**********************************************************************/
				 // Post-Processing //		
/**********************************************************************/
//rhoComp1.save("rhoComp1Final.txt",raw_ascii);
//rhoComp2.save("rhoComp2Final.txt",raw_ascii);
//pressure.save("pressureFinal.txt",raw_ascii);
rhoComp1.print();

// To find radius of the droplet
double finalRadius;
double rhoAverage = accu(rhoComp1); // 
rhoAverage = rhoAverage/ (NX*NY);
cout<< "Average Density of bubble fluid (fluid1 is  "<< rhoAverage << endl;

double rhoMax = rhoComp1.max();
double rhoMin = rhoComp1.min();

finalRadius = (rhoAverage*NX*NY/(rhoMax*PI) - NX*NY*rhoMin/(rhoMax*PI))/(1.0-(rhoMin/rhoMax));
finalRadius = sqrt(finalRadius);

cout<<"Final radius = "<< finalRadius<<endl;
cout << "Max Density is "<<rhoMax <<"Minimum Density is  "<< rhoMin <<endl;

//// Density profile plot
//vec rhoComp1_profile(NX);
//vec rhoComp2_profile(NX);

//rhoComp1_profile = rhoComp1.col(NX/2);
//rhoComp2_profile = rhoComp2.col(NX/2);

//rhoComp1_profile.save("rhoComp1Profile.txt", raw_ascii);
//rhoComp2_profile.save("rhoComp2Profile.txt",raw_ascii);
//uxMix.save("uxMix.txt",raw_ascii);
//uyMix.save("uyMix.txt",raw_ascii);

//rhoCombined.save("rhoCombinedFinal.txt",raw_ascii);
//pressure.print();

	  return 0;
  }
