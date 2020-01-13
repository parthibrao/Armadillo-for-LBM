/******************************************************************/
/* Multicomponent, Multiphase Psuedopotential LBM code with BGK collsion 
/* Static Bubble Test without Gravity /*
/* Version 1 with (NY,NX) instead of (NX,NY); why because
/* (NY,NX) is more intuitive to program/understand than (NX,NY)
/* Copyright (c) Parthib R. Rao; Rice University, Houston TX; parthib.rao@rice.edu
 * Last major update: Aug 22, 2017
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
		// Fluid and simulation parameters //
/******************************************************************/
				
double tau1 = 0.9; 		double omega1 = 1.0/tau1;
double tau2 = 1.0;  	double omega2 = 1.0/tau2;

double rho1 = 5.0; 	//rho1>rho2		      
double rho2 = 1.0; 
double rhoMini = 0.0005;			  

		// Simulation related constants //
double G11 = 0.0;
double G22 = 0.0;
double G12 = 1.35;

double G1s = 0.0;
double G2s = 0.0;

float gx = 0.0; 
float gy = 0.0; // Force in x and y direction	

int NY = 72;
int NX = 72;  
int initRadius = 10;
int time_save = 100;

int NT = 25000; // Number of iterations
/**********************************************************************/
		// VARIABLE DECLARATION //
/**********************************************************************/
mat solid_map = zeros(NY,NX);	
mat rhoCombined (NY,NX);

cube fComp1 = zeros<cube>(NY,NX,NV); 		cube fComp2 = zeros<cube>(NY,NX,NV);
cube fTempComp1 = zeros<cube>(NY,NX,NV); 	cube fTempComp2 = zeros<cube>(NY,NX,NV);
cube SComp1 = zeros <cube> (NY,NX,NV); 		cube SComp2 = zeros<cube>(NY,NX,NV);
vec fComp1_temp(NV); 						vec fComp2_temp(NV); 
vec fEqComp1(NV); 							vec fEqComp2(NV);

mat rhoComp1 = zeros(NY,NX);		mat rhoComp2 = zeros(NY,NX);
mat psiComp1 = zeros(NY,NX); 		mat psiComp2 = zeros (NY,NX);
mat pressure = zeros(NY,NX);

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
			}
		}
	}

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
			
			// Component specific velocities					
			uxComp1(j,i) = sum_f1cx/rhoComp1(j,i);
			uyComp1(j,i) = sum_f1cy/rhoComp1(j,i);		
			uxComp2(j,i) = sum_f2cx/rhoComp2(j,i); 
			uyComp2(j,i) = sum_f2cy/rhoComp2(j,i);
				
			// Mixture velocity that is used for computing EDF
			uxMix(j,i) = (sum_f1cx + 0.5*FXComp1(j,i) + sum_f2cx + 0.5*FXComp2(j,i))/(rhoComp1(j,i) + rhoComp2(j,i));
			uyMix(j,i) = (sum_f1cy + 0.5*FYComp1(j,i) + sum_f2cy + 0.5*FYComp2(j,i))/(rhoComp1(j,i) + rhoComp2(j,i));			
		}
	}

/**********************************************************************/
					// SHAN-CHEN FORCES //	
			//  Only fluid-fluid forces; no fluid-solid forces //
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
	// Compute Source term due to Forces //
/**********************************************************************/
double inv_cs_sq = 1.0/cs_sq; // 1/c_s^2

for (int k=0; k<NV; ++k) {
	for (int i=0; i<NX; ++i) {
		for (int j=0; j<NY; ++j) {
			
	SComp1(j,i,k) = w(k)*FXComp1(j,i)*inv_cs_sq * (cx(k) + inv_cs_sq*(cx(k)*cy(k)*uyMix(j,i)+ uxMix(j,i)*(cx(k)*cx(k)-cs_sq))) + w(k)*FYComp1(j,i)*inv_cs_sq*(cy(k) + inv_cs_sq*(cx(k)*cy(k)*uxMix(j,i) + uyMix(j,i)*(cy(k)*cy(k)-cs_sq)));

	SComp2(j,i,k) = w(k)*cx(k)*FXComp2(j,i)/cs_sq + w(k)*cy(k)*FYComp2(j,i)/cs_sq;
	//SComp2(j,i,k) =  w(k)*FXComp2(j,i)*inv_cs_sq * (cx(k) + inv_cs_sq*(cx(k)*cy(k)*uyMix(j,i)+ uxMix(j,i)*(cx(k)*cx(k)-cs_sq))) + w(k)*FYComp2(j,i)*inv_cs_sq*(cy(k) + inv_cs_sq*(cx(k)*cy(k)*uxMix(j,i) + uyMix(j,i)*(cy(k)*cy(k)-cs_sq)));
		}
	}
}

///**********************************************************************/	
//		// COMBINED COLLISION + FORCING //
///**********************************************************************/
double uMixSq;
vec cdotuMix = zeros<vec>(NV);
	
for (int i=0; i<NX; ++i) {
	for (int j=0; j<NY; ++j) {

		uMixSq = uxMix(j,i)*uxMix(j,i) + uyMix(j,i)*uyMix(j,i);
			
		for(int k=0; k<NV; ++k) {
				
			cdotuMix(k) = cx(k)*uxMix(j,i) + cy(k)*uyMix(j,i);
						
			fEqComp1(k) = w(k)*rhoComp1(j,i) + (w(k)*rhoComp1(j,i) * (3.0*cdotuMix(k) + 4.5*(cdotuMix(k)*cdotuMix(k)) - 1.5*uMixSq));
			fEqComp2(k) = w(k)*rhoComp2(j,i) + (w(k)*rhoComp2(j,i) * (3.0*cdotuMix(k) + 4.5*(cdotuMix(k)*cdotuMix(k)) - 1.5*uMixSq)); 
				
			// Collision for fluid nodes. Now fComp1 and fComp2 are post-collided distributions	
			fComp1(j,i,k) = (1.0-omega1)* fComp1(j,i,k) + omega1*fEqComp1(k) + (1.0 -(0.5*omega1))*SComp1(j,i,k);
			fComp2(j,i,k) = (1.0-omega2)* fComp2(j,i,k) + omega2*fEqComp2(k) + (1.0 -(0.5*omega2))*SComp2(j,i,k);
		} //end of k loop
	}
}
/**********************************************************************/
				 // STREAMING //
		// Periodic in both direcions
/**********************************************************************/
for (int i=0; i<NX; ++i) {
	for (int j=0; j<NY; ++j) 	{
		for (int k=0; k<NV; ++k)	{
			
			int nextI = (i + cx(k) + NX) % NX; 
			int nextJ = (j + cy(k) + NY) % NY;
			fTempComp1(nextJ,nextI,k) = fComp1(j,i,k);
			fTempComp2(nextJ,nextI,k) = fComp2(j,i,k);	
		}	 
	}
}

	// Update distributions for the next timestep
for (int k=0; k<NV; ++k) 	{
	for (int i=0; i<NX; ++i) 	{
		for (int j=0; j<NY; ++j)	{				
				fComp1(j,i,k) = fTempComp1(j,i,k);
				fComp2(j,i,k) = fTempComp2(j,i,k);	
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
rhoComp1.save("rhoComp1Final.txt",raw_ascii);
rhoComp2.save("rhoComp2Final.txt",raw_ascii);
pressure.save("pressureFinal.txt",raw_ascii);

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

// Density profile plot
vec rhoComp1_profile(NX);
vec rhoComp2_profile(NX);

rhoComp1_profile = rhoComp1.col(NX/2);
rhoComp2_profile = rhoComp2.col(NX/2);

//rhoComp1_profile.save("rhoComp1Profile.txt", raw_ascii);
//rhoComp2_profile.save("rhoComp2Profile.txt",raw_ascii);
//uxMix.save("uxMix.txt",raw_ascii);
//uyMix.save("uyMix.txt",raw_ascii);

//rhoCombined.save("rhoCombinedFinal.txt",raw_ascii);
//pressure.print();

	  return 0;
  }
