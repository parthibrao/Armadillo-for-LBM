/******************************************************************/
/* Singlephase LBM code for 3D Lid-driven cavity
/* Copyright (c) Parthib R. Rao; 
/* Rice University, Houston TX; parthib.rao@rice.edu
/* The following code is for demonstration of Armadillo's capabilties 
/* for LBM prototyping only.
/******************************************************************/	
#include <iostream>
#include <armadillo>
#include <cmath>
#include <fstream>
#include <cstdlib>
using namespace std;
using namespace arma;

// DOMAIN RELATED
#define NX 32 
#define NY 32
#define NZ 32
#define NY1 (NY+1)
#define NX1 (NX+1)
#define NZ1 (NZ+1)
#define L (NY+1)
#define NV 19	
#define rho0 1.0
#define uw 0.01 // Initial velocity of the lid
#define Re 	100
#define maxT 3000

              //1,2,3,4,5,6, 7,8,9,10,11,12,13,14,15,16,17,18
//int opp (NV)={0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15};

cube f0 = zeros<cube>(NY1,NX1,NZ1);		cube fPostColl0 = zeros<cube>(NY1,NX1,NZ1);
cube f1 = zeros<cube>(NY1,NX1,NZ1);		cube fPostColl1 = zeros<cube>(NY1,NX1,NZ1);
cube f2 = zeros<cube>(NY1,NX1,NZ1);		cube fPostColl2 = zeros<cube>(NY1,NX1,NZ1);
cube f3 = zeros<cube>(NY1,NX1,NZ1);		cube fPostColl3 = zeros<cube>(NY1,NX1,NZ1);		
cube f4 = zeros<cube>(NY1,NX1,NZ1);		cube fPostColl4 = zeros<cube>(NY1,NX1,NZ1);		
cube f5 = zeros<cube>(NY1,NX1,NZ1);		cube fPostColl5 = zeros<cube>(NY1,NX1,NZ1);		
cube f6 = zeros<cube>(NY1,NX1,NZ1);		cube fPostColl6 = zeros<cube>(NY1,NX1,NZ1);		
cube f7 = zeros<cube>(NY1,NX1,NZ1);		cube fPostColl7 = zeros<cube>(NY1,NX1,NZ1);		
cube f8 = zeros<cube>(NY1,NX1,NZ1);		cube fPostColl8 = zeros<cube>(NY1,NX1,NZ1);	
cube f9 = zeros<cube>(NY1,NX1,NZ1);		cube fPostColl9 = zeros<cube>(NY1,NX1,NZ1);	
cube f10 = zeros<cube>(NY1,NX1,NZ1);	cube fPostColl10 = zeros<cube>(NY1,NX1,NZ1);	
cube f11 = zeros<cube>(NY1,NX1,NZ1);	cube fPostColl11 = zeros<cube>(NY1,NX1,NZ1);	
cube f12 = zeros<cube>(NY1,NX1,NZ1);	cube fPostColl12 = zeros<cube>(NY1,NX1,NZ1);	
cube f13 = zeros<cube>(NY1,NX1,NZ1);	cube fPostColl13 = zeros<cube>(NY1,NX1,NZ1);	
cube f14 = zeros<cube>(NY1,NX1,NZ1);	cube fPostColl14 = zeros<cube>(NY1,NX1,NZ1);	
cube f15 = zeros<cube>(NY1,NX1,NZ1);	cube fPostColl15 = zeros<cube>(NY1,NX1,NZ1);	
cube f16 = zeros<cube>(NY1,NX1,NZ1);	cube fPostColl16 = zeros<cube>(NY1,NX1,NZ1);	
cube f17 = zeros<cube>(NY1,NX1,NZ1);	cube fPostColl17 = zeros<cube>(NY1,NX1,NZ1);	
cube f18 = zeros<cube>(NY1,NX1,NZ1);	cube fPostColl18 = zeros<cube>(NY1,NX1,NZ1);	 

vec cx;
vec cy;
vec cz;
vec weights; 

ivec opp;
	
vec FEQ(NV); //

cube rho = zeros<cube>(NY1,NX1,NZ1);
cube uxNS = zeros<cube>(NY1,NX1,NZ1);
cube uyNS = zeros<cube>(NY1,NX1,NZ1);
cube uzNS = zeros<cube>(NY1,NX1,NZ1);
cube uNS = zeros<cube>(NY1,NX1,NZ1);

double ux0=0.0; 
double uy0=0.0;
double uz0=0.0;

double u0=0.0; 
double v0=0.0;
double w0=0.0;  
 
double tau;

ivec jd(NV);
ivec id(NV);
ivec kd(NV);

/******************************************************************/
// Function Declarations
void InitAtEqui(void);
double fEq(double RHO, double U, double V, double W, int m);
void BGKCollision(void);
void Streaming(void);
void FluidVars(void);
void BounceBack(void);
/******************************************************************/
int main()
{
	 		
cx<<0<<1<<-1<<0<<0<<0<<0<<1<<1<<-1<<-1<<1<<-1<<1<<-1<<0<<0<<0<<0; 
cy<<0<<0<<0<<1<<-1<<0<<0<<1<<-1<<1<<-1<<0<<0<<0<<0<<1<<1<<-1<<-1; 
cz<<0<<0<<0<<0<<0<<1<<-1<<0<<0<<0<<0<<1<<1<<-1<<-1<<1<<-1<<1<<-1;
weights<<1./3.<<1./18.<<1./18.<<1./18.<<1./18.<<1./18.<<1./18.<<1./36.<<1./36.<<1./36.<<1./36.<<1./36.<<1./36.<<1./36.<<1./36.<<1./36.<<1./36.<<1./36.<<1./36.;
	
	tau = 1.0;
	InitAtEqui();
	int centerY = NY/2; 
	int centerX = NX/2;
	int centerZ = NZ/2;

	for (int iT=0;iT<maxT;iT++) 
		{
			BGKCollision();
			Streaming();
			BounceBack();
			FluidVars();
	
	if (iT%100==0)
		{
			cout<< "Centre X-Velocity "<<uxNS(centerY,centerX,centerZ)<<endl;
			cout<< "Centre Y-Velocity "<<uyNS(centerY,centerX,centerZ)<<endl;
			cout<< "Centre Z-Velocity "<<uzNS(centerY,centerX,centerZ)<<endl;
			cout<< iT<<endl;
		}
	} // End of time loop

	//uxNS.print();
	uNS.save("uNS.txt",raw_ascii);
	
return 0;
} // End of main()
/**********************************************************************/

/**********************************************************************/
/**			 INITIALIZATION                               			 **/	
/**********************************************************************/
void InitAtEqui()
{
	for (int k=0; k<=NZ; k++) {
		for (int i=0; i<=NX; i++) {
			for (int j=0; j<=NY; j++) {						
				rho(j,i,k) = rho0;
				uxNS(j,i,k) = ux0;
				uyNS(j,i,k) = uy0;
				uzNS(j,i,k) = uz0;	
					for (int m=0; m<NV; m++) {
						FEQ(m) = fEq(rho(j,i,k),uxNS(j,i,k),uyNS(j,i,k),uzNS(j,i,k),m);
					}
				f0(j,i,k) = FEQ(0); f1(j,i,k) = FEQ(1); f2(j,i,k) = FEQ(2);
				f3(j,i,k) = FEQ(3);	f4(j,i,k) = FEQ(4); f5(j,i,k) = FEQ(5);
				f6(j,i,k) = FEQ(6); f7(j,i,k) = FEQ(7); f8(j,i,k) = FEQ(8);
				f9(j,i,k) = FEQ(9);	f10(j,i,k) = FEQ(10); f11(j,i,k) = FEQ(11);
				f12(j,i,k) = FEQ(12); f13(j,i,k) = FEQ(13); f14(j,i,k) = FEQ(14);
				f15(j,i,k) = FEQ(15); f16(j,i,k) = FEQ(16); f17(j,i,k) = FEQ(17);
				f18(j,i,k) = FEQ(18);
			}
		}
	}
}

/**********************************************************************/
/**			 COMPUTE EQUILIBRIUM                          			 **/	
/**********************************************************************/
double fEq (double RHO, double U, double V, double W, int m)
{
	double cdotu;
	double uSq;
	cdotu = cx(m)*U+cy(m)*V+cz(m)*W;
	uSq = U*U+V*V+W*W;
	return weights(m)*RHO*(1.0+3.0*cdotu+4.5*cdotu*cdotu-1.5*uSq);
}

/**********************************************************************/
/**			 BGK COLLISION                          			 	**/	
/**********************************************************************/	
void BGKCollision()
{
	for (int k=0; k<=NZ; k++) {
		for (int i=0; i<=NX; i++) {
			for (int j=0; j<=NY; j++) {			
				for (int m=0; m<NV; m++) {
					FEQ(m) = fEq(rho(j,i,k),uxNS(j,i,k),uyNS(j,i,k),uzNS(j,i,k),m);
					}
				fPostColl0(j,i,k) = f0(j,i,k) - (f0(j,i,k)-FEQ(0))/tau;
				fPostColl1(j,i,k) = f1(j,i,k) - (f1(j,i,k)-FEQ(1))/tau;
				fPostColl2(j,i,k) = f2(j,i,k) - (f2(j,i,k)-FEQ(2))/tau;
				fPostColl3(j,i,k) = f3(j,i,k) - (f3(j,i,k)-FEQ(3))/tau;
				fPostColl4(j,i,k) = f4(j,i,k) - (f4(j,i,k)-FEQ(4))/tau;
				fPostColl5(j,i,k) = f5(j,i,k) - (f5(j,i,k)-FEQ(5))/tau;
				fPostColl6(j,i,k) = f6(j,i,k) - (f6(j,i,k)-FEQ(6))/tau;
				fPostColl7(j,i,k) = f7(j,i,k) - (f7(j,i,k)-FEQ(7))/tau;
				fPostColl8(j,i,k) = f8(j,i,k) - (f8(j,i,k)-FEQ(8))/tau;
				fPostColl9(j,i,k) = f9(j,i,k) - (f9(j,i,k)-FEQ(9))/tau;
				fPostColl10(j,i,k) = f10(j,i,k) - (f10(j,i,k)-FEQ(10))/tau;
				fPostColl11(j,i,k) = f11(j,i,k) - (f11(j,i,k)-FEQ(11))/tau;
				fPostColl12(j,i,k) = f12(j,i,k) - (f12(j,i,k)-FEQ(12))/tau;
				fPostColl13(j,i,k) = f13(j,i,k) - (f13(j,i,k)-FEQ(13))/tau;
				fPostColl14(j,i,k) = f14(j,i,k) - (f14(j,i,k)-FEQ(14))/tau;
				fPostColl15(j,i,k) = f15(j,i,k) - (f15(j,i,k)-FEQ(15))/tau;
				fPostColl16(j,i,k) = f16(j,i,k) - (f16(j,i,k)-FEQ(16))/tau;
				fPostColl17(j,i,k) = f17(j,i,k) - (f17(j,i,k)-FEQ(17))/tau;
				fPostColl18(j,i,k) = f18(j,i,k) - (f18(j,i,k)-FEQ(18))/tau;
				}
			}
		}
	}
/**********************************************************************/
				 /** STREAMING **/
/**********************************************************************/
void Streaming()
{
	for (int k=0;k<=NZ;k++)	{
		for (int i=0;i<=NX;i++)	{
			for (int j=0;j<=NY;j++)	{
				for (int m=0;m<NV;m++)	{
					jd(m)=j-cy(m); id(m)=i-cx(m); kd(m)=k-cz(m); // indices of upwind nodes
						}					 
					if (jd(0)>=0 && jd(0)<=NY && id(0)>=0 && id(0)<=NX && kd(0)>=0 && kd(0)<=NZ){
							f0(j,i,k) = fPostColl0(jd(0),id(0),kd(0));
						}
					if (jd(1)>=0 && jd(1)<=NY && id(1)>=0 && id(1)<=NX && kd(1)>=0 && kd(1)<=NZ){
							f1(j,i,k) = fPostColl1(jd(1),id(1),kd(1));
						}
					if (jd(2)>=0 && jd(2)<=NY && id(2)>=0 && id(2)<=NX && kd(2)>=0 && kd(2)<=NZ){
							f2(j,i,k) = fPostColl2(jd(2),id(2),kd(2));
						}
					if (jd(3)>=0 && jd(3)<=NY && id(3)>=0 && id(3)<=NX && kd(3)>=0 && kd(3)<=NZ){
							f3(j,i,k) = fPostColl3(jd(3),id(3),kd(3));
						}
					if (jd(4)>=0 && jd(4)<=NY && id(4)>=0 && id(4)<=NX && kd(4)>=0 && kd(4)<=NZ){
							f4(j,i,k) = fPostColl4(jd(4),id(4),kd(4));
						}
					if (jd(5)>=0 && jd(5)<=NY && id(5)>=0 && id(5)<=NX && kd(5)>=0 && kd(5)<=NZ){
							f5(j,i,k) = fPostColl5(jd(5),id(5),kd(5));
						}
					if (jd(6)>=0 && jd(6)<=NY && id(6)>=0 && id(6)<=NX && kd(6)>=0 && kd(6)<=NZ){
							f6(j,i,k) = fPostColl6(jd(6),id(6),kd(6));
						}
					if (jd(7)>=0 && jd(7)<=NY && id(7)>=0 && id(7)<=NX && kd(7)>=0 && kd(7)<=NZ){
							f7(j,i,k) = fPostColl7(jd(7),id(7),kd(7));
						}
					if (jd(8)>=0 && jd(8)<=NY && id(8)>=0 && id(8)<=NX && kd(8)>=0 && kd(8)<=NZ){
							f8(j,i,k) = fPostColl8(jd(8),id(8),kd(8));
						}
					if (jd(9)>=0 && jd(9)<=NY && id(9)>=0 && id(9)<=NX && kd(9)>=0 && kd(9)<=NZ){
							f9(j,i,k) = fPostColl9(jd(9),id(9),kd(9));
						}
					if (jd(10)>=0 && jd(10)<=NY && id(10)>=0 && id(10)<=NX && kd(10)>=0 && kd(10)<=NZ){
							f10(j,i,k) = fPostColl10(jd(10),id(10),kd(10));
						}
					if (jd(11)>=0 && jd(11)<=NY && id(11)>=0 && id(11)<=NX && kd(11)>=0 && kd(11)<=NZ){
							f11(j,i,k) = fPostColl11(jd(11),id(11),kd(11));
						}
					if (jd(12)>=0 && jd(12)<=NY && id(12)>=0 && id(12)<=NX && kd(12)>=0 && kd(12)<=NZ){
							f12(j,i,k) = fPostColl12(jd(12),id(12),kd(12));
						}
					if (jd(13)>=0 && jd(13)<=NY && id(13)>=0 && id(13)<=NX && kd(13)>=0 && kd(13)<=NZ){
							f13(j,i,k) = fPostColl13(jd(13),id(13),kd(13));
						}
					if (jd(14)>=0 && jd(14)<=NY && id(14)>=0 && id(14)<=NX && kd(14)>=0 && kd(14)<=NZ){
							f14(j,i,k) = fPostColl14(jd(14),id(14),kd(14));
						}
					if (jd(15)>=0 && jd(15)<=NY && id(15)>=0 && id(15)<=NX && kd(15)>=0 && kd(15)<=NZ){
							f15(j,i,k) = fPostColl15(jd(15),id(15),kd(15));
						}
					if (jd(16)>=0 && jd(16)<=NY && id(16)>=0 && id(16)<=NX && kd(16)>=0 && kd(16)<=NZ){
							f16(j,i,k) = fPostColl16(jd(16),id(16),kd(16));
						}
					if (jd(17)>=0 && jd(17)<=NY && id(17)>=0 && id(17)<=NX && kd(17)>=0 && kd(17)<=NZ){
							f17(j,i,k) = fPostColl17(jd(17),id(17),kd(17));
						}
					if (jd(18)>=0 && jd(18)<=NY && id(18)>=0 && id(18)<=NX && kd(18)>=0 && kd(18)<=NZ){
							f18(j,i,k) = fPostColl18(jd(18),id(18),kd(18));				
						}						
				}
			}
		}
}
/**********************************************************************/
				 ///** BOUNCEBACK **/
/**********************************************************************/
void BounceBack()
{
	int i,j,k;
	//k=NZ: Top plate moves in the positive x-direction
	for (i=0;i<=NX;i++)	{
		for (j=0;j<=NY;j++)	{		
			f6(j,i,NZ) = fPostColl5(j,i,NZ);
			f13(j,i,NZ) = fPostColl12(j,i,NZ) - 6*rho(j,i,NZ)*weights(13)*cx(13)*uw;
			f14(j,i,NZ) = fPostColl11(j,i,NZ) - 6*rho(j,i,NZ)*weights(14)*cx(14)*uw;
			f16(j,i,NZ) = fPostColl17(j,i,NZ) - 6*rho(j,i,NZ)*weights(16)*cx(16)*uw;
			f18(j,i,NZ) = fPostColl15(j,i,NZ) - 6*rho(j,i,NZ)*weights(18)*cx(18)*uw;
		}
	}
	
	//k=0: Stationary Bottom wall
	for (j=0;j<=NY;j++)	{
		for (i=0;i<=NX;i++)	{
			f5(j,i,0) = fPostColl6(j,i,0);
			f12(j,i,0) = fPostColl13(j,i,0);
			f11(j,i,0) = fPostColl14(j,i,0);
			f17(j,i,0) = fPostColl16(j,i,0);
			f15(j,i,0) = fPostColl18(j,i,0);
		}
	}
	
	// i=0; Left wall
	for (j=0;j<=NY;j++)	{
		for (k=0;k<=NZ;k++)	{
			f1(j,0,k) = fPostColl2(j,0,k);
			f7(j,0,k) = fPostColl10(j,0,k);
			f8(j,0,k) = fPostColl9(j,0,k);
			f11(j,0,k) = fPostColl14(j,0,k);
			f13(j,0,k) = fPostColl12(j,0,k);
			}
		}
	
	// i=NX; Right wall
	for (j=0;j<=NY;j++)	{
		for ( k=0;k<=NZ;k++)	{
			f2(j,NX,k) = fPostColl1(j,NX,k);
			f10(j,NX,k) = fPostColl7(j,NX,k);
			f9(j,NX,k) = fPostColl8(j,NX,k);
			f14(j,NX,k) = fPostColl11(j,NX,k);
			f12(j,NX,k) = fPostColl13(j,NX,k);
			}
		}
	
	// j=0; Front wall
	for (k=0;k<=NZ;k++)	{
		for (i=0;i<=NX;i++)	{
			f3(0,i,k) = fPostColl4(0,i,k);
			f7(0,i,k) = fPostColl10(0,i,k);
			f9(0,i,k) = fPostColl8(0,i,k);
			f15(0,i,k) = fPostColl18(0,i,k);
			f16(0,i,k) = fPostColl17(0,i,k);
			}
		}
	
	// j=NY; Rear wall
	for (k=0;k<=NZ;k++)	{
		for (i=0;i<=NX;i++)	{
			f4(NY,i,k) = fPostColl3(NY,i,k);
			f10(NY,i,k) = fPostColl7(NY,i,k);
			f8(NY,i,k) = fPostColl9(NY,i,k);
			f18(NY,i,k) = fPostColl15(NY,i,k);
			f17(NY,i,k) = fPostColl16(NY,i,k);
			}
		}
	
}
/**********************************************************************/
				/** FLUID VARIABLES **/
/**********************************************************************/
void FluidVars()
{
	double Xmomentum, Ymomentum, Zmomentum;
	
	for (int k=0;k<=NZ;k++)	{
		for (int i=0;i<=NX;i++)	{
			for (int j=0;j<=NY;j++)	{
			
				rho(j,i,k)=0.;
				uxNS(j,i,k)=0.;
				uyNS(j,i,k)=0.;
				uzNS(j,i,k)=0.;
				 
	Xmomentum=Ymomentum=Zmomentum=0.0;
						
	rho(j,i,k) = f0(j,i,k)+f1(j,i,k)+f2(j,i,k)+f3(j,i,k)+f4(j,i,k)+f5(j,i,k)+f6(j,i,k)+f7(j,i,k)+f8(j,i,k)+f9(j,i,k)
				+f10(j,i,k)+f11(j,i,k)+f12(j,i,k)+f13(j,i,k)+f14(j,i,k)+f15(j,i,k)+f16(j,i,k)+f17(j,i,k)+f18(j,i,k);
	Xmomentum=f1(j,i,k)-f2(j,i,k)+f7(j,i,k)+f8(j,i,k)-f9(j,i,k)-f10(j,i,k)+f11(j,i,k)-f12(j,i,k)+f13(j,i,k)-f14(j,i,k);
	Ymomentum=f3(j,i,k)-f4(j,i,k)+f7(j,i,k)-f8(j,i,k)+f9(j,i,k)-f10(j,i,k)+f15(j,i,k)+f16(j,i,k)-f17(j,i,k)-f18(j,i,k);
	Zmomentum=f5(j,i,k)-f6(j,i,k)+f11(j,i,k)+f12(j,i,k)-f13(j,i,k)-f14(j,i,k)+f15(j,i,k)-f16(j,i,k)+f17(j,i,k)-f18(j,i,k);
				uxNS(j,i,k) = Xmomentum/rho(j,i,k);
				uyNS(j,i,k) = Ymomentum/rho(j,i,k);
				uzNS(j,i,k) = Zmomentum/rho(j,i,k);
			uNS(j,i,k) = sqrt(uxNS(j,i,k)*uxNS(j,i,k) + uyNS(j,i,k)*uyNS(j,i,k));
				
			}
		}
	}
}
/**********************************************************************/
