#include <iostream>
#include <armadillo>
#include <cmath>
using namespace std;
using namespace arma;

int main(int argc, char** argv)	
{
	unsigned int nx, ny, nz;
	nx = 64; ny=64; nz = 48;
	
	/** Create and populate a 3D array named 'NodeID' by reading data from a .dat/.txt file
	 NodeID is a flag variable. For example: if NodeID(j,i,k)==1: fluid if NodeID(j,ik)==0:solid **/
	ucube NodeID; 
	bool ok=NodeID.load("twoSpheres.dat"); // In Matlab: NodeID=importdata(twoSpheres.dat);
		if (ok=false) 	{		
			cout<<"Problem with loading input file"<<endl; 	}		
	NodeID.reshape(ny,nx,nz);  // Matlab: NodeID=reshape(nodeID,[ny,nx,nz]);
	cout<<"Rows: "<<NodeID.n_rows<<" Columns: "<<NodeID.n_cols<<" Slices: "<<NodeID.n_slices<<endl;	
	
	/** Extract linear indices of solid and fluid nodes **/
	uvec solidNodes = find(NodeID==0);  // Matlab: solidNodes = find(NodeID==0);
	cout<<"Number of solid nodes: "<<size(solidNodes)<<endl;
	uvec fluidNodes = find(NodeID==1);  // Matlab: fluidNodes = find(NodeID==1);
  
	/** Create a 'NodeID' type flag variable or simple Cartesian, grid-aligned geometries
	 instead of importing from a input geometry file **/
	ucube NodeID1(ny,nx,nz);	
	NodeID1.fill(1); /** Initialize all elements of NodeID1 array with 1 **/
	/** Set bottom and top walls nodes to zero **/	
	NodeID1(span(0),span::all, span::all)=zeros<ucube>(1,nx,nz); 
	NodeID1(span(ny-1),span::all, span::all)=zeros<ucube>(1,nx,nz);
		
	/** Save fluid variables as a textfile. The output file is a whitespace seperated ASCII 
	 format without metadata or a header **/
	cube ux = randu<cube>(ny,nx,nz);
	/** Extract a slice (2D array/matrix) from a 3D array ux **/
	mat uxOut=ux.slice(nx/2); 
	uxOut.save("ux.txt",raw_ascii); // Matlab: dlmwrite('ux.txt',uxOut,'delimiter',' ');
	cout<<"Max ux is "<<ux.max()<<endl;
	
	/** Submatrix/subcube read/write access of contigous elements for simulation setup, data analysis, etc**/
	mat rho = randu<mat>(ny,nx);
	rho.col(1) = zeros <vec>(nx);            // In Matlab: rho(:,1)=0.0;
	rho(span::all, 1) = zeros<mat>(nx,1);    // rho(:,1)= 0.0;
	mat rhoPart = rho.submat(1,2,3,4);       // rhoPart= rho(1:2,3:4) 
	mat rhoPart1 = rho(span(0,2),span(1,3)); // rhoPart= rho(1:2,3:4) 

	/** Submatrix read/write access of non contigous elements described via linear indices
	/** For example, obtain linear indices of elements of the 'rho' array that are greater than 0.5 **/
	vec rhoHalf = rho.elem(find(rho > 0.5)); // rhoHalf=(rhoHalf>0.5)
	rhoHalf.save("rhoHalf.txt", arma_ascii);
	/** Set ux at all solidNodes to be zero **/
	int nf = solidNodes.n_elem;
	ux(solidNodes) = zeros<vec>(nf); // In Matlab: ux(solidNodes)=0;

	/** set some specific elements of rho to 1 **/
	uvec indices; indices << 2 << 3 << 6 << 8; 	
	int numberOfIndices = indices.n_elem;
	rho.elem(indices) = ones<vec>(numberOfIndices); // rho(indices) = 1.0;
	
	/** Demonstration of Meshgrid **/
	rowvec xSpan(nx); colvec ySpan(ny); 
	mat xCord(ny,nx); mat yCord(ny,nx);
	xSpan = linspace<rowvec>(0,1,nx); // In Matlab: xSpan=linspace(0,1,nx)
	ySpan = linspace<colvec>(0,1,ny); // ySpan = linspace(0,1,ny)

	/** Create xCord matrix from xSpan vector **/
	for (int i=0;i<nx;i++) 	{
			xCord.row(i) = xSpan; }	// In Matlab: [xCord,yCord]=meshgrid(xSpan,ySpan)
	/** Create yCord matrix from ySpan vector **/
	for (int i=0;i<ny;i++) 	{
			yCord.col(i) = ySpan; }

return 0;
}
