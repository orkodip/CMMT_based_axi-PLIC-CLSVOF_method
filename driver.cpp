//SOLUTION OF NON-DIMENSIONAL NAVIER STOKES EQUATION IN AN AXISYMMETRIC FINITE VOLUME DISCRETIZATION TECHNIQUE IN A COLOCATED GRID USING A BALANCED FORCE METHOD.
//THE ADVECTION SCHEME IS ENO2. FDQGMRES IS USED TO SOLVE THE LINEAR SYSTEMS.
#include<iostream>
#include<fstream>
#include<cmath>
#include<string>
#include<csignal>
#include<unistd.h>
#include<ctime>
#define EPS 1e-12
#define cc_NRBS 1e-16
#define TRUNC 1e-5
#define SMALL 1e-8
#define TOL 1e-6	//convergence criteria for FDQGMRES solver
using namespace std;
int LC=1;	//loop condition
void control(int n)	//signal control function
{
	cout<<"PROGRAM INTERUPTED!"<<endl;
	cout<<"Enter 0 to stop: ";
	cin>>LC;
	if(LC==0) cout<<"STOP SIGNAL RECIEVED. PROGRAM WILL STOP AFTER SIMULATING CURRENT TIME-STEP."<<endl;
	else {LC=1; cout<<"SIMULATION CONTINUES."<<endl;}
}
const int I=384,J=768+1;	//no of grids in i and j directions
const double RADIUS=3.0,HEIGHT=6.0;	//domain size
#include "cfd_solvers.h"
#include "mbase.cpp"
#include "GMG2.cpp"
#include "MG_FDQGMRES.cpp"
#include "axi_ivf_ell.cpp"
#include "clsvof.cpp"
#include "NS_g.cpp"
int main()
{
	signal(SIGINT,control);	//define signal and its function (FOR MANUALLY CONTROLLING THE PROGRAM DURING RUNTIME)
	int CNT=0;
	clock_t start=clock();
	NS ms;
	ms.MBASE::ini(CNT,8e-3,0.13625,1.225e-3,8.385e-3,2e-5);	//count,OH,BO,Den_r,Vis_r,Dt
	ms.CLSVOF::ini(4.505,0.5,0.225,2.495,1.5);	//ymc,rm_o,rm_i,yfc,rf
	//ALSO CHANGE THE INITIAL VELOCITY FIELD IN CLSVOF CLASS.
	ms.grid_write();
	ms.lsvf_write(0);
	ms.write_den_vis(0);
	ms.write(0);
	ms.max_GFN();
	for(int COUNT=CNT;(LC && (COUNT<64000));COUNT++)	//manually controlled loop
	{
		ms.CLSVOF::solve(COUNT);	//CLSVOF advection
		ms.NS::solve();
		if((((COUNT+1)%500)==0)||(COUNT==0))
		{
			ms.write(COUNT+1);
			ms.lsvf_write(COUNT+1);
			ms.write_bin(COUNT+1);
		}
	}
	clock_t end=clock();
	cout<<"Simulation Time : "<<(double)(end-start)/CLOCKS_PER_SEC<<" seconds"<<endl;
	return 0;
}
