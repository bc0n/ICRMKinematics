#include <iostream>
#include <fstream>
#include "kinematicsDLL.h"


int main(int argc, char **argv) {
	printf("Starting runQpsGrid_funxyzpp\n");
	for (int i = 0; i < argc; i++) {
		printf("%d = %s\n", i, argv[i]);
	}


	// read qps file	
	char *fname = argv[1];
	printf("\n\nOpening %s\n", fname);
	FILE *pf = fopen(fname, "r");

	long nrows = 0, ncols = 5, irow = 0;
	double xyzG[3] = { 0,0,0 }, ppG[2] = { 0,0 };
	//fscanf(pf, "%d\n", &nrows);
	//printf("There are %d rows\n", nrows);
	fscanf(pf, "%d,%lf,%lf,%lf,%lf,%lf\n", &nrows, &xyzG[0], &xyzG[1], &xyzG[2], &ppG[0], &ppG[1]); //%lf for long float = double
	printf("There are %d rows with xyzG[%f,%f,%f] ppG[%f,%f]\n", nrows, xyzG[0],xyzG[1],xyzG[2],ppG[0],ppG[1]);

	double **qps = new double*[ncols];
	for (int i = 0; i < ncols; i++) {
		qps[i] = new double[nrows];
	}
	float flt0, flt1, flt2, flt3, flt4;
	for (irow = 0; irow < nrows; irow++) {
		fscanf(pf, "%f,%f,%f,%f,%f\n", &flt0, &flt1, &flt2, &flt3, &flt4);
		qps[0][irow] = (double)flt0;
		qps[1][irow] = (double)flt1;
		qps[2][irow] = (double)flt2;
		qps[3][irow] = (double)flt3;
		qps[4][irow] = (double)flt4;
		//printf("%d: [%f %f %f %f %f]\n", irow, qps[0][irow], qps[1][irow], qps[2][irow], qps[3][irow], qps[4][irow]);
	}
	fclose(pf);
	//irow = 33;
	//printf("%d: [%f %f %f %f %f]\n", irow, qps[0][irow], qps[1][irow], qps[2][irow], qps[3][irow], qps[4][irow]);



	//setup
	double kn11A[11];
	double nlArray[6]; //nonlinear optimization params
	double jArray[10]; //joint limits
	kn11A[0] = 806; kn11A[1] = -66.0; kn11A[2] = -28.0; kn11A[3] = 0; kn11A[4] = -.24; //TxyzRyz01
	kn11A[5] = 8.2; //Tx23
	kn11A[6] = 0; kn11A[7] = 0; kn11A[8] = 98; //Ryz34CathL
	kn11A[9] = 0; kn11A[10] = 0; //Ryz45
	nlArray[0] = 1e6; // maxIts
	nlArray[1] = 60; // max time sec
	nlArray[3] = 1e-9; // min fun val
	nlArray[4] = 1e-11; // tol fun
	nlArray[5] = 1e-9; //tol x
	//nlArray[2] = 00; // GN_DIRECT
	nlArray[2] = 8; // MLSL w/ NM
	//nlArray[2] = 14; //LN_NelderMead
	jArray[0] = -3; // joint minima
	jArray[1] = -.8;
	jArray[2] = -1;
	jArray[3] = 1e-3;
	jArray[4] = 0;
	jArray[5] = 3; // maxima
	jArray[6] = .8;
	jArray[7] = 1;
	jArray[8] = 5;
	jArray[9] = 500;

	//setup file for results
	char resname[100];
	sprintf(resname, "%s_res.csv", fname); printf("Writing to %s\n", resname);
	pf = fopen(resname, "w+");


	//run grid
	int ret = 0;
	double qpsi[5], xyzi[3], ppi[2], fval;
	for (irow = 0; irow < nrows; irow++) {
		for (int i = 0; i < 5; i++) {
			qpsi[i] = qps[i][irow];
		}
		xyzi[0] = xyzG[0]; xyzi[1] = xyzG[1]; xyzi[2] = xyzG[2]; ppi[0] = ppG[0]; ppi[1] = ppG[1]; //set goal
		ret = getFunVal_xyzpp11A(qpsi, kn11A, nlArray, jArray, xyzi, ppi, &fval);
		//            i,q0,q1,q2,q3,q4, x, y, z, p, p, f
		fprintf(pf, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", irow, qpsi[0],qpsi[1],qpsi[2],qpsi[3],qpsi[4], xyzi[0],xyzi[1],xyzi[2], ppi[0],ppi[1], fval);
	}
	fclose(pf);
	
	delete[]qps;
	
	printf("Done, press Enter\n");
	std::getchar();
	
	return 0;

}