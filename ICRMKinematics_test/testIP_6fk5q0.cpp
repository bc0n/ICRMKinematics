#include <fstream>
#include <iostream>
#include "kinematicsDLL.h"
#include "ip_nlopt.h"




void main() {
	//load test file
	char fname[] = "testSquareQps.dat";
	//char fname[] = "D:\\20160401_inter2d_estimate\\02_optimalHome\\testSquareQps.dat";
	int ncols = 38; //it, q0,q1,q2,q3,q4, H1, H2

	//open
	printf("\n\nOpening %s\n", fname);
	std::fstream fs;
	fs.open(fname, std::fstream::in | std::fstream::binary); //name, read/write mode|binary
	if (!fs) { std::cerr << "warning, could not open " << fname << std::endl; return; };

	//fileSize & initialize variables
	fs.seekg(0, fs.end);
	long fileSize = (long)fs.tellg();
	long nrows = fileSize / sizeof(double) / ncols;
	fs.seekg(0, fs.beg);
//nrows = 1;
	printf("fileSize %d with %d rows\n", fileSize, nrows);

	char *buffer = new char[fileSize];
	double *stackedQ = new double[nrows * 5];
	double *stackedU = new double[nrows * 3];
	double *stackedX = new double[nrows * 3];
	double *index = new double[nrows];

	//read
	fs.read(buffer, fileSize);
	printf("read %d chars\n", (int)fs.gcount());
	fs.close();

	//extract -- Hs in row-major
	double *pd = reinterpret_cast<double*>(buffer); //buffer is a byte array, but we know it to be double data
	//0i 1q0 2q1 3q2 4q3 5q4 : 6ux 7vx 8wx 9x 10uy 11vy 12wy 13y 14uz 15vz 16wz 17z 180 190 200 211
	for (int irow = 0; irow < nrows; irow ++){
		index[irow] = *(pd + 0 + irow*ncols);
		stackedQ[irow * 5  + 0] = *(pd + 1 + irow*ncols);
		stackedQ[irow * 5 + 1] = *(pd + 2 + irow*ncols);
		stackedQ[irow * 5 + 2] = *(pd + 3 + irow*ncols);
		stackedQ[irow * 5 + 3] = *(pd + 4 + irow*ncols);
		stackedQ[irow * 5 + 4] = *(pd + 5 + irow*ncols);

		stackedU[irow * 3 + 0] = *(pd + 6  + irow*ncols);
		stackedU[irow * 3 + 1] = *(pd + 10 + irow*ncols);
		stackedU[irow * 3 + 2] = *(pd + 14 + irow*ncols);

		stackedX[irow * 3 + 0] = *(pd + 9  + irow*ncols);
		stackedX[irow * 3 + 1] = *(pd + 13 + irow*ncols);
		stackedX[irow * 3 + 2] = *(pd + 17 + irow*ncols);

		//printf("%+5.4f: q[%+5.4f %+5.4f %+5.4f %+5.4f %+5.4f] ", index[irow], stackedQ[irow * 5 + 0], stackedQ[irow * 5 + 1], stackedQ[irow * 5 + 2], stackedQ[irow * 5 + 3], stackedQ[irow * 5 + 4]);
		//printf("u[%+5.4f %+5.4f %+5.4f]  ", stackedU[irow * 3 + 0], stackedU[irow * 3 + 1], stackedU[irow * 3 + 2]);
		//printf("x[%+5.4f %+5.4f %+5.4f]\n", stackedX[irow * 3 + 0], stackedX[irow * 3 + 1], stackedX[irow * 3 + 2]);
	}
	/*for (int i = 0; i < 800; i = i+10) {
		//printf("i[%d]: q[%+5.4f] u[%+5.4f] x[%+5.4f]\n", i, stackedQ[i], stackedU[i], stackedX[i]);
		//printf("i[%d]: q[%+5.4f] u[%+5.4f] x[%+5.4f]\n", i, stackedQ[(int)fmod(i, 5)], stackedU[i], stackedX[i]);
		printf("%+5.4f: q[%+5.4f %+5.4f %+5.4f %+5.4f %+5.4f] ", index[i], stackedQ[i * 5 + 0], stackedQ[i * 5 + 1], stackedQ[i * 5 + 2], stackedQ[i * 5 + 3], stackedQ[i * 5 + 4]);
		printf("u[%+5.4f %+5.4f %+5.4f]  ", stackedU[i * 3 + 0], stackedU[i * 3 + 1], stackedU[i * 3 + 2]);
		printf("x[%+5.4f %+5.4f %+5.4f]\n", stackedX[i * 3 + 0], stackedX[i * 3 + 1], stackedX[i * 3 + 2]);
	}//*/

	//test IP
	double k11up[11], k11dn[11], k11pm[11], k110[11], k50[5],k5pm[5],k5up[5],k5dn[5];
	//       tx01,            ty01,            tz01,        ry01,           rz01,          tx23,        ry34,          rz34,        cathL,        ry45,         rz45
	k110[0] = 806; k110[1] = -66.0; k110[2] = -28.0; k110[3] = 0; k110[4] = -.24; k110[5] = 8.2; k110[6] = 0;   k110[7] = 0; k110[8] = 98; k110[9] = 0; k110[10] = 0;
	k11pm[0] = 20; k11pm[1] = 20;   k11pm[2] = 20;   k11pm[3] = .2; k11pm[4] = .2; k11pm[5] = 1; k11pm[6] = .2; k11pm[7] = .2; k11pm[8] = 10; k11pm[9] = .2; k11pm[10] = .2;
	//
	k50[0] = 806; k50[1] = -66; k50[2] = -28; k50[3] = -.24; k50[4] = 95;
	k5pm[0] = 20; k5pm[1] = 20; k5pm[2] = 20; k5pm[3] = .2; k5pm[4] = 10;

	double nlArray[6]; //http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
		nlArray[0] = 1e5; // maxIts
		nlArray[1] = 6; // max time sec
		nlArray[3] = 1e-9; // min fun val
		nlArray[4] = 1e-9; // tol fun
		nlArray[5] = 1e-9; //tol x
		//nlArray[2] = 00; // GN_DIRECT
		//nlArray[2] = 01; // GN_DIRECT_L --locally biased
		//nlArray[2] = 03; // GN_DIRECT_L_RAND
		//nlArray[2] = 04; // GN_ESCH
		//nlArray[2] = 05; // GN_ISRES
		//nlArray[2] = 06; // GN_MLSL -- slow due to local searches
		//nlArray[2] = 07; // GN_MLSL_LDS
		//nlArray[2] = 12; // LN_BOBYQA
		nlArray[2] = 13; // LN_COBYLA
		//nlArray[2] = 14; // LN_NelderMead
		//nlArray[2] = 17; // LN_PRAXIS
		//nlArray[2] = 18; // LN_SUBPLX

	double q0Array[10];
		q0Array[0] = -.4; // joint minima
		q0Array[1] = -.2;
		q0Array[2] = -.2;	
		q0Array[3] = -.2;
		q0Array[4] = -10;
		q0Array[5] =  .4; // maxima
		q0Array[6] =  .2;
		q0Array[7] =  .2;
		q0Array[8] =  .2;
		q0Array[9] =  10;

	double fmin = 100;
	int ret = 0;
	double qps0[5];for (int i = 0; i < 5; i++) { qps0[i] = 0; }

	//ret = funIP_xyzdotu11A(nrows, stackedQ, stackedU, stackedX, qps0, k110, &fmin);	printf("val0 = %f\n", fmin);
	//for (int i = 0; i < 5; i++) { qps0[i] = 0; }
	//for (int i = 0; i < 11; i++) { k11up[i] = k110[i] + k11pm[i]; k11dn[i] = k110[i] - k11pm[i]; }
	//ret = estimatePmsQ_IPNLOpt_xyzdotu11A(nrows, stackedQ, stackedU, stackedX, k11up, k11dn, q0Array, nlArray, qps0, k110, &fmin);
	//printf("ret %d  fmin %f q0[", ret, fmin);
	//for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); }
	//printf("] k0[");
	//for (int i = 0; i < 11; i++) { printf("%8.3f ", k110[i]); }

	//ret = funIP_xyzpp11A(nrows, stackedQ, stackedU, stackedX, qps0, k110, &fmin); printf("val0 = %f\n", fmin);
	//for (int i = 0; i < 5; i++) { qps0[i] = 0; }
	//for (int i = 0; i < 11; i++) { k11up[i] = k110[i] + k11pm[i]; k11dn[i] = k110[i] - k11pm[i]; }
	//ret = estimatePmsQ_IPNLOpt_xyzpp11A(nrows, stackedQ, stackedU, stackedX, k11up, k11dn, q0Array, nlArray, qps0, k110, &fmin);
	//printf("ret %d  fmin %f\n", ret, fmin);


	ret = funIP_kn0_xyz5A(nrows, stackedQ, stackedX, k110, &fmin); printf("val0 = %f\n", fmin);
	for (int i = 0; i < 5; i++) { k5up[i] = k50[i] + k5pm[i]; k5dn[i] = k50[i] - k5pm[i]; }
	ret = estimate_kn0_xyz5A(nrows, stackedQ, stackedX, k50, k5up, k5dn, nlArray, &fmin);
	printf("ret %d  fmin %f\n", ret, fmin);
	printf("k50[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k50[i]); } printf("]\n");
	printf("k5u[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5up[i]); } printf("]\n");
	printf("k5d[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5dn[i]); } printf("]\n");

	std::cout << "\n\nPress Enter";
	std::getchar();
	printf("done\n");

	delete[] buffer; //everybody do your share
	delete[] stackedQ;
}


