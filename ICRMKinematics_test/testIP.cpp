#include <fstream>
#include <iostream>
#include "kinematicsDLL.h"


char *buffer;
double *stackedQ;
double *stackedU;
double *stackedX;
double *index;
double nlArray[6], lnlArray[6];
double k11u[11], k11d[11], k110[11], k11[11];
double k60[6], k6u[6], k6d[6], k6[6];
double k50[5], k5u[5], k5d[5], k5[5];
double qpup[5], qpdn[5], qps0[5], qpupdn[10];
int ret, nrows;


void init() {
	//       tx01             ty01             tz01           ry01            rz01           ry34           rz34         kAlpha         eAlpha          lCath            ry45
	k11u[0] = 826; k11u[1] = -46.0; k11u[2] =  -8.0; k11u[3] =  .2; k11u[4] = -.04; k11u[5] =  .2; k11u[6] =  .2; k11u[7] = 1.2; k11u[8] = 1.2; k11u[9] = 110; k11u[10] =  .2;
	k110[0] = 806; k110[1] = -66.0; k110[2] = -28.0; k110[3] =   0; k110[4] = -.24; k110[5] =   0; k110[6] =   0; k110[7] =   1; k110[8] =   1; k110[9] =  95; k110[10] =   0;
	k11d[0] = 786; k11d[1] = -86.0; k11d[2] = -48.0; k11d[3] = -.2; k11d[4] = -.44; k11d[5] = -.2; k11d[6] = -.2; k11d[7] =  .8; k11d[8] =  .8; k11d[9] =  90; k11d[10] = -.2;	

	//          tx01              ty01              tz01              ry01              ry34             lCath
	k6u[0] = k11u[0]; k6u[1] = k11u[1]; k6u[2] = k11u[2]; k6u[3] = k11u[3]; k6u[4] = k11u[5]; k6u[5] = k11u[9];
	k60[0] = k110[0]; k60[1] = k110[1]; k60[2] = k110[2]; k60[3] = k110[3]; k60[4] = k110[5]; k60[5] = k110[9];
	k6d[0] = k11d[0]; k6d[1] = k11d[1]; k6d[2] = k11d[2]; k6d[3] = k11d[3]; k6d[4] = k11d[5]; k6d[5] = k11d[9];

	//          tx01              ty01              tz01              ry01             lCath
	k5u[0] = k11u[0]; k5u[1] = k11u[1]; k5u[2] = k11u[2]; k5u[3] = k11u[3]; k5u[4] = k11u[9];
	k50[0] = k110[0]; k50[1] = k110[1]; k50[2] = k110[2]; k50[3] = k110[3]; k50[4] = k110[9];
	k5d[0] = k11d[0]; k5d[1] = k11d[1]; k5d[2] = k11d[2]; k5d[3] = k11d[3]; k5d[4] = k11d[9];

	qpup[0] = .5; qpup[1] = .3; qpup[2] = .3; qpup[3] = 1; qpup[4] = 10;
	qpdn[0] = -.5; qpdn[1] = -.3; qpdn[2] = -.3; qpdn[3] = 1e-3; qpdn[4] = -10;
	qpupdn[0] = -.5; // joint minima
	qpupdn[1] = -.3;
	qpupdn[2] = -.3;
	qpupdn[3] = 1e-3;
	qpupdn[4] = -10;
	qpupdn[5] = .5; // maxima
	qpupdn[6] = .3;
	qpupdn[7] = .3;
	qpupdn[8] =  1;
	qpupdn[9] = 10;

	//http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
	nlArray[0] = 1e5; // maxIts
	nlArray[1] = 60; // max time sec
	nlArray[3] = 1e-9; // min fun val
	nlArray[4] = 1e-9; // tol fun
	nlArray[5] = 1e-9; //tol x
	//nlArray[2] = 00; // GN_DIRECT
	//nlArray[2] = 01; // GN_DIRECT_L --locally biased
	//nlArray[2] = 03; // GN_DIRECT_L_RAND
	//nlArray[2] = 06; // GN_ESCH
	//nlArray[2] = 07; // GN_ISRES
	//nlArray[2] = 08; // GN_MLSL
	nlArray[2] = 9; // GN_MLSL_LDS -- slow due to local searches
	//nlArray[2] = 12; // LN_BOBYQA
	//nlArray[2] = 13; // LN_COBYLA
	//nlArray[2] = 14; // LN_NelderMead
	//nlArray[2] = 17; // LN_PRAXIS
	//nlArray[2] = 18; // LN_SUBPLX

	lnlArray[0] = 1e9; // maxIts
	lnlArray[1] = 6; // max time sec
	lnlArray[3] = 1e-9; // min fun val
	lnlArray[4] = 5; // tol fun
	lnlArray[5] = .1; //tol x
	//lnlArray[2] = 12; // LN_BOBYQA
	//lnlArray[2] = 13; // LN_COBYLA
	//lnlArray[2] = 14; // LN_NelderMead
	//lnlArray[2] = 17; // LN_PRAXIS
	lnlArray[2] = 18; // LN_SUBPLX
}

//load data
void load38Col(char *fname) {
	int ncols = 38; //it, q0,q1,q2,q3,q4, H1, H2

	//open
	printf("\n\nOpening %s\n", fname);
	std::fstream fs;
	fs.open(fname, std::fstream::in | std::fstream::binary); //name, read/write mode|binary
	if (!fs) { std::cerr << "warning, could not open " << fname << std::endl; return; };

	//fileSize & initialize variables
	fs.seekg(0, fs.end);
	long fileSize = (long)fs.tellg();
	nrows = fileSize / sizeof(double) / ncols;
	fs.seekg(0, fs.beg);
	printf("fileSize %d with %d rows\n", fileSize, nrows);

	buffer = new char[fileSize];
	stackedQ = new double[nrows * 5];
	stackedU = new double[nrows * 3];
	stackedX = new double[nrows * 3];
	index = new double[nrows];

	//read
	fs.read(buffer, fileSize);
	printf("read %d chars\n\n", (int)fs.gcount());
	fs.close();

	//extract -- Hs in row-major
	double *pd = reinterpret_cast<double*>(buffer); //buffer is a byte array, but we know it to be double data
													//0i 1q0 2q1 3q2 4q3 5q4 : 6ux 7vx 8wx 9x 10uy 11vy 12wy 13y 14uz 15vz 16wz 17z 180 190 200 211
	for (int irow = 0; irow < nrows; irow++) {
		index[irow] = *(pd + 0 + irow*ncols);
		stackedQ[irow * 5 + 0] = *(pd + 1 + irow*ncols);
		stackedQ[irow * 5 + 1] = *(pd + 2 + irow*ncols);
		stackedQ[irow * 5 + 2] = *(pd + 3 + irow*ncols);
		stackedQ[irow * 5 + 3] = *(pd + 4 + irow*ncols);
		stackedQ[irow * 5 + 4] = *(pd + 5 + irow*ncols);

		stackedU[irow * 3 + 0] = *(pd + 6 + irow*ncols);
		stackedU[irow * 3 + 1] = *(pd + 10 + irow*ncols);
		stackedU[irow * 3 + 2] = *(pd + 14 + irow*ncols);

		stackedX[irow * 3 + 0] = *(pd + 9 + irow*ncols);
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
}
void loadColumnDat(char *fname) {
	//analogue of readDatXml.m version b34190bdaf4ac4bc3d9bd6a6c7349f01b7593d61, not parsing the xml
	printf("\n\nOpening %s\n", fname);
	std::fstream fs;
	fs.open(fname, std::fstream::in | std::fstream::binary);
	if (!fs) { std::cerr << "warning, could not open " << fname << std::endl; return; };

	//fileSize & initialize variables
	fs.seekg(0, fs.end);
	long fileSize = (long)fs.tellg();
	fs.seekg(0, fs.beg);
	printf("fileSize %d\n", fileSize);
	buffer = new char[fileSize];
	
	//read
	fs.read(buffer, fileSize);
	printf("read %d chars\n\n", (int)fs.gcount());
	fs.close();

	//extract -- Hs in row-major
	double *pd = reinterpret_cast<double*>(buffer); //buffer is a byte array, but we know it to be double data
	int na = (int)pd[0]; //number of points in the square
	int nb = (int)pd[1]; //number of points after IK and ramp fills
	nrows = nb; //as typically used
	//std::cout << na << " " << nb << std::endl;

	//commanded
	int ahsq = 2; //for (int i = Hsqa; i < Hsqa + 10; i++) { printf("%8.3f ", pd[i]); }
	int bhsq = ahsq + na * 16 -1;//for (int i = Hsqb-10; i <= Hsqb; i++) { printf("%8.3f ", pd[i]); }

	//converged
	int amap = bhsq + 1;
	int bmap = amap + na;
	//for (int i = amap; i < amap + 10; i++) { printf("%8.3f ", pd[i]); } printf("\n");
	//for (int i = bmap-10; i <= bmap; i++) { printf("%8.3f ", pd[i]); }printf("\n");
	int aqcv = bmap + 1;
	int bqcv = aqcv + 5 * na - 1;
	//for (int i = aqcv; i < aqcv + 10; i++) { printf("%8.3f ", pd[i]); } printf("\n");
	//for (int i = bqcv-10; i <= bqcv; i++) { printf("%8.3f ", pd[i]); }printf("\n");
	int axcm = bqcv + 1;
	int bxcm = axcm + 3 * na - 1;
	//for (int i = axcm; i < axcm + 10; i++) { printf("%8.3f ", pd[i]); } printf("\n");
	//for (int i = bxcm-10; i <= bxcm; i++) { printf("%8.3f ", pd[i]); }printf("\n");
	int axcv = bxcm + 1;
	int bxcv = axcv + 3 * na - 1;
	//for (int i = axcv; i < axcv + 10; i++) { printf("%8.3f ", pd[i]); } printf("\n");
	//for (int i = bxcv-10; i <= bxcv; i++) { printf("%8.3f ", pd[i]); }printf("\n");
	int anrm = bxcv + 1;
	int bnrm = anrm + na - 1;
	//for (int i = anrm; i < anrm + 10; i++) { printf("%8.3f ", pd[i]); } printf("\n");
	//for (int i = bnrm-10; i <= bnrm; i++) { printf("%8.3f ", pd[i]); }printf("\n");
	int aret = bnrm + 1;
	int bret = aret + na - 1;
	//for (int i = aret; i < aret + 10; i++) { printf("%8.3f ", pd[i]); } printf("\n");
	//for (int i = bret-10; i <= bret; i++) { printf("%8.3f ", pd[i]); }printf("\n");

	//measured
	stackedQ = new double[nb * 5];
	stackedU = new double[nb * 3];
	stackedX = new double[nb * 3];
	int aqms = bret + 1;
	int bqms = aqms + 5 * nb -1;
	for (int i = 0; i < 5*nb; i++) {
		stackedQ[i] = pd[i+aqms];
	}
	//for (int i = 5*9; i < 5*13; i++) { printf("%8.3f ", stackedQ[i]); } printf("\n");
	//for (int i = (5*nb)-100; i < 5*nb; i++) { printf("%8.3f\n", stackedQ[i]); } printf("\n");
	int ahms = bqms + 1;
	int bhms = ahms + 16 * nb - 1;
	//Hms is saved row-major
	//for (int i = ahms; i < ahms + 16; i++) { printf("%8.5f ", pd[i]); } printf("\n");
	for (int ib = 0; ib < nb; ib++) {
		for (int irow = 0; irow < 3; irow++) {
			//printf("%d %d %8.5f\n", ib, irow, pd[ahms + ib * 16 + irow*4]);
			stackedU[ib * 3 + irow] = pd[ahms + ib * 16 + irow * 4 + 0];
			stackedX[ib * 3 + irow] = pd[ahms + ib * 16 + irow * 4 + 3];
		}
	}
	//for (int i = 0; i < 10; i++) { printf("%8.5f ", stackedU[i]); } printf("\n");
	//for (int i = 0; i < 10; i++) { printf("%8.5f ", stackedX[i]); } printf("\n");
	//for (int i = 3*nb-10; i < 3*nb; i++) { printf("%8.5f ", stackedU[i]); } printf("\n"); //d.Hms(end-3:end,1:3,1)'
	//for (int i = 3*nb-10; i < 3*nb; i++) { printf("%8.5f ", stackedX[i]); } printf("\n");
}
void loadNQXU(char *fname) {
	printf("\n\nOpening %s\n", fname);
	std::fstream fs;
	fs.open(fname, std::fstream::in | std::fstream::binary);
	if (!fs) { std::cerr << "warning, could not open " << fname << std::endl; return; };

	//fileSize & initialize variables
	fs.seekg(0, fs.end);
	long fileSize = (long)fs.tellg();
	fs.seekg(0, fs.beg);
	printf("fileSize %d\n", fileSize);
	buffer = new char[fileSize];

	//read
	fs.read(buffer, fileSize);
	printf("read %d chars\n\n", (int)fs.gcount());
	fs.close();

	//extract-- format is fwrite( [n, reshape(qps,1,n*5), reshape(Hms(:,1:3,4),1,n*3), reshape(Hms(:,1:3,1),1,n*3) ] )
	double *pd = reinterpret_cast<double*>(buffer); //buffer is a byte array, but we know it to be double data
	nrows = (int)pd[0]; //number of points in each of qps, xyz, uxyz	

	//for (int i = 0; i < 100; i++) { std::cout << *(pd + i + nrows*5) << "\n"; } std::cout << std::endl;

	stackedQ = new double[nrows * 5];
	stackedU = new double[nrows * 3];
	stackedX = new double[nrows * 3];

	for (int i = 0; i < nrows * 5; i++) {
		stackedQ[i] = -7.7;
	}
	for (int i = 0; i < nrows * 3; i++) {
		stackedX[i] = -8.8;
	}
	for (int i = 0; i < nrows * 3; i++) {
		stackedU[i] = -9.9;
	}

	for (int i = 0; i < nrows; i++) {
		stackedQ[i * 5 + 0] = *(pd + 1+0 + i*5);
		stackedQ[i * 5 + 1] = *(pd + 1+1 + i*5);
		stackedQ[i * 5 + 2] = *(pd + 1+2 + i*5);
		stackedQ[i * 5 + 3] = *(pd + 1+3 + i*5);
		stackedQ[i * 5 + 4] = *(pd + 1+4 + i*5);

		stackedX[i * 3 + 0] = *(pd + 1+nrows*5+0 + i*3);
		stackedX[i * 3 + 1] = *(pd + 1+nrows*5+1 + i*3);
		stackedX[i * 3 + 2] = *(pd + 1+nrows*5+2 + i*3);

		stackedU[i * 3 + 0] = *(pd + 1+nrows*5+nrows*3+0 + i*3);
		stackedU[i * 3 + 1] = *(pd + 1+nrows*5+nrows*3+1 + i*3);
		stackedU[i * 3 + 2] = *(pd + 1+nrows*5+nrows*3+2 + i*3);
	}
	//for (int i = 0; i < 100; i++) {
	//	printf("i[ %d ] q[%8.4f %8.4f %8.4f %8.4f %8.4f] x[%8.4f %8.4f %8.4f] u[%8.4f %8.4f %8.4f]\n", i, stackedQ[i * 5 + 0], stackedQ[i * 5 + 1], stackedQ[i * 5 + 2], stackedQ[i * 5 + 3], stackedQ[i * 5 + 4], stackedX[i * 3 + 0], stackedX[i * 3 + 1], stackedX[i * 3 + 2], stackedU[i * 3 + 0], stackedU[i * 3 + 1], stackedU[i * 3 + 2]);
	//}
	//for (int i = nrows-100; i < nrows; i++) {
	//	printf("i[ %d ] q[%8.4f %8.4f %8.4f %8.4f %8.4f] x[%8.4f %8.4f %8.4f] u[%8.4f %8.4f %8.4f]\n", i, stackedQ[i * 5 + 0], stackedQ[i * 5 + 1], stackedQ[i * 5 + 2], stackedQ[i * 5 + 3], stackedQ[i * 5 + 4], stackedX[i * 3 + 0], stackedX[i * 3 + 1], stackedX[i * 3 + 2], stackedU[i * 3 + 0], stackedU[i * 3 + 1], stackedU[i * 3 + 2]);
	//}
}

// find qp0
void check_qp0_xyzuxuyuz5a() {
	double fmin = 22;
	printf("\nqp0_xyzuxuyuz5a\n");
	for (int i = 0; i < 5; i++) { qps0[i] = .2; }
	
	ret = fun_qp0_xyzuxuyuz5A(nrows, stackedQ, stackedX, stackedU, qps0, k50, &fmin); printf("funQp = %f\n", fmin);
	ret = estimate_qp0_xyzuxuyuz5A(nrows, stackedQ, stackedX, stackedU, qps0, qpupdn, k50, nlArray, &fmin);

	printf("ret %d  fmin %f\n", ret, fmin);
	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpup[i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpdn[i]); } printf("]\n");
}
void check_qp0_xyzuxuyuz6a() {
	double fmin = 22;
	printf("\nqp0_xyzuxuyuz6a\n");
	for (int i = 0; i < 5; i++) { qps0[i] = .2; }

	ret = fun_qp0_xyzuxuyuz6A(nrows, stackedQ, stackedX, stackedU, qps0, k60, &fmin); printf("funQp = %f\n", fmin);
	ret = estimate_qp0_xyzuxuyuz6A(nrows, stackedQ, stackedX, stackedU, qps0, qpupdn, k60, nlArray, &fmin);

	printf("ret %d  fmin %f\n", ret, fmin);
	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpup[i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpdn[i]); } printf("]\n");
}
void check_qp0_xyzuxuyuz11a() {
	double fmin = 22;
	printf("\nqp0_xyzuxuyuz11a\n");
	for (int i = 0; i < 5; i++) { qps0[i] = .2; }

	ret = fun_qp0_xyzuxuyuz11A(nrows, stackedQ, stackedX, stackedU, qps0, k110, &fmin); printf("funQp = %f\n", fmin);
	ret = estimate_qp0_xyzuxuyuz11A(nrows, stackedQ, stackedX, stackedU, qps0, qpupdn, k110, nlArray, &fmin);

	printf("ret %d  fmin %f\n", ret, fmin);
	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpup[i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpdn[i]); } printf("]\n");
}
void check_qp0_xyzuxuyuz11a_mlsl() {
	double fmin = 22;
	printf("\nqp0_xyzuxuyuz11a_mlsl\n");
	for (int i = 0; i < 5; i++) { qps0[i] = .2; }

	ret = fun_qp0_xyzuxuyuz11A(nrows, stackedQ, stackedX, stackedU, qps0, k110, &fmin); printf("funQp = %f\n", fmin);
	ret = estimate_qp0_xyzuxuyuz11A_mlsllds(nrows, stackedQ, stackedX, stackedU, qps0, qpupdn, k110, nlArray, lnlArray, &fmin);

	printf("ret %d  fmin %f\n", ret, fmin);
	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpup[i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpdn[i]); } printf("]\n");
}

// find kn0
void check_kn0_xyzuxuyuz5a() {
	double fmin = 22;
	printf("\nkn0_xyzuxuyuz5a\n");
	for (int i = 0; i < 5; i++) { k5[i] = k50[i]; }

	ret = fun_kn0_xyzuxuyuz5A(nrows, stackedQ, stackedX, stackedU, k5, &fmin); printf("val0 = %f\n", fmin);
	ret = estimate_kn0_xyzuxuyuz5A(nrows, stackedQ, stackedX, stackedU, k5, k5u, k5d, nlArray, &fmin);
	
	printf("ret %d  fmin %f\n", ret, fmin);
	printf("k5u[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5u[i]); } printf("]\n");
	printf("k50[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5[i]); } printf("]\n");
	printf("k5d[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5d[i]); } printf("]\n");
}
void check_kn0_xyzuxuyuz6a() {
	double fmin = 22;
	printf("\nkn0_xyzuxuyuz6a\n");
	for (int i = 0; i < 6; i++) { k6[i] = k60[i]; }

	ret = fun_kn0_xyzuxuyuz6A(nrows, stackedQ, stackedX, stackedU, k6, &fmin); printf("val0 = %f\n", fmin);

	ret = estimate_kn0_xyzuxuyuz6A(nrows, stackedQ, stackedX, stackedU, k6, k6u, k6d, nlArray, &fmin);

	printf("ret %d  fmin %f\n", ret, fmin);
	printf("k6u[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6u[i]); } printf("]\n");
	printf("k6 [");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6[i]); } printf("]\n");
	printf("k6d[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6d[i]); } printf("]\n");
}
void check_kn0_xyzuxuyuz11a() {
	double fmin = 22;
	printf("\nkn0_xyzuxuyuz11a\n");
	for (int i = 0; i < 11; i++) { k11[i] = k110[i]; }

	ret = fun_kn0_xyzuxuyuz11A(nrows, stackedQ, stackedX, stackedU, k11, &fmin); printf("val0 = %f\n", fmin);

	ret = estimate_kn0_xyzuxuyuz11A(nrows, stackedQ, stackedX, stackedU, k11, k11u, k11d, nlArray, &fmin);

	printf("ret %d  fmin %f\n", ret, fmin);
	printf("k11u[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11u[i]); } printf("]\n");
	printf("k110[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11[i]); } printf("]\n");
	printf("k11d[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11d[i]); } printf("]\n");
}

void check_kn0_xyzuxuyuz5a_sub() {
	double fmin = 220000;
	printf("\nkn0_xyzuxuyuz5a_sub\n");
	for (int i = 0; i < 5; i++) { k5[i] = k50[i]; }
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }

	double knSubset[] = { 1, 1, 1, 1, 1 };
	//double knSubset[] = { 0, 0, 0, 1, 1 };
	//double knSubset[] = { 0, 0, 0, 1, 1}; //true = estimate

	ret = fun_kn0_xyzuxuyuz5A(nrows, stackedQ, stackedX, stackedU, k5, &fmin); printf("fun_kn0 = %f\n", fmin);
	ret = estimate_kn0_xyzuxuyuz5A_subset(nrows, stackedQ, stackedX, stackedU, k5, k5u, k5d, knSubset, nlArray, &fmin);

	printf("k5u[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5u[i]); } printf("]\n");
	printf("k50[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k50[i]); } printf("]\n");
	printf("k5 [");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5[i]); } printf("]\n");
	printf("k5d[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5d[i]); } printf("]\n");
	printf("ret %d  fmin %f\n", ret, fmin);
}
void check_kn0_xyzuxuyuz6a_sub() {
	double fmin = 220000;
	printf("\nkn0_xyzuxuyuz6a_sub\n");
	for (int i = 0; i < 6; i++) { k6[i] = k60[i]; }

	double knSubset[] = { 1, 1, 1, 1, 1, 1 };
	//double knSubset[] = { 0, 0, 0, 1, 1, 0 };
	//double knSubset[] = { 0, 0, 0, 1, 1, 0 }; //true = estimate

	ret = fun_kn0_xyzuxuyuz6A(nrows, stackedQ, stackedX, stackedU, k6, &fmin); printf("fun_kn0 = %f\n", fmin);
	ret = estimate_kn0_xyzuxuyuz6A_subset(nrows, stackedQ, stackedX, stackedU, k6, k6u, k6d, knSubset, nlArray, &fmin);

	printf("k6u[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6u[i]); } printf("]\n");
	printf("k60[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k60[i]); } printf("]\n");
	printf("k6 [");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6[i]); } printf("]\n");
	printf("k6d[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6d[i]); } printf("]\n");
	printf("ret %d  fmin %f\n", ret, fmin);
}
void check_kn0_xyzuxuyuz11a_sub() {
	double fmin = 220000;
	printf("\nkn0_xyzuxuyuz11a_sub\n");
	for (int i = 0; i < 11; i++) { k11[i] = k110[i]; }

	//double knSubset[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	//double knSubset[] = { 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0 };
	double knSubset[] = { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0 }; //true = estimate

	ret = fun_kn0_xyzuxuyuz11A(nrows, stackedQ, stackedX, stackedU, k11, &fmin); printf("fun_kn0 = %f\n", fmin);
	ret = estimate_kn0_xyzuxuyuz11A_subset(nrows, stackedQ, stackedX, stackedU, k11, k11u, k11d, knSubset, nlArray, &fmin);

	printf("k11u[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11u[i]); } printf("]\n");
	printf("k110[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k110[i]); } printf("]\n");
	printf("k11 [");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11[i]); } printf("]\n");
	printf("k11d[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11d[i]); } printf("]\n");
	printf("ret %d  fmin %f\n", ret, fmin);
}
void check_kn0_xyzuxuyuz11a_sub_mlsl() {
	double fmin = 220000;
	printf("\nkn0_xyzuxuyuz11a_sub_mlsl\n");
	for (int i = 0; i < 11; i++) { k11[i] = k110[i]; }

	//double knSubset[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	//double knSubset[] = { 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0 };
	double knSubset[] = { 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0 }; //true = estimate

	ret = fun_kn0_xyzuxuyuz11A(nrows, stackedQ, stackedX, stackedU, k11, &fmin); printf("fun_kn0 = %f\n", fmin);
	ret = estimate_kn0_xyzuxuyuz11A_subset_mlsllds(nrows, stackedQ, stackedX, stackedU, k11, k11u, k11d, knSubset, nlArray, lnlArray, &fmin);

	printf("k11u[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11u[i]); } printf("]\n");
	printf("k110[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k110[i]); } printf("]\n");
	printf("k11 [");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11[i]); } printf("]\n");
	printf("k11d[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11d[i]); } printf("]\n");
	printf("ret %d  fmin %f\n", ret, fmin);
}

//simultaneous qp0 & kn0
void check_qp0kn0_xyz5A_sub() {
	double fmin = 220000;
	printf("\nqp0kn0_xyz5a_sub\n");
	for (int i = 0; i < 5; i++) { k5[i] = k50[i]; }
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }
	double knSubset[] = { 1, 1, 1, 1, 1 }; //true = estimate
	ret = fun_qp0kn0_xyz5A(nrows, stackedQ, stackedX, qps0, k5, &fmin); printf("fun_qp0kn0 = %f\n", fmin);
	ret = estimate_qp0kn0_xyz5A_subset(nrows, stackedQ, stackedX, qps0, qpupdn, k5, k5u, k5d, knSubset, nlArray, &fmin);

	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[5 + i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[i]); } printf("]\n");
	printf("k5u[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5u[i]); } printf("]\n");
	printf("k50[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k50[i]); } printf("]\n");
	printf("k5 [");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5[i]); } printf("]\n");
	printf("k5d[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5d[i]); } printf("]\n");
	printf("ret %d  fmin %f\n", ret, fmin);
}
void check_qp0kn0_xyz6A_sub() {
	double fmin = 220000;
	printf("\nqp0kn0_xyz6a_sub\n");
	for (int i = 0; i < 6; i++) { k6[i] = k60[i]; }
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }
	double knSubset[] = { 1, 1, 1, 1, 1, 1}; //true = estimate
	ret = fun_qp0kn0_xyz6A(nrows, stackedQ, stackedX, qps0, k6, &fmin); printf("fun_qp0kn0 = %f\n", fmin);
	ret = estimate_qp0kn0_xyz6A_subset(nrows, stackedQ, stackedX, qps0, qpupdn, k6, k6u, k6d, knSubset, nlArray, &fmin);

	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[5 + i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[i]); } printf("]\n");
	printf("k11u[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6u[i]); } printf("]\n");
	printf("k110[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k60[i]); } printf("]\n");
	printf("k11 [");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6[i]); } printf("]\n");
	printf("k11d[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6d[i]); } printf("]\n");
	printf("ret %d  fmin %f\n", ret, fmin);
}
void check_qp0kn0_xyz11A_sub() {
	double fmin = 220000;
	printf("\nqp0kn0_xyz11a_sub\n");
	for (int i = 0; i < 11; i++) { k11[i] = k110[i]; }
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }
	double knSubset[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	//double knSubset[] = {0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0};
	//double knSubset[] = { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0 }; //true = estimate
	ret = fun_qp0kn0_xyz11A(nrows, stackedQ, stackedX, qps0, k11, &fmin); printf("fun_qp0kn0 = %f\n", fmin);
	ret = estimate_qp0kn0_xyz11A_subset(nrows, stackedQ, stackedX, qps0, qpupdn, k11, k11u, k11d, knSubset, nlArray, &fmin);

	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[5 + i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[i]); } printf("]\n");
	printf("k11u[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11u[i]); } printf("]\n");
	printf("k110[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k110[i]); } printf("]\n");
	printf("k11 [");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11[i]); } printf("]\n");
	printf("k11d[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11d[i]); } printf("]\n");
	printf("ret %d  fmin %f\n", ret, fmin);
}
void check_qp0kn0_xyz11A_sub_mlsl() {
	double fmin = 220000;
	printf("\nqp0kn0_xyz11a_sub\n");
	for (int i = 0; i < 11; i++) { k11[i] = k110[i]; }
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }
	double knSubset[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	//double knSubset[] = {0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0};
	//double knSubset[] = { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0 }; //true = estimate

	ret = fun_qp0kn0_xyz11A(nrows, stackedQ, stackedX, qps0, k11, &fmin); printf("fun_qp0kn0 = %f\n", fmin);
	ret = estimate_qp0kn0_xyz11A_subset_mlsllds(nrows, stackedQ, stackedX, qps0, qpupdn, k11, k11u, k11d, knSubset, nlArray, lnlArray, &fmin);

	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[5 + i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[i]); } printf("]\n");
	printf("k11u[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11u[i]); } printf("]\n");
	printf("k110[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k110[i]); } printf("]\n");
	printf("k11 [");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11[i]); } printf("]\n");
	printf("k11d[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11d[i]); } printf("]\n");
	printf("ret %d  fmin %f\n", ret, fmin);
}

void check_qp0kn0_xyzuxuyuz5A() {
	double fmin = 220000;
	printf("\nqp0kn0_xyzuxuyuz5a\n");
	for (int i = 0; i < 5; i++) { k5[i] = k50[i]; }
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }

	ret = fun_qp0kn0_xyzuxuyuz5A     (nrows, stackedQ, stackedX, stackedU, qps0,         k5,                    &fmin); printf("fun_qp0kn0 = %f\n", fmin );
	ret = estimate_qp0kn0_xyzuxuyuz5A(nrows, stackedQ, stackedX, stackedU, qps0, qpupdn, k5, k5u, k5d, nlArray, &fmin);

	printf("ret %d  fmin %f\n", ret, fmin);
	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[5 + i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[i]); } printf("]\n");
	printf("k5u[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5u[i]); } printf("]\n");
	printf("k50[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5[i]); } printf("]\n");
	printf("k5d[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5d[i]); } printf("]\n");
}
void check_qp0kn0_xyzuxuyuz6A() {
	double fmin = 220000;
	printf("\nqp0kn0_xyzuxuyuz6a\n");
	for (int i = 0; i < 6; i++) { k6[i] = k60[i]; }
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }
	
	ret = fun_qp0kn0_xyzuxuyuz6A(     nrows, stackedQ, stackedX, stackedU, qps0,         k6,                    &fmin); printf("fun_qp0kn0 = %f\n", fmin);
	ret = estimate_qp0kn0_xyzuxuyuz6A(nrows, stackedQ, stackedX, stackedU, qps0, qpupdn, k6, k6u, k6d, nlArray, &fmin);

	printf("ret %d  fmin %f\n", ret, fmin);
	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[5 + i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[i]); } printf("]\n");
	printf("k6u[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6u[i]); } printf("]\n");
	printf("k6 [");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6[i]); } printf("]\n");
	printf("k6d[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6d[i]); } printf("]\n");
}
void check_qp0kn0_xyzuxuyuz11A() {
	double fmin = 220000;
	printf("\nqp0kn0_xyzuxuyuz11a\n");
	for (int i = 0; i < 11; i++) { k11[i] = k110[i]; }
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }

	ret = fun_qp0kn0_xyzuxuyuz11A(     nrows, stackedQ, stackedX, stackedU, qps0,         k11,                      &fmin); printf("fun_qp0kn0 = %f\n", fmin);
	ret = estimate_qp0kn0_xyzuxuyuz11A(nrows, stackedQ, stackedX, stackedU, qps0, qpupdn, k11, k11u, k11d, nlArray, &fmin);

	printf("ret %d  fmin %f\n", ret, fmin);
	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[5 + i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[i]); } printf("]\n");
	printf("k11u[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11u[i]); } printf("]\n");
	printf("k11 [");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11[i]); } printf("]\n");
	printf("k11d[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11d[i]); } printf("]\n");
}

void check_qp0kn0_xyzuxuyuz5A_sub() {
	double fmin = 220000;
	printf("\nqp0kn0_xyzuxuyuz5a_sub\n");
	for (int i = 0; i < 5; i++) { k5[i] = k50[i]; }
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }

	//double knSubset[] = { 1, 1, 1, 1, 1 };
	double knSubset[] = { 0, 0, 0, 1, 1};
	//double knSubset[] = { 1, 1, 1, 1, 1}; //true = estimate

	ret = fun_qp0kn0_xyzuxuyuz5A(nrows, stackedQ, stackedX, stackedU, qps0, k5, &fmin); printf("fun_qp0kn0 = %f\n", fmin);
	ret = estimate_qp0kn0_xyzuxuyuz5A_subset(nrows, stackedQ, stackedX, stackedU, qps0, qpupdn, k5, k5u, k5d, knSubset, nlArray, &fmin);

	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[5 + i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[i]); } printf("]\n");
	printf("k5u[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5u[i]); } printf("]\n");
	printf("k50[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k50[i]); } printf("]\n");
	printf("k5 [");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5[i]); } printf("]\n");
	printf("k5d[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", k5d[i]); } printf("]\n");
	printf("ret %d  fmin %f\n", ret, fmin);
}
void check_qp0kn0_xyzuxuyuz6A_sub() {
	double fmin = 220000;
	printf("\nqp0kn0_xyzuxuyuz6a_sub\n");
	for (int i = 0; i < 6; i++) { k6[i] = k60[i]; }
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }

	//double knSubset[] = { 1, 1, 1, 1, 1, 1};
	double knSubset[] = { 0, 0, 0, 1, 1, 1};
	//double knSubset[] = { 1, 1, 0, 0, 1, 1}; //true = estimate

	ret = fun_qp0kn0_xyzuxuyuz6A(nrows, stackedQ, stackedX, stackedU, qps0, k6, &fmin); printf("fun_qp0kn0 = %f\n", fmin);
	ret = estimate_qp0kn0_xyzuxuyuz6A_subset(nrows, stackedQ, stackedX, stackedU, qps0, qpupdn, k6, k6u, k6d, knSubset, nlArray, &fmin);

	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[5 + i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[i]); } printf("]\n");
	printf("k6u[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6u[i]); } printf("]\n");
	printf("k60[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k60[i]); } printf("]\n");
	printf("k6 [");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6[i]); } printf("]\n");
	printf("k6d[");  for (int i = 0; i < 6; i++) { printf("%8.3f ", k6d[i]); } printf("]\n");
	printf("ret %d  fmin %f\n", ret, fmin);
}
void check_qp0kn0_xyzuxuyuz11A_sub() {
	double fmin = 220000;
	printf("\nqp0kn0_xyzuxuyuz11a_sub\n");
	for (int i = 0; i < 11; i++) { k11[i] = k110[i]; }
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }

	double knSubset[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	//double knSubset[] = {0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0};
	//double knSubset[] = { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0 }; //true = estimate

	ret = fun_qp0kn0_xyzuxuyuz11A(nrows, stackedQ, stackedX, stackedU, qps0, k11, &fmin); printf("fun_qp0kn0 = %f\n", fmin);
	ret = estimate_qp0kn0_xyzuxuyuz11A_subset(nrows, stackedQ, stackedX, stackedU, qps0, qpupdn, k11, k11u, k11d, knSubset, nlArray, &fmin);
	
	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[5 + i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[i]); } printf("]\n");
	printf("k11u[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11u[i]); } printf("]\n");
	printf("k110[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k110[i]); } printf("]\n");
	printf("k11 [");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11[i]); } printf("]\n");
	printf("k11d[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11d[i]); } printf("]\n");
	printf("ret %d  fmin %f\n", ret, fmin);
}
void check_qp0kn0_xyzuxuyuz11A_sub_mlsl() {
	double fmin = 220000;
	printf("\nqp0kn0_xyzuxuyuz11a_sub_mlsl\n");
	for (int i = 0; i < 11; i++) { k11[i] = k110[i]; }
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }

	double knSubset[] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	//double knSubset[] = {0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0};
	//double knSubset[] = { 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0 }; //true = estimate

	ret = fun_qp0kn0_xyzuxuyuz11A(nrows, stackedQ, stackedX, stackedU, qps0, k11, &fmin); printf("fun_qp0kn0 = %f\n", fmin);
	ret = estimate_qp0kn0_xyzuxuyuz11A_subset_mlsllds(nrows, stackedQ, stackedX, stackedU, qps0, qpupdn, k11, k11u, k11d, knSubset, nlArray, lnlArray, &fmin);

	printf("qpup[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[5 + i]); } printf("]\n");
	printf("qps0[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qps0[i]); } printf("]\n");
	printf("qpdn[");  for (int i = 0; i < 5; i++) { printf("%8.3f ", qpupdn[i]); } printf("]\n");
	printf("k11u[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11u[i]); } printf("]\n");
	printf("k110[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k110[i]); } printf("]\n");
	printf("k11 [");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11[i]); } printf("]\n");
	printf("k11d[");  for (int i = 0; i < 11; i++) { printf("%8.3f ", k11d[i]); } printf("]\n");
	printf("ret %d  fmin %f\n", ret, fmin);
}

void main() {
	init();
	//load38Col("testSquareQps.dat");
	//loadColumnDat("testSquareXYZ_i0_n2942_160622_191046.dat");
	loadNQXU("forTestIPCpp.dat");
	
	//check_qp0_xyzuxuyuz5a();
	//check_qp0_xyzuxuyuz6a();
	//check_qp0_xyzuxuyuz11a();
	check_qp0_xyzuxuyuz11a_mlsl();

	//check_kn0_xyzuxuyuz5a();
	//check_kn0_xyzuxuyuz5a_sub();
	//check_kn0_xyzuxuyuz6a();
	//check_kn0_xyzuxuyuz6a_sub();
	//check_kn0_xyzuxuyuz11a();
	//check_kn0_xyzuxuyuz11a_sub();
	check_kn0_xyzuxuyuz11a_sub_mlsl();

	//check_qp0kn0_xyz5A_sub();
	//check_qp0kn0_xyz6A_sub();
	//check_qp0kn0_xyz11A_sub();
	check_qp0kn0_xyz11A_sub_mlsl();
	//check_qp0kn0_xyzuxuyuz5A();
	//check_qp0kn0_xyzuxuyuz5A_sub();
	//check_qp0kn0_xyzuxuyuz6A();
	//check_qp0kn0_xyzuxuyuz6A_sub();
	//check_qp0kn0_xyzuxuyuz11A();
	//check_qp0kn0_xyzuxuyuz11A_sub();
	check_qp0kn0_xyzuxuyuz11A_sub_mlsl();

	//estimate joint angles by iterative IK(Hms) & IP(qps^, Hms)?
	
	std::cout << "\n\nPress Enter";
	std::getchar();
	printf("done\n");

	delete[] buffer; //everybody do your share
	delete[] stackedQ;
}


