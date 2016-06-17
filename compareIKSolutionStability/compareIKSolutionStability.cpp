#include <fstream>
#include <iostream>
#include "kinematicsDLL.h"

double kn6A[6], kn11A[11], nrArray[5], nlArray[6], jArray[10];
long nrows;
double *stackedQ, *stackedX, *stackedU;
double *stackedQResult, *stackedXResult, *stackedUResult;

void readQpsFile(char *fname){
	long ncols = 38, fileSize;
	double xyzG[3] = { 0,0,0 }, ppG[2] = { 0,0 };

	// read qps|H|H file: 38 doubles arranged in [0i 1q0 2q1 3q2 4q3 5q4 : 6ux 7vx 8wx 9x 10uy 11vy 12wy 13y 14uz 15vz 16wz 17z 180 190 200 211 : 22ux2...
	printf("\n\nOpening %s\n", fname);
	std::fstream fs;
	fs.open(fname, std::fstream::in | std::fstream::binary); //name, read/write mode|binary
	if (!fs) { std::cerr << "warning, could not open " << fname << std::endl; return; };

	//fileSize & initialize variables
	fs.seekg(0, fs.end);
	fileSize = (long)fs.tellg();
	nrows = fileSize / sizeof(double) / ncols;
	fs.seekg(0, fs.beg);
	//nrows = 1;
	printf("fileSize %d with %d rows\n", fileSize, nrows);

	char *buffer = new char[fileSize];
	stackedQ = new double[nrows * 5];
	stackedU = new double[nrows * 3];
	stackedX = new double[nrows * 3];
	//double *index = new double[nrows];

	//read
	fs.read(buffer, fileSize);
	printf("read %d chars\n", (int)fs.gcount());
	fs.close();

	//extract -- Hs in row-major
	double *pd = reinterpret_cast<double*>(buffer); //buffer is a byte array, but we know it to be double data
	//0i 1q0 2q1 3q2 4q3 5q4 : 6ux 7vx 8wx 9x 10uy 11vy 12wy 13y 14uz 15vz 16wz 17z 180 190 200 211
	for (int irow = 0; irow < nrows; irow++) {
		//index[irow] = *(pd + 0 + irow*ncols);
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
	//for (int i = 0; i < 800; i = i + 10) {
	//	//printf("i[%d]: q[%+5.4f] u[%+5.4f] x[%+5.4f]\n", i, stackedQ[i], stackedU[i], stackedX[i]);
	//	//printf("i[%d]: q[%+5.4f] u[%+5.4f] x[%+5.4f]\n", i, stackedQ[(int)fmod(i, 5)], stackedU[i], stackedX[i]);
	//	printf("%+5.4f: q[%+5.4f %+5.4f %+5.4f %+5.4f %+5.4f] ", index[i], stackedQ[i * 5 + 0], stackedQ[i * 5 + 1], stackedQ[i * 5 + 2], stackedQ[i * 5 + 3], stackedQ[i * 5 + 4]);
	//	printf("u[%+5.4f %+5.4f %+5.4f]  ", stackedU[i * 3 + 0], stackedU[i * 3 + 1], stackedU[i * 3 + 2]);
	//	printf("x[%+5.4f %+5.4f %+5.4f]\n", stackedX[i * 3 + 0], stackedX[i * 3 + 1], stackedX[i * 3 + 2]);
	//}
}

void writeResultsToFile(char *taskFKname, int im){

	//open file
	char fname[100];
	sprintf_s(fname, "%s_n%d_m%d_%fs.dat", taskFKname, nrows, im, nlArray[1]);
	printf("Opening %s\n", fname);
	std::fstream fs;
	fs.open(fname, std::fstream::out | std::fstream::binary | std::fstream::trunc); //name, write mode|binary
	if (!fs) { std::cerr << "warning, could not open " << fname << std::endl; return; };
	
	//write stackedQ and stackedX to file
	fs.write(reinterpret_cast<char *>(stackedQResult), nrows * 5 * 8);
	fs.write(reinterpret_cast<char *>(stackedXResult), nrows * 3 * 8);
	fs.write(reinterpret_cast<char *>(stackedUResult), nrows * 3 * 8);
	
	fs.close();
	printf("Closed\n");
}


int main() {

	kn6A[0] = 98; //[mm] cathL
	kn6A[1] = -.24; //Rz01
	kn6A[2] = 806; //Tx01
	kn6A[3] = -66; //Ty01
	kn6A[4] = -28; //Tz01
	kn6A[5] = 8.2; //Tx23

	kn11A[0] = 806; kn11A[1] = -66.0; kn11A[2] = -28.0; kn11A[3] = 0; kn11A[4] = -.24; //TxyzRyz01
	kn11A[5] = 8.2; //Tx23
	kn11A[6] = 0; kn11A[7] = 0; kn11A[8] = 98; //Ryz34CathL
	kn11A[9] = 0; kn11A[10] = 0; //Ryz45

	jArray[0] = -2; // joint minima
	jArray[1] = -.8;
	jArray[2] = -1;
	jArray[3] = 1e-3;
	jArray[4] = 0;
	jArray[5] = 2; // maxima
	jArray[6] = .8;
	jArray[7] = 1;
	jArray[8] = 5;
	jArray[9] = 500;

	nlArray[0] = 1e9; // maxIts
	nlArray[1] = 1; // max time sec
	nlArray[3] = 1e-9; // min fun val
	nlArray[4] = 1e-11; // tol fun
	nlArray[5] = 1e-9; //tol x
	nlArray[2] = 00; //GN_DIRECT
	//nlArray[2] = 01; // GN_DIRECT_L --locally biased
	//nlArray[2] = 03; // GN_DIRECT_L_RAND
	//nlArray[2] = 04; // GN_ESCH
	//nlArray[2] = 05; // GN_ISRES
	//nlArray[2] = 06; // GN_MLSL -- slow due to local searches
	//nlArray[2] = 07; // GN_MLSL_LDS
	//nlArray[2] = 12; // LN_BOBYQA
	//nlArray[2] = 13; // LN_COBYLA
	//nlArray[2] = 14; // LN_NelderMead
	//nlArray[2] = 17; // LN_PRAXIS
	//nlArray[2] = 18; // LN_SUBPLX


	char fname[] = "testSquareQps.dat";
	double a = -1;
	stackedQ = &a; stackedX = &a; stackedU = &a;

	//read saved data
	readQpsFile(fname);

	//run IKs on data, varying IKmethod and Task function
	int imethods[11] = { 0,1,2,3,8,12,13,14,18 }; //IKmethods to run
	int nmethods = 1;

	//setup after reading file
	double qps[5], xyz[3], uxyz[3];
	int i = 0;
	stackedQResult = new double[nrows * 5];
	stackedXResult = new double[nrows * 3];
	stackedUResult = new double[nrows * 3];
	int *retResult = new int[nrows * 3];

	//xyz
	//*
	for (int im = 0; im < nmethods; im++) {
		nlArray[2] = imethods[im];
		printf("Starting xyz6A method %d\n", imethods[im]);

		for (long irow = 0; irow < nrows; irow++) {
			if (fmod((double)irow, 100)) { printf("%d of %d\n", irow, nrows); }
			//extract target
			for (i = 0; i < 3; i++) { 
				xyz[i] = stackedX[irow * 3 + i];
			}

			//run ik; subsequent search begins at previous qps
			retResult[irow] = getQps_IKnlopt_xyz6A(qps, kn6A, nlArray, jArray, xyz);

			//save for writing
			for (i = 0; i < 3; i++) {
				stackedXResult[irow * 3 + i] = xyz[i];
			}
			for (i = 0; i < 5; i++) {
				stackedQResult[irow * 5 + i] = qps[i];
			}
			//printf("qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", stackedQResult[irow * 5 + 0], stackedQResult[irow * 5 + 1], stackedQResult[irow * 5 + 2], stackedQResult[irow * 5 + 3], stackedQResult[irow * 5 + 4]);
			
		}//irow
		writeResultsToFile("xyz6A", imethods[im]);
	}//im

	for (int im = 0; im < nmethods; im++) {
		nlArray[2] = imethods[im];
		printf("Starting xyz11A method %d\n", imethods[im]);

		for (long irow = 0; irow < nrows; irow++) {
			if (fmod((double)irow, 100)) { printf("%d of %d\n", irow, nrows); }
			//extract target
			for (i = 0; i < 3; i++) {
				xyz[i] = stackedX[irow * 3 + i];
			}

			//run ik; subsequent search begins at previous qps
			retResult[irow] = getQps_IKnlopt_xyz6A(qps, kn6A, nlArray, jArray, xyz);

			//save for writing
			for (i = 0; i < 3; i++) {
				stackedXResult[irow * 3 + i] = xyz[i];
			}
			for (i = 0; i < 5; i++) {
				stackedQResult[irow * 5 + i] = qps[i];
			}
			//printf("qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", stackedQResult[irow * 5 + 0], stackedQResult[irow * 5 + 1], stackedQResult[irow * 5 + 2], stackedQResult[irow * 5 + 3], stackedQResult[irow * 5 + 4]);

		}//irow
		writeResultsToFile("xyz11A", imethods[im]);
	}//im
	//*/

	//xyzuxuyuz
	//*
	for (int im = 0; im < nmethods; im++) {
		nlArray[2] = imethods[im];
		printf("Starting xyzuxuyuz6A method %d\n", im);

		for (long irow = 0; irow < nrows; irow++) {
			if (fmod((double)irow, 100)) { printf("%d of %d\n", irow, nrows); }
			//extract target
			for (i = 0; i < 3; i++) {
				xyz[i] = stackedX[irow * 3 + i];
				uxyz[i] = stackedU[irow * 3 + i];
			}

			//run ik; subsequent search begins at previous qps
			retResult[irow] = getQps_IKnlopt_xyzuxuyuz6A(qps, kn6A, nlArray, jArray, xyz, uxyz);

			//save for writing
			for (i = 0; i < 3; i++) {
				stackedXResult[irow * 3 + i] = xyz[i];
				stackedUResult[irow * 3 + i] = uxyz[i];
			}
			for (i = 0; i < 5; i++) {
				stackedQResult[irow * 5 + i] = qps[i];
			}
			//printf("qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", stackedQResult[irow * 5 + 0], stackedQResult[irow * 5 + 1], stackedQResult[irow * 5 + 2], stackedQResult[irow * 5 + 3], stackedQResult[irow * 5 + 4]);
		}//irow
		writeResultsToFile("xyzuxuyuz6A", imethods[im]);
	}//im

	for (int im = 0; im < nmethods; im++) {
		nlArray[2] = imethods[im];
		printf("Starting xyzuxuyuz11A method %d\n", im);

		for (long irow = 0; irow < nrows; irow++) {
			if (fmod((double)irow, 100)) { printf("%d of %d\n", irow, nrows); }
			//extract target
			for (i = 0; i < 3; i++) {
				xyz[i] = stackedX[irow * 3 + i];
				uxyz[i] = stackedU[irow * 3 + i];
			}

			//run ik; subsequent search begins at previous qps
			retResult[irow] = getQps_IKnlopt_xyzuxuyuz11A(qps, kn6A, nlArray, jArray, xyz, uxyz);

			//save for writing
			for (i = 0; i < 3; i++) {
				stackedXResult[irow * 3 + i] = xyz[i];
				stackedUResult[irow * 3 + i] = uxyz[i];
			}
			for (i = 0; i < 5; i++) {
				stackedQResult[irow * 5 + i] = qps[i];
			}
			//printf("qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", stackedQResult[irow * 5 + 0], stackedQResult[irow * 5 + 1], stackedQResult[irow * 5 + 2], stackedQResult[irow * 5 + 3], stackedQResult[irow * 5 + 4]);
		}//irow
		writeResultsToFile("xyzuxuyuz11A", imethods[im]);
	}//im
	//*/


	delete[] stackedQ,stackedX,stackedU;
	delete[] stackedQResult, stackedXResult, stackedUResult;
	std::cout << "\n\nPress Enter";
	std::getchar();
	return 0;
}