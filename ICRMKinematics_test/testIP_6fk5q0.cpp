
#include "testIP_6fk5q0.h"

void main() {
	//load test file
	char fname[] = "D:\\20160401_inter2d_estimate\\06_parameterLoop\\testSquareXYZ_n6208_160621_160312.dat";

	//open
	printf("\n\nOpening %s\n", fname);
	std::fstream fs;
	fs.open(fname, std::fstream::in | std::fstream::binary); //name, read/write mode|binary
	if (!fs) { std::cerr << "warning, could not open " << fname << std::endl; return; };

	//fileSize & initialize variables
	fs.seekg(0, fs.end);
	long fileSize = (long)fs.tellg();
	fs.seekg(0, fs.beg);
	printf("fileSize %d\n", fileSize);
	
	//read
	char *buffer = new char[fileSize];
	fs.read(buffer, fileSize);
	printf("read %d chars\n", (int)fs.gcount());
	fs.close();

	//extract -- Hs in row-major
	double *pd = reinterpret_cast<double*>(buffer); //buffer is a byte array, but we know it to be double data
	int na = (int)pd[0];
	int nb = (int) pd[1];
	printf("na=%d, nb=%d\n", na, nb);

	double *stackedQ = new double[nb * 5];
	double *stackedU = new double[nb * 3];
	double *stackedX = new double[nb * 3];

	//data is ordered [na,nb,Hsquare, indexMap,qCnv,xyzCmd,xyzCnv,nrmCnv,retCnv, qMeasured,Hmeasured]
	//we want to invP on the measureds which are nb long
	int iqm = 1 + 1 + na*(16 + 1 + 5 + 3 + 3 + 1 + 1) + 1; //last +1 for 0indexing
	int ihm = iqm + nb * 5;
	//printf("iqm = %d, ihm = %d\n", iqm, ihm);
	//for (int i = iqm - 10; i < iqm + 30; i++) {
	//	printf("%d: %f\n", i, pd[i]);
	//}
	//for (int i = ihm - 10; i < ihm + 30; i++) {
	//	printf("%d: %f\n", i, pd[i]);
	//}
	for (int i = 0; i < nb; i++) {
		stackedQ[i * 5 + 0] = pd[iqm + i * 5 + 0];
		stackedQ[i * 5 + 1] = pd[iqm + i * 5 + 1];
		stackedQ[i * 5 + 2] = pd[iqm + i * 5 + 2];
		stackedQ[i * 5 + 3] = pd[iqm + i * 5 + 3];
		stackedQ[i * 5 + 4] = pd[iqm + i * 5 + 4];
		stackedU[i * 3 + 0] = pd[ihm + i * 16 + 0];
		stackedU[i * 3 + 1] = pd[ihm + i * 16 + 4];
		stackedU[i * 3 + 2] = pd[ihm + i * 16 + 8];
		stackedX[i * 3 + 0] = pd[ihm + i * 16 + 3];
		stackedX[i * 3 + 1] = pd[ihm + i * 16 + 7];
		stackedX[i * 3 + 2] = pd[ihm + i * 16 + 11];
	}
	//for (int i = nb-20; i < nb; i++) {
	//	printf("q[%5.4f %5.4f %5.4f %5.4f %5.4f] x[%5.4f %5.4f %5.4f] u[%5.4f %5.4f %5.4f]\n", stackedQ[i*5+0], stackedQ[i * 5 + 1], stackedQ[i * 5 + 2], stackedQ[i * 5 + 3], stackedQ[i * 5 + 4], stackedX[i * 3 + 0], stackedX[i * 3 + 1], stackedX[i * 3 + 2], stackedU[i * 3 + 0], stackedU[i * 3 + 1], stackedU[i * 3 + 2]);
	//}

	//test IP
	double k11up[11], k11dn[11], k11pm[11], k110[11];
	//       tx01,            ty01,            tz01,        ry01,           rz01,          tx23,        ry34,          rz34,        cathL,        ry45,         rz45
	k110[0] = 806; k110[1] = -66.0; k110[2] = -28.0; k110[3] = 0; k110[4] = -.24; k110[5] = 8.2; k110[6] = 0;   k110[7] = 0; k110[8] = 98; k110[9] = 0; k110[10] = 0;
	k11pm[0] = 20; k11pm[1] = 20;   k11pm[2] = 20;   k11pm[3] = .2; k11pm[4] = .2; k11pm[5] = 1; k11pm[6] = .2; k11pm[7] = .2; k11pm[8] = 10; k11pm[9] = .2; k11pm[10] = .2;
	double nlArray[6]; //http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
		nlArray[0] = 1e5; // maxIts
		nlArray[1] = 60; // max time sec
		nlArray[3] = 1e-9; // min fun val
		nlArray[4] = 1e-9; // tol fun
		nlArray[5] = 1e-9; //tol x
		nlArray[2] = 00; // GN_DIRECT
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
	double qps0[5];
	//for (int i = 0; i < 5; i++) { qps0[i] = 0; }
	qps0[0] = 0; qps0[1] = 0; qps0[2] = 0; qps0[3] = 1; qps0[4] = 100;

	ret = funIP_xyzdotu11A(nb, stackedQ, stackedU, stackedX, qps0, k110, &fmin);	printf("val0 = %f\n", fmin);
	ret = funIP_xyzpp11A(nb, stackedQ, stackedU, stackedX, qps0, k110, &fmin);	printf("val0 = %f\n", fmin);

	//ret = estimatePmsQ_IPNLOpt_xyzdotu11A_assumeX0(nb, stackedQ, stackedU, stackedX, k11up, k11dn, q0Array, nlArray, &fmin);
	//printf("ret %d  fmin %f\n", ret, fmin);
	for (int i = 0; i < 5; i++) { qps0[i] = 0; }
	for (int i = 0; i < 11; i++) { k11up[i] = k110[i] + k11pm[i]; k11dn[i] = k110[i] - k11pm[i]; }
	ret = estimatePmsQ_IPNLOpt_xyzdotu11A(nb, stackedQ, stackedU, stackedX, k11up, k11dn, q0Array, nlArray, qps0, k110, &fmin);
	printf("ret %d  fmin %f\n\n", ret, fmin);
	//for (int i = 0; i < 5; i++) { qps0[i] = 0; }
	//for (int i = 0; i < 11; i++) { k11up[i] = k110[i] + k11pm[i]; k11dn[i] = k110[i] - k11pm[i]; }
	//ret = estimatePmsQ_IPNLOpt_xyzpp11A(nb, stackedQ, stackedU, stackedX, k11up, k11dn, q0Array, nlArray, qps0, k110, &fmin);
	//printf("ret %d  fmin %f\n", ret, fmin);


	delete[] buffer; //everybody do your share
	delete[] stackedQ;
	//*/

	std::cout << "\n\nPress Enter";
	std::getchar();
	printf("done\n");

	
}


