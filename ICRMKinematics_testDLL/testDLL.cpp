#include <iostream>
#include "kinematicsDLL.h"

double qps0[5];    //joint angles at start of solver
double qpsGoal[5]; //joint angles associated with the xyzGoal (one possible solution to achieving xyzGoal)
double xyzGoal[3]; //the goal position
double kn5A[5], kn6A[6], kn11A[11]; //kinematics params
double nlArray[6]; //nonlinear optimization params
double jArray[10]; //joint limits

void check_fk5A() {
	printf("kn5A = ["); for (int i = 0; i < 5; i++) { printf("%5.4f ", kn5A[i]); } printf("]\n");

	qps0[0] = .5; qps0[1] = .3; qps0[2] = -.4; qps0[3] = 2; qps0[4] = 50;
	printf("qps0 = ["); for (int i = 0; i < 5; i++) { printf("%5.4f ", qps0[i]); } printf("]\n");

	double arrayH05[12];
	int ret = get5AH05(qps0, kn5A, arrayH05); //arrayH05 is column-major

	Eigen::Matrix4d H;
	H.setIdentity();
	for (int i = 0; i < 4; i++) {//column
		for (int j = 0; j < 3; j++) {//row
			H(j, i) = arrayH05[i * 3 + j];
		}
	}
	printf("ret %d\n", ret);
	std::cout << H << std::endl;
}

void check_nl_xyz(){
	double qpsGoal[5], qps[5], xyzGoal[3], xyz[3], fmin;
	int ret = 99;
	
	qpsGoal[0] = .5; qpsGoal[1] = .5; qpsGoal[2] = -.3; qpsGoal[3] = 2; qpsGoal[4] = 30;
	
	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0; fmin = 1e3;
	ret = getTask5A_xyz(qpsGoal, kn5A, xyzGoal);
	printf("NLOpt XYZ 5A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	ret = getQps_IKnlopt_xyz5A(qps, kn5A, nlArray, jArray, xyz, &fmin);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] fmin=%8.3f ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], fmin, ret);
	std::cout << std::endl;

	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0; fmin = 1e3;
	ret = getTask6A_xyz(qpsGoal, kn6A, xyzGoal);
	printf("NLOpt XYZ 6A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	ret = getQps_IKnlopt_xyz6A(qps, kn6A, nlArray, jArray, xyz, &fmin);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] fmin=%8.3f ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], fmin, ret);
	std::cout << std::endl;


	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0; fmin = 1e3;
	ret = getTask11A_xyz(qpsGoal, kn11A, xyzGoal);
	printf("NLOpt XYZ 11A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	ret = getQps_IKnlopt_xyz11A(qps, kn11A, nlArray, jArray, xyz, &fmin);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] fmin=%8.3f ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], fmin, ret);
	std::cout << std::endl;
}

void check_nl_xyzuxuyuz() {
	double qpsGoal[5], qps[5], xyzGoal[3], xyz[3], uxyzGoal[3], uxyz[3], fmin;
	int ret = 99;

	qpsGoal[0] = .5; qpsGoal[1] = .5; qpsGoal[2] = -.3; qpsGoal[3] = 2; qpsGoal[4] = 30;

	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0; fmin = 1e3;
	ret = getTask6A_xyzuxuyuz(qpsGoal, kn6A, xyzGoal, uxyzGoal);
	printf("NLOpt XYZUxUyUz 6A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] \n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2], uxyzGoal[0], uxyzGoal[1], uxyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	uxyz[0] = uxyzGoal[0]; uxyz[1] = uxyzGoal[1]; uxyz[2] = uxyzGoal[2];
	ret = getQps_IKnlopt_xyzuxuyuz6A(qps, kn6A, nlArray, jArray, xyz, uxyz, &fmin);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] fmin=%8.3f ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], uxyz[0], uxyz[1], uxyz[2], fmin, ret);
	std::cout << std::endl;


	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0; fmin = 1e3;
	ret = getTask11A_xyzuxuyuz(qpsGoal, kn11A, xyzGoal, uxyzGoal);
	printf("NLOpt XYZUxUyUz 11A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] \n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2], uxyzGoal[0], uxyzGoal[1], uxyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	uxyz[0] = uxyzGoal[0]; uxyz[1] = uxyzGoal[1]; uxyz[2] = uxyzGoal[2];
	ret = getQps_IKnlopt_xyzuxuyuz11A(qps, kn11A, nlArray, jArray, xyz, uxyz, &fmin);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] fmin=%8.3f ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], uxyz[0], uxyz[1], uxyz[2], fmin, ret);
	std::cout << std::endl;
}


int main(){
	kn11A[0] = 806; kn11A[1] = -66.0; kn11A[2] = -28.0; kn11A[3] = 0; kn11A[4] = -.24; //TxyzRyz01
	kn11A[5] = 8.2; //Tx23
	kn11A[6] = 0; kn11A[7] = 0; kn11A[8] = 98; //Ryz34lCath
	kn11A[9] = 0; 
	kn11A[10] = 0; //Ryz45

	kn6A[0] = kn11A[8]; //[mm] lCath
	kn6A[1] = kn11A[4]; //Rz01
	kn6A[2] = kn11A[0]; //Tx01
	kn6A[3] = kn11A[1]; //Ty01
	kn6A[4] = kn11A[2]; //Tz01
	kn6A[5] = kn11A[5]; //Tx23

	kn5A[0] = 806; //Tx01
	kn5A[1] = -66; //Ty01
	kn5A[2] = -28; //Tz01
	kn5A[3] = -.24; //Rz01
	kn5A[4] = 95; //[mm] lCath

	nlArray[0] = 1e9; // maxIts
	nlArray[1] = 60; // max time sec
	nlArray[3] = 1e-6; // min fun val
	nlArray[4] = 1e-6; // tol fun
	nlArray[5] = 1e-6; //tol x
	//nlArray[2] = 00; // GN_DIRECT
	//nlArray[2] = 01; // GN_DIRECT_L --locally biased
	//nlArray[2] = 03; // GN_DIRECT_L_RAND
	//nlArray[2] = 04; // GN_ESCH
	//nlArray[2] = 05; // GN_ISRES
	//nlArray[2] = 06; // GN_MLSL -- slow due to local searches
	//nlArray[2] = 07; // GN_MLSL_LDS
	//nlArray[2] = 12; // LN_BOBYQA
	//nlArray[2] = 13; // LN_COBYLA
	nlArray[2] = 14; // LN_NelderMead
	//nlArray[2] = 17; // LN_PRAXIS
	//nlArray[2] = 18; // LN_SUBPLX
	//nlArray[2] = 19; // GD_MLSL = 19,
	//nlArray[2] = 20; // GD_MLSL_LDS = 20,
	//nlArray[2] = 21; // GD_STOGO = 21,
	//nlArray[2] = 22; // GD_STOGO_RAND = 22,
	//nlArray[2] = 23; // LD_CCSAQ = 23,
	//nlArray[2] = 24; // LD_LBFGS = 24,
	//nlArray[2] = 25; // LD_LBFGS_NOCEDAL = 25,
	//nlArray[2] = 26; // LD_MMA = 26,
	//nlArray[2] = 27; // LD_TNEWTON = 27,
	//nlArray[2] = 28; // LD_TNEWTON_RESTART = 28,
	//nlArray[2] = 29; // LD_TNEWTON_PRECOND = 29,
	//nlArray[2] = 30; // LD_TNEWTON_PRECOND_RESTART = 30,
	//nlArray[2] = 31; // LD_VAR1 = 31,
	//nlArray[2] = 32; // LD_VAR2 = 32,
	//nlArray[2] = 33; // LD_SLSQP = 33

	
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
	
	//check_fk5A();
	check_nl_xyz();
	//check_nl_xyzuxuyuz();
	
	
	std::cout << "\n\nPress Enter";
	std::getchar();
	return 0;
}
