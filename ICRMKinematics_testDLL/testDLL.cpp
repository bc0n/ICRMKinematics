#include <iostream>
#include <Eigen/Eigen>
#include "kinematicsDLL.h"

double qps0[5];    //joint angles at start of solver
double qpsGoal[5]; //joint angles associated with the xyzGoal (one possible solution to achieving xyzGoal)
double xyzGoal[3]; //the goal position
double kn05[5], kn06[6], kn11[11]; //kinematics params
double nlArray[6]; //nonlinear optimization params
double jArray[10]; //joint limits

void check_fk5A() {
	printf("kn05 = ["); for (int i = 0; i < 5; i++) { printf("%5.4f ", kn05[i]); } printf("]\n");

	qps0[0] = .5; qps0[1] = .3; qps0[2] = -.4; qps0[3] = 2; qps0[4] = 50;
	printf("qps0 = ["); for (int i = 0; i < 5; i++) { printf("%5.4f ", qps0[i]); } printf("]\n");

	double arrayH05[12];
	int ret = get5AH05(qps0, kn05, arrayH05); //arrayH05 is column-major

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
	ret = getTask5A_xyz(qpsGoal, kn05, xyzGoal);
	printf("NLOpt XYZ 5A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	ret = estimate_qps_xyz5A(qps, kn05, nlArray, jArray, xyz, &fmin);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] fmin=%8.3f ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], fmin, ret);
	std::cout << std::endl;

	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0; fmin = 1e3;
	ret = getTask6A_xyz(qpsGoal, kn06, xyzGoal);
	printf("NLOpt XYZ 6A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	ret = estimate_qps_xyz6A(qps, kn06, nlArray, jArray, xyz, &fmin);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] fmin=%8.3f ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], fmin, ret);
	std::cout << std::endl;


	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0; fmin = 1e3;
	ret = getTask11A_xyz(qpsGoal, kn11, xyzGoal);
	printf("NLOpt XYZ 11A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	ret = estimate_qps_xyz11A(qps, kn11, nlArray, jArray, xyz, &fmin);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] fmin=%8.3f ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], fmin, ret);
	std::cout << std::endl;
}

void check_nl_xyzuxuyuz() {
	double qpsGoal[5], qps[5], xyzGoal[3], xyz[3], uxyzGoal[3], uxyz[3], fmin;
	int ret = 99;

	qpsGoal[0] = .5; qpsGoal[1] = .5; qpsGoal[2] = -.3; qpsGoal[3] = 2; qpsGoal[4] = 30;

	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0; fmin = 1e3;
	ret = getTask5A_xyzuxuyuz(qpsGoal, kn05, xyzGoal, uxyzGoal);
	printf("NLOpt XYZUxUyUz 5A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] \n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2], uxyzGoal[0], uxyzGoal[1], uxyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	uxyz[0] = uxyzGoal[0]; uxyz[1] = uxyzGoal[1]; uxyz[2] = uxyzGoal[2];
	ret = estimate_qps_xyzuxuyuz5A(qps, kn06, nlArray, jArray, xyz, uxyz, &fmin);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] fmin=%8.3f ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], uxyz[0], uxyz[1], uxyz[2], fmin, ret);
	std::cout << std::endl;

	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0; fmin = 1e3;
	ret = getTask6A_xyzuxuyuz(qpsGoal, kn06, xyzGoal, uxyzGoal);
	printf("NLOpt XYZUxUyUz 6A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] \n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2], uxyzGoal[0], uxyzGoal[1], uxyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	uxyz[0] = uxyzGoal[0]; uxyz[1] = uxyzGoal[1]; uxyz[2] = uxyzGoal[2];
	ret = estimate_qps_xyzuxuyuz6A(qps, kn06, nlArray, jArray, xyz, uxyz, &fmin);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] fmin=%8.3f ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], uxyz[0], uxyz[1], uxyz[2], fmin, ret);
	std::cout << std::endl;

	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0; fmin = 1e3;
	ret = getTask11A_xyzuxuyuz(qpsGoal, kn11, xyzGoal, uxyzGoal);
	printf("NLOpt XYZUxUyUz 11A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] \n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2], uxyzGoal[0], uxyzGoal[1], uxyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	uxyz[0] = uxyzGoal[0]; uxyz[1] = uxyzGoal[1]; uxyz[2] = uxyzGoal[2];
	ret = estimate_qps_xyzuxuyuz11A(qps, kn11, nlArray, jArray, xyz, uxyz, &fmin);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] fmin=%8.3f ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], uxyz[0], uxyz[1], uxyz[2], fmin, ret);
	std::cout << std::endl;
}


int main(){
	//       tx01             ty01             tz01         ry01            rz01         ry34         rz34       kAlpha      eAlpha          lCath            ry45
	kn11[0] = 806; kn11[1] = -66.0; kn11[2] = -28.0; kn11[3] = 0; kn11[4] = -.24; kn11[5] = 0; kn11[6] = 0; kn11[7] = 1; kn11[8] = 1; kn11[9] = 95; kn11[10] = 0;
	//           tx01               ty01               tz01               rz01               ry34              lCath
	kn06[0] = kn11[0]; kn06[1] = kn11[1]; kn06[2] = kn11[2]; kn06[3] = kn11[4]; kn06[4] = kn11[5]; kn06[5] = kn11[9];
	//           tx01               ty01               tz01               rz01              lCath
	kn05[0] = kn11[0]; kn05[1] = kn11[1]; kn05[2] = kn11[2]; kn05[3] = kn11[4]; kn05[4] = kn11[9];


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
	check_nl_xyzuxuyuz();
	
	
	std::cout << "\n\nPress Enter";
	std::getchar();
	return 0;
}
