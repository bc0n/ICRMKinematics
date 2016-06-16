#include <iostream>
#include "kinematicsDLL.h"

double qps0[5];    //joint angles at start of solver
double qpsGoal[5]; //joint angles associated with the xyzGoal (one possible solution to achieving xyzGoal)
double xyzGoal[3]; //the goal position
double kn6A[6]; //kinematics params
double kn11A[11];
double nrArray[5]; //newton-raphson params
double nlArray[6]; //nonlinear optimization params
double jArray[10]; //joint limits

void check_nl_xyz(){
	double qpsGoal[5], qps[5], xyzGoal[3], xyz[3];
	int ret = 99;
	
	qpsGoal[0] = .5; qpsGoal[1] = .5; qpsGoal[2] = -.3; qpsGoal[3] = 2; qpsGoal[4] = 30;
	
	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0;
	ret = getTask6A_xyz(qpsGoal, kn6A, xyzGoal);

	printf("NLOpt XYZ 6A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	ret = getQps_IKnlopt_xyz6A(qps, kn6A, nlArray, jArray, xyz);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], ret);
	std::cout << std::endl;


	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0;
	ret = getTask11A_xyz(qpsGoal, kn11A, xyzGoal);

	printf("NLOpt XYZ 11A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	ret = getQps_IKnlopt_xyz11A(qps, kn11A, nlArray, jArray, xyz);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], ret);
	std::cout << std::endl;
}

void check_nl_xyzuxuyuz() {
	double qpsGoal[5], qps[5], xyzGoal[3], xyz[3], uxyzGoal[3],uxyz[3];
	int ret = 99;

	qpsGoal[0] = .5; qpsGoal[1] = .5; qpsGoal[2] = -.3; qpsGoal[3] = 2; qpsGoal[4] = 30;

	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0;
	ret = getTask6A_xyzuxuyuz(qpsGoal, kn6A, xyzGoal, uxyzGoal);

	printf("NLOpt XYZUxUyUz 6A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] \n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2], uxyzGoal[0], uxyzGoal[1], uxyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	uxyz[0] = uxyzGoal[0]; uxyz[1] = uxyzGoal[1]; uxyz[2] = uxyzGoal[2];
	ret = getQps_IKnlopt_xyzuxuyuz6A(qps, kn6A, nlArray, jArray, xyz, uxyz);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], uxyz[0], uxyz[1], uxyz[2], ret);
	std::cout << std::endl;


	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0;
	ret = getTask11A_xyzuxuyuz(qpsGoal, kn11A, xyzGoal, uxyzGoal);

	printf("NLOpt XYZUxUyUz 11A\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] \n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2], uxyzGoal[0], uxyzGoal[1], uxyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	uxyz[0] = uxyzGoal[0]; uxyz[1] = uxyzGoal[1]; uxyz[2] = uxyzGoal[2];
	ret = getQps_IKnlopt_xyzuxuyuz11A(qps, kn11A, nlArray, jArray, xyz, uxyz);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] uxyz[%8.3f %8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], uxyz[0], uxyz[1], uxyz[2], ret);
	std::cout << std::endl;
}

int main(){
	
	kn6A[0] = 103; //[mm] cathL
	kn6A[1] = -.3; //Rz01
	kn6A[2] = 814; //Tx01
	kn6A[3] = -159; //Ty01
	kn6A[4] = -8.3; //Tz01
	kn6A[5] = 8.2; //Tx23

	kn11A[0] = 806; kn11A[1] = -66.0; kn11A[2] = -28.0; kn11A[3] = 0; kn11A[4] = -.24; //TxyzRyz01
	kn11A[5] = 8.2; //Tx23
	kn11A[6] = 0; kn11A[7] = 0; kn11A[8] = 98; //Ryz34CathL
	kn11A[9] = 0; kn11A[10] = 0; //Ryz45

	nrArray[0] = 1e-3; //epsilon
	nrArray[1] = 1e-3; //errTol
	nrArray[2] = 100; //maxIts
	nrArray[3] = .3; //primary factor
	nrArray[4] = .01; //secondary factor;
	
	nlArray[0] = 1e9; // maxIts
	nlArray[1] = 60; // max time sec
	nlArray[3] = 1e-9; // min fun val
	nlArray[4] = 1e-11; // tol fun
	nlArray[5] = 1e-9; //tol x
	nlArray[2] = 00; // GN_DIRECT
	//nlArray[2] = 8; // MLSL w/ NM
	//nlArray[2] = 14; //LN_NelderMead
	
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
	
	check_nl_xyz();
	check_nl_xyzuxuyuz();
	
	
	std::cout << "\n\nPress Enter";
	std::getchar();
	return 0;
}
