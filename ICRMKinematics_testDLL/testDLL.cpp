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
////////////////////////////////////////////////////////////fontSizeChecker///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void checkFwdK6A(double *qps, double *knArray) {
	int ret = 0;
	double arrayH01[12], arrayH02[12], arrayH03[12], arrayH04[12], arrayH05[12];
	//H01 = obj.Tx(766.6)*obj.Ty(-112)*obj.Tz(14.7)*obj.Rz(-27*pi/180) * obj.Rx(qp(1)); #ok<*PROP> %prox roll
	ret = get6AH01(qps, knArray, arrayH01);
	std::cout << "H01 = " << std::endl;
	for(int i = 0; i < 12; i++) { printf("%3.4f ", arrayH01[i]); }
	std::cout << std::endl << std::endl;

	//H12 = obj.Rz(qp(2)); %pitch
	ret = get6AH02(qps, knArray, arrayH02);
	std::cout << "H02 = " << std::endl;
	for (int i = 0; i < 12; i++) { printf("%3.4f ", arrayH02[i]); }
	std::cout << std::endl << std::endl;

	ret = get6AH03(qps, knArray, arrayH03);
	std::cout << "H03 = " << std::endl;
	for (int i = 0; i < 12; i++) { printf("%3.4f ", arrayH03[i]); }
	std::cout << std::endl << std::endl;

	ret = get6AH04(qps, knArray, arrayH04);
	std::cout << "H04 = " << std::endl;
	for (int i = 0; i < 12; i++) { printf("%3.4f ", arrayH04[i]); }
	std::cout << std::endl << std::endl;

	ret = get6AH05(qps, knArray, arrayH05);
	std::cout << "H05 = " << std::endl;
	for (int i = 0; i < 12; i++) { printf("%3.4f ", arrayH05[i]); }
	std::cout << std::endl << std::endl;
}
void checkFwdK11A(double *qps, double *knArray) {
	int ret = 0;
	double arrayH01[12], arrayH02[12], arrayH03[12], arrayH04[12], arrayH05[12];
	//ret = get11AH01(qps, knArray, arrayH01);
	//std::cout << "H01 = " << std::endl;
	//for (int i = 0; i < 12; i++) { printf("%3.4f ", arrayH01[i]); }
	//std::cout << std::endl << std::endl;

	////H12 = obj.Rz(qp(2)); %pitch
	//ret = get11AH02(qps, knArray, arrayH02);
	//std::cout << "H02 = " << std::endl;
	//for (int i = 0; i < 12; i++) { printf("%3.4f ", arrayH02[i]); }
	//std::cout << std::endl << std::endl;

	//ret = get11AH03(qps, knArray, arrayH03);
	//std::cout << "H03 = " << std::endl;
	//for (int i = 0; i < 12; i++) { printf("%3.4f ", arrayH03[i]); }
	//std::cout << std::endl << std::endl;

	//ret = get11AH04(qps, knArray, arrayH04);
	//std::cout << "H04 = " << std::endl;
	//for (int i = 0; i < 12; i++) { printf("%3.4f ", arrayH04[i]); }
	//std::cout << std::endl << std::endl;

	ret = get11AH05(qps, knArray, arrayH05);
	std::cout << "H05 = " << std::endl;
	for (int i = 0; i < 12; i++) { printf("%3.4f ", arrayH05[i]); }
	std::cout << std::endl << std::endl;
}


void check_nl_xyz(){
	double qpsGoal[5], qps[5], xyzGoal[3], xyz[3];
	int ret = 99;

	qps[0] = 0; qps[1] = 0; qps[2] = 0; qps[3] = 0; qps[4] = 0;
	qpsGoal[0] = .5; qpsGoal[1] = .5; qpsGoal[2] = -.3; qpsGoal[3] = 2; qpsGoal[4] = 30;

	ret = get6A_xyz(qpsGoal, kn6A, xyzGoal);

	printf("NLOpt Point\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	ret = getQps_IKnlopt_xyz6A(qps, kn6A, nlArray, jArray, xyz);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], ret);
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
	
	nlArray[0] = 1e6; // maxIts
	nlArray[1] = 60; // max time sec
	nlArray[3] = 1e-9; // min fun val
	nlArray[4] = 1e-11; // tol fun
	nlArray[5] = 1e-9; //tol x

	//nlArray[2] = 00; // GN_DIRECT
	nlArray[2] = 8; // MLSL w/ NM
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
	
	
	std::cout << "\n\nPress Enter";
	std::getchar();
	return 0;
}
