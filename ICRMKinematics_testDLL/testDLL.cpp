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

void check_nrSinglePoint(double *qpsGoal, double *xyzGoal, double *qps, double *knArray, double *nrArray, double *jArray, double *xyz) {
	int ret = 99;
	printf("Newton Point\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	ret = getQps_IKNewtonSinglePoint(qps, knArray, nrArray, jArray, xyz);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], ret);
	std::cout << std::endl;
}
void check_nrPriPointPhiPsi(double *qpsGoal, double *xyzGoal, double *ppGoal, double *qps, double *knArray, double *nrArray, double *jArray, double *xyz, double *pp) {
	int ret = 99;
	printf("Newton Priority Point>PhiPsi\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] = pp[%8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2], ppGoal[0], ppGoal[1]);
	ret = getQps_IKNewtonPriorityPointPhiPsi(qps, knArray, nrArray, jArray, xyz, pp);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] = pp[%8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], pp[0], pp[1], ret);
	std::cout << std::endl;
}
void check_nrAugPointPhiPsi(double *qpsGoal, double *xyzGoal, double *ppGoal, double *qps, double *knArray, double *nrArray, double *jArray, double *xyz, double *pp) {
	int ret = 99;
	printf("Newton Augmented PointPhiPsi\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] = pp[%8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2], ppGoal[0], ppGoal[1]);
	ret = getQps_IKNewtonAugmentedPointPhiPsiStops(qps, knArray, nrArray, jArray, xyz, pp);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] = pp[%8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], pp[0], pp[1], ret);
	std::cout << std::endl;
}
void check_nrPriPointCtrSum(double *qpsGoal, double *xyzGoal, double *qps, double *knArray, double *nrArray, double *jArray, double *xyz){
	int ret = 99;
	printf("Newton Priority PointCenterSum\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	ret = getQps_IKNewtonPriorityPointCenterSum(qps, knArray, nrArray, jArray, xyz);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], ret);
	std::cout << std::endl;
}
void check_nrPriPointCtrSumStops(double *qpsGoal, double *xyzGoal, double *qps, double *knArray, double *nrArray, double *jArray, double *xyz) {
	int ret = 99;
	printf("Newton Priority PointCenterSumStops\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	ret = getQps_IKNewtonPriorityPointCenterSumStops(qps, knArray, nrArray, jArray, xyz);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], ret);
	std::cout << std::endl;
}
void check_nrAugPointCtrSumStops(double *qpsGoal, double *xyzGoal, double *qps, double *knArray, double *nrArray, double *jArray, double *xyz) {
	int ret = 99;
	printf("Newton Augmented PointCenterSumStops\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	ret = getQps_IKNewtonAugmentedPointCenterSumStops(qps, knArray, nrArray, jArray, xyz);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], ret);
	std::cout << std::endl;
}

void check_nlSinglePoint(double *qpsGoal, double *xyzGoal, double *qps, double *knArray, double *nlArray, double *jArray, double *xyz) {
	int ret = 99;
	printf("NLOpt Point\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	ret = getQps_IKNLOpt_xyz6A(qps, knArray, nlArray, jArray, xyz);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], ret);
	std::cout << std::endl;
}
void check_nlPointPhiPsi(double *qpsGoal, double *xyzGoal, double *ppGoal, double *qps, double *knArray, double *nlArray, double *jArray, double *xyz, double *pp) {
	int ret = 99;
	//printf("NLOPT PointPhiPsi using %2.f||||", nlArray[2]);
	printf("NLOPT PointPhiPsi using %2.f\n", nlArray[2]);
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] = pp[%8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2], ppGoal[0], ppGoal[1]);
	ret = getQps_IKNLOpt_xyzpp6A(qps, knArray, nlArray, jArray, xyz, pp);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] = pp[%8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], pp[0], pp[1], ret);
}
void check_nlPointCtrSumStops(double *qpsGoal, double *xyzGoal, double *qps, double *knArray, double *nlArray, double *jArray, double *xyz) {
	int ret = 99;
	printf("NLOpt PointCenterSum using %2.f\n", nlArray[2]);
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	ret = getQps_IKNLOpt_xyzCtrSum6A(qps, knArray, nlArray, jArray, xyz);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], ret);
	std::cout << std::endl;
}
void check_nlPointPhiPsiCtrSum(double *qpsGoal, double *xyzGoal, double *ppGoal, double *qps, double *knArray, double *nlArray, double *jArray, double *xyz, double *pp) {
	int ret = 99;
	printf("NLOPT PointPhiPsiCtrSum using %2.f\n", nlArray[2]);
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] = pp[%8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2], ppGoal[0], ppGoal[1]);
	ret = getQps_IKNLOpt_xyzCtrSum6A(qps, knArray, nlArray, jArray, xyz);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] = pp[%8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], pp[0], pp[1], ret);

	std::cout << std::endl;
}
void check_nlPointMinChange() {
	int ret = 99;
	qps0[0] = -.548044; qps0[1] = -.0912528; qps0[2] = .407782; qps0[3] = .988095; qps0[4] = 207.976;
	//qps0[0] = -.548044; qps0[1] = -.0912528; qps0[2] = .407782; qps0[3] = .988095; qps0[4] = 210;
	xyzGoal[0] = 1090.61; xyzGoal[1] = -33.8123; xyzGoal[2] = -28.7369;
	qpsGoal[0] = 0; qpsGoal[1] = 0; qpsGoal[2] = 0; qpsGoal[3] = 0; qpsGoal[4] = 0; // not  used
	double qpsLast0[5] = { -.548044, -.0912528, .407782, .988095, 207.976 };

	//output from test_kalmanBounded.m
	//qps0[0] = 0.00000000; qps0[1] = 0.00000000; qps0[2] = 0.00000000; qps0[3] = 0.00100000; qps0[4] = 10.00000000;
	//xyzGoal[0] = 1224.49587541; xyzGoal[1] = 33.88769847; xyzGoal[2] = 72.67085838;
	//qpsLast0[0] = 0.02000000; qpsLast0[1] = 0.00000000; qpsLast0[2] = 0.00000000; qpsLast0[3] = 01.00100000; qpsLast0[4] = 10.00000000;


	double qps[5]; qps[0] = qps0[0]; qps[1] = qps0[1]; qps[2] = qps0[2]; qps[3] = qps0[3]; qps[4] = qps0[4];
	double xyz[3]; xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	double qpsLast[5]; qpsLast[0] = qpsLast0[0]; qpsLast[1] = qpsLast0[1]; qpsLast[2] = qpsLast0[2]; qpsLast[3] = qpsLast0[3]; qpsLast[4] = qpsLast0[4];


	printf("NLOpt Point MinChange\n");
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	printf("Last: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qpsLast0[0], qpsLast0[1], qpsLast0[2], qpsLast0[3], qpsLast0[4]);
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps0[0], qps0[1], qps0[2], qps0[3], qps0[4]);
	//      getQps_IKNLOptPointMinChange(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *qpsLast)
	ret = getQps_IKNLOpt_xyzMinChange6A(qps, xyz, qpsLast, kn6A, nlArray, jArray);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], ret);
	


	std::cout << std::endl;
}
void check_nlPoint_kalmanX() {
	int ret = 99;

	qps0[0] = -.548044; qps0[1] = -.0912528; qps0[2] = .407782; qps0[3] = .988095; qps0[4] = 207.976;
	//qps0[0] = -.548044; qps0[1] = -.0912528; qps0[2] = .407782; qps0[3] = .988095; qps0[4] = 210;
	xyzGoal[0] = 1090.61; xyzGoal[1] = -33.8123; xyzGoal[2] = -28.7369;
	qpsGoal[0] = 0; qpsGoal[1] = 0; qpsGoal[2] = 0; qpsGoal[3] = 0; qpsGoal[4] = 0; // not  used
	double K11s[5] = { .0666424, .0267724, .0285114, .0226446, .001841150 };
	double K21s[5] = { .2149370, .0353657, .0400750, .0253523, .000169336 };
	double qpsLast0[5] = { -.548044, -.0912528, .407782, .988095, 207.976 };
	double qdsLast0[5] = { -.453517, -.0949798, .387482, .979786, 26.261 };
	double tsSec = .01;
	
	//output from test_kalmanBounded.m
	qps0[0] = 0.00000000; qps0[1] = 0.00000000; qps0[2] = 0.00000000; qps0[3] = 0.00100000; qps0[4] = 10.00000000;
	xyzGoal[0] = 1224.49587541; xyzGoal[1] = 33.88769847; xyzGoal[2] = 72.67085838;
	K11s[0] = 0.54852819; K11s[1] = 0.54852819; K11s[2] = 0.54852819; K11s[3] = 0.54852819; K11s[4] = 0.02483491;
	K21s[0] = 2.12478200; K21s[1] = 2.12478200; K21s[2] = 2.12478200; K21s[3] = 2.12478200; K21s[4] = 0.00312272;
	qpsLast0[0] = 0.00000000; qpsLast0[1] = 0.00000000; qpsLast0[2] = 0.00000000; qpsLast0[3] = 0.00100000; qpsLast0[4] = 10.00000000;
	qdsLast0[0] = 0.01000000; qdsLast0[1] = -0.01000000; qdsLast0[2] = 0.01000000; qdsLast0[3] = 0.05000000; qdsLast0[4] = 10.00000000;
	//tsSec = .1;


	double qps[5]; qps[0] = qps0[0]; qps[1] = qps0[1]; qps[2] = qps0[2]; qps[3] = qps0[3]; qps[4] = qps0[4];
	double xyz[3]; xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	double qpsLast[5]; qpsLast[0] = qpsLast0[0]; qpsLast[1] = qpsLast0[1]; qpsLast[2] = qpsLast0[2]; qpsLast[3] = qpsLast0[3]; qpsLast[4] = qpsLast0[4];
	double qdsLast[5]; qdsLast[0] = qdsLast0[0]; qdsLast[1] = qdsLast0[1]; qdsLast[2] = qdsLast0[2]; qdsLast[3] = qdsLast0[3]; qdsLast[4] = qdsLast0[4];

	printf("NLOpt Point MinChangeKalman\n");
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps0[0], qps0[1], qps0[2], qps0[3], qps0[4]);
	for(int i = 0; i < 100; i++){
		xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
		//for (int j = 0; j < 5; j++) { std::cout << K21s[j] << " "; } std::cout << std::endl;
		ret = getQps_IKNLOpt_xyzKalX6A(
			qps,     //nl opt starts from here | converged qps
			xyz,     //goal in xyz | converged xyz
			qpsLast, //previous qps Kalman state | Kalman state for converged qps, == qps
			qdsLast, //previous qds Kalman state | Kalman state for converged qps
			kn6A,
			nlArray,
			jArray,
			K11s,
			K21s,
			tsSec);  //Kalman timestep
		printf("Cnvg: %d qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] ret=%d\n", i, qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], ret);
	}
	
	printf("Last: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qpsLast0[0], qpsLast0[1], qpsLast0[2], qpsLast0[3], qpsLast0[4]);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qpsLast[0], qpsLast[1], qpsLast[2], qpsLast[3], qpsLast[4]);
	printf("Last: qds[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qdsLast0[0], qdsLast0[1], qdsLast0[2], qdsLast0[3], qdsLast0[4]);
	printf("Cnvg: qds[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qdsLast[0], qdsLast[1], qdsLast[2], qdsLast[3], qdsLast[4]);
}
void check_nlPoint_kalmanStdBounds() {
	int ret = 99;

	//qps0[0] = -.548044; qps0[1] = -.0912528; qps0[2] = .407782; qps0[3] = .988095; qps0[4] = 207.976;
	xyzGoal[0] = 1290.61; xyzGoal[1] = -33.8123; xyzGoal[2] = -28.7369;
	qpsGoal[0] = 0; qpsGoal[1] = 0; qpsGoal[2] = 0; qpsGoal[3] = 0; qpsGoal[4] = 0; // not  used
	double K11s[5] = { .0666424, .0267724, .0285114, .0226446, .001841150 };
	double K21s[5] = { .2149370, .0353657, .0400750, .0253523, .000169336 };
	double stds[5] = { 10, 10, 10, 10, 1000 };
	double qpsLast0[5] = { -.548044, -.0912528, .407782, .988095, 207.976 };
	double qdsLast0[5] = { -.453517, -.0949798, .387482, .979786, 26.261 };
	double tsSec = .1;

	//output from test_kalmanBounded.m
	qps0[0] = 0.00000000; qps0[1] = 0.00000000; qps0[2] = 0.00000000; qps0[3] = 0.00100000; qps0[4] = 10.00000000;
	xyzGoal[0] = 1024.49587541; xyzGoal[1] = 33.88769847; xyzGoal[2] = 72.67085838;
	K11s[0] = 0.54852819; K11s[1] = 0.54852819; K11s[2] = 0.54852819; K11s[3] = 0.54852819; K11s[4] = 0.2483491;
	K21s[0] = 2.12478200; K21s[1] = 2.12478200; K21s[2] = 2.12478200; K21s[3] = 2.12478200; K21s[4] = 0.312272;
	stds[0] = 0.10000000; stds[1] = 0.10000000; stds[2] = 0.10000000; stds[3] = 0.10000000; stds[4] = 1000.00000000;
	qpsLast0[0] = 0.00000000; qpsLast0[1] = 0.00000000; qpsLast0[2] = 0.00000000; qpsLast0[3] = 0.00100000; qpsLast0[4] = 10.00000000;
	qdsLast0[0] = 0.01000000; qdsLast0[1] = -0.01000000; qdsLast0[2] = 0.01000000; qdsLast0[3] = 0.05000000; qdsLast0[4] = 10.00000000;
	//tsSec = .1;

	double qps[5]; qps[0] = qps0[0]; qps[1] = qps0[1]; qps[2] = qps0[2]; qps[3] = qps0[3]; qps[4] = qps0[4];
	double xyz[3]; xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
	double qpsLast[5]; qpsLast[0] = qpsLast0[0]; qpsLast[1] = qpsLast0[1]; qpsLast[2] = qpsLast0[2]; qpsLast[3] = qpsLast0[3]; qpsLast[4] = qpsLast0[4];
	double qdsLast[5]; qdsLast[0] = qdsLast0[0]; qdsLast[1] = qdsLast0[1]; qdsLast[2] = qdsLast0[2]; qdsLast[3] = qdsLast0[3]; qdsLast[4] = qdsLast0[4];
	
	printf("NLOpt Point Kalman Standard Deviation Bounds\n");
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps0[0], qps0[1], qps0[2], qps0[3], qps0[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	for (int i = 0; i < 100; i++) {
		xyz[0] = xyzGoal[0]; xyz[1] = xyzGoal[1]; xyz[2] = xyzGoal[2];
 		ret = getQps_IKNLOpt_xyzKalBounds6A(
			qps,     //nl opt starts from here | converged qps
			xyz,     //goal in xyz | converged xyz
			qpsLast, //previous qps Kalman state | Kalman state for converged qps, == qps
			qdsLast, //previous qds Kalman state | Kalman state for converged qps
			kn6A,
			nlArray,
			jArray,
			K11s,
			K21s,
			stds,
			tsSec);  //Kalman timestep

		printf("%4d: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] ret=%d\n", i, qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], ret);
		//printf("Cnvg: qpL[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qpsLast[0], qpsLast[1], qpsLast[2], qpsLast[3], qpsLast[4]);
	}
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2]);
	printf("Last: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qpsLast0[0], qpsLast0[1], qpsLast0[2], qpsLast0[3], qpsLast0[4]);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qpsLast[0], qpsLast[1], qpsLast[2], qpsLast[3], qpsLast[4]);
	printf("Last: qds[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qdsLast0[0], qdsLast0[1], qdsLast0[2], qdsLast0[3], qdsLast0[4]);
	printf("Cnvg: qds[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qdsLast[0], qdsLast[1], qdsLast[2], qdsLast[3], qdsLast[4]);
}

void check_iknlopt_xyzpp11A(){
	int ret = 99;

	double qps0[5];
	double qpsGoal[5], xyzGoal[3], ppGoal[2];
	double qps[5], xyz[3], pp[2];

	//initialize
	for (int i = 0; i < 5; i++) {
		qps0[i] = 0;
		qpsGoal[i] = 0;
		qps[i] = 0;
	}
	for (int i = 0; i < 3; i++) {
		xyz[i] = 0;
		xyzGoal[i] = 0;
	}
	ppGoal[0] = 0; ppGoal[1] = 0; pp[0] = 0; pp[1] = 0;

	//set goals
	qpsGoal[0] = .1; qpsGoal[1] = -.5; qpsGoal[2] = .3; qpsGoal[3] = 1.5; qpsGoal[4] = 43;
	ret = get11APoint(qpsGoal, kn11A, xyzGoal);
	ret = get11APhiPsi(qpsGoal, kn11A, ppGoal);

	//get fval at initial spot
	//double fval = ret = getFunVal_xyzpp11A(qpsGoal, kn11A, nlArray, jArray, xyzGoal, ppGoal, &fval);
	//printf("fval @qpsGoal = %f\n", fval);
	
	//solve IK
	printf("NLOPT xyzpp11A using %2.f\n", nlArray[2]);
	printf("Init: qps[%8.3f %8.3f %8.3f %8.3f %8.3f]\n", qps0[0], qps0[1], qps0[2], qps0[3], qps0[4]);
	printf("Goal: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] = pp[%8.3f %8.3f]\n", qpsGoal[0], qpsGoal[1], qpsGoal[2], qpsGoal[3], qpsGoal[4], xyzGoal[0], xyzGoal[1], xyzGoal[2], ppGoal[0], ppGoal[1]);
	for (int i = 0; i < 5; i++) { qps[i] = qps0[i]; }
	for (int i = 0; i < 3; i++) { xyz[i] = xyzGoal[i]; }
	pp[0] = ppGoal[0]; pp[1] = ppGoal[1];
	//printf("    : qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] = pp[%8.3f %8.3f]\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], pp[0], pp[1]);
	//printf("NL  : "); for (int i = 0; i < 6; i++) { printf("%f ", nlArray[i]); } printf("\n");
	//printf("JN  : "); for (int i = 0; i < 10; i++) { printf("%f ", jArray[i]); } printf("\n");
	ret = getQps_IKNLOpt_xyzpp11A(qps, kn11A, nlArray, jArray, xyz, pp);
	printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] = pp[%8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], pp[0], pp[1], ret);


	for (int i = 0; i < 100; i++) {
		xyz[0] += (double) i * 10;
		ret = getQps_IKNLOpt_xyzpp11A(qps, kn11A, nlArray, jArray, xyz, pp);
		printf("Cnvg: qps[%8.3f %8.3f %8.3f %8.3f %8.3f] = xyz[%8.3f %8.3f %8.3f] = pp[%8.3f %8.3f] ret=%d\n", qps[0], qps[1], qps[2], qps[3], qps[4], xyz[0], xyz[1], xyz[2], pp[0], pp[1], ret);
	}
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
	
	
	//check_nlPointMinChange();
	//check_nlPoint_kalmanX();
	//check_nlPoint_kalmanStdBounds();
	check_iknlopt_xyzpp11A();

	std::cout << "\n\nPress Enter";
	std::getchar();
	return 0;
}
