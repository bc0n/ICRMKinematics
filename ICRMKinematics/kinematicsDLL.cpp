#include "kinematicsDLL.h"
#include "kinematics_structs.h"
#include "forwardKinematics.h"
#include "taskDefinitions.h"
#include "ik_nlopt.h"
#include "ip_nlopt.h"
#include <iostream>


//if in the header these confuse matlab's dll import
KINEMATICPARAMS5A kinArray2Struct5A(double *kinArray);
KINEMATICPARAMS6A kinArray2Struct6A(double *kinArray);
KINEMATICPARAMS11A kinArray2Struct11A(double *kinArray);
NLOPTPARAMS nlArray2Struct(double *nlArray);
JOINTLIMITS jntArray2Struct(double *jntArray);


// wrap forwardK
DLLIMPORT int get5AH01(double *qps, double *kinArray, double *arrayH01) {
	int ret = -1;

	FwdK5A fk(kinArray2Struct5A(kinArray));
	Eigen::Matrix4d H01 = fk.qps2H01(qps);
	arrayH01[0] = H01(0, 0);
	arrayH01[1] = H01(1, 0);
	arrayH01[2] = H01(2, 0);
	arrayH01[3] = H01(0, 1);
	arrayH01[4] = H01(1, 1);
	arrayH01[5] = H01(2, 1);
	arrayH01[6] = H01(0, 2);
	arrayH01[7] = H01(1, 2);
	arrayH01[8] = H01(2, 2);
	arrayH01[9] = H01(0, 3);
	arrayH01[10] = H01(1, 3);
	arrayH01[11] = H01(2, 3);

	ret = 0;

	return ret;
}
DLLIMPORT int get5AH02(double *qps, double *kinArray, double *arrayH02) {
	int ret = -1;

	FwdK5A fk(kinArray2Struct5A(kinArray));
	Eigen::Matrix4d H02 = fk.qps2H02(qps);

	arrayH02[0] = H02(0, 0);
	arrayH02[1] = H02(1, 0);
	arrayH02[2] = H02(2, 0);
	arrayH02[3] = H02(0, 1);
	arrayH02[4] = H02(1, 1);
	arrayH02[5] = H02(2, 1);
	arrayH02[6] = H02(0, 2);
	arrayH02[7] = H02(1, 2);
	arrayH02[8] = H02(2, 2);
	arrayH02[9] = H02(0, 3);
	arrayH02[10] = H02(1, 3);
	arrayH02[11] = H02(2, 3);

	ret = 0;

	return ret;
}
DLLIMPORT int get5AH03(double *qps, double *kinArray, double *arrayH03) {
	int ret = -1;

	FwdK5A fk(kinArray2Struct5A(kinArray));
	Eigen::Matrix4d H03 = fk.qps2H03(qps);

	arrayH03[0] = H03(0, 0);
	arrayH03[1] = H03(1, 0);
	arrayH03[2] = H03(2, 0);
	arrayH03[3] = H03(0, 1);
	arrayH03[4] = H03(1, 1);
	arrayH03[5] = H03(2, 1);
	arrayH03[6] = H03(0, 2);
	arrayH03[7] = H03(1, 2);
	arrayH03[8] = H03(2, 2);
	arrayH03[9] = H03(0, 3);
	arrayH03[10] = H03(1, 3);
	arrayH03[11] = H03(2, 3);

	ret = 0;

	return ret;
}
DLLIMPORT int get5AH04(double *qps, double *kinArray, double *arrayH04) {
	int ret = -1;

	FwdK5A fk(kinArray2Struct5A(kinArray));
	Eigen::Matrix4d H04 = fk.qps2H04(qps);

	arrayH04[0] = H04(0, 0);
	arrayH04[1] = H04(1, 0);
	arrayH04[2] = H04(2, 0);
	arrayH04[3] = H04(0, 1);
	arrayH04[4] = H04(1, 1);
	arrayH04[5] = H04(2, 1);
	arrayH04[6] = H04(0, 2);
	arrayH04[7] = H04(1, 2);
	arrayH04[8] = H04(2, 2);
	arrayH04[9] = H04(0, 3);
	arrayH04[10] = H04(1, 3);
	arrayH04[11] = H04(2, 3);

	ret = 0;

	return ret;
}
DLLIMPORT int get5AH05(double *qps, double *kinArray, double *arrayH05) {
	int ret = -1;

	FwdK5A fk(kinArray2Struct5A(kinArray));
	Eigen::Matrix4d H05 = fk.qps2H05(qps);

	arrayH05[0] = H05(0, 0);
	arrayH05[1] = H05(1, 0);
	arrayH05[2] = H05(2, 0);
	arrayH05[3] = H05(0, 1);
	arrayH05[4] = H05(1, 1);
	arrayH05[5] = H05(2, 1);
	arrayH05[6] = H05(0, 2);
	arrayH05[7] = H05(1, 2);
	arrayH05[8] = H05(2, 2);
	arrayH05[9] = H05(0, 3);
	arrayH05[10] = H05(1, 3);
	arrayH05[11] = H05(2, 3);

	ret = 0;
	
	return ret;
}

DLLIMPORT int get6AH01(double *qps, double *kinArray, double *arrayH01) {
	int ret = -1;
	
	FwdK6A fk(kinArray2Struct6A(kinArray));
	Eigen::Matrix4d H01 = fk.qps2H01(qps); // non-static usage;
	//Eigen::Matrix4d H01;
	//H01 = FwdK::qps2H01(qps, params); //now static 

	//arrayH01 = H01.data(); lv can't read it..
	arrayH01[0] = H01(0, 0);
	arrayH01[1] = H01(1, 0);
	arrayH01[2] = H01(2, 0);
	arrayH01[3] = H01(0, 1);
	arrayH01[4] = H01(1, 1);
	arrayH01[5] = H01(2, 1);
	arrayH01[6] = H01(0, 2);
	arrayH01[7] = H01(1, 2);
	arrayH01[8] = H01(2, 2);
	arrayH01[9] = H01(0, 3);
	arrayH01[10] = H01(1, 3);
	arrayH01[11] = H01(2, 3);

	ret = 0;

	return ret;
}
DLLIMPORT int get6AH02(double *qps, double *kinArray, double *arrayH02) {
	int ret = -1;

	FwdK6A fk(kinArray2Struct6A(kinArray));
	Eigen::Matrix4d H02 = fk.qps2H02(qps);

	arrayH02[0] = H02(0, 0);
	arrayH02[1] = H02(1, 0);
	arrayH02[2] = H02(2, 0);
	arrayH02[3] = H02(0, 1);
	arrayH02[4] = H02(1, 1);
	arrayH02[5] = H02(2, 1);
	arrayH02[6] = H02(0, 2);
	arrayH02[7] = H02(1, 2);
	arrayH02[8] = H02(2, 2);
	arrayH02[9] = H02(0, 3);
	arrayH02[10] = H02(1, 3);
	arrayH02[11] = H02(2, 3);

	ret = 0;

	return ret;
}
DLLIMPORT int get6AH03(double *qps, double *kinArray, double *arrayH03) {
	int ret = -1;

	FwdK6A fk(kinArray2Struct6A(kinArray));
	Eigen::Matrix4d H03 = fk.qps2H03(qps);

	arrayH03[0] = H03(0, 0);
	arrayH03[1] = H03(1, 0);
	arrayH03[2] = H03(2, 0);
	arrayH03[3] = H03(0, 1);
	arrayH03[4] = H03(1, 1);
	arrayH03[5] = H03(2, 1);
	arrayH03[6] = H03(0, 2);
	arrayH03[7] = H03(1, 2);
	arrayH03[8] = H03(2, 2);
	arrayH03[9] = H03(0, 3);
	arrayH03[10] = H03(1, 3);
	arrayH03[11] = H03(2, 3);

	ret = 0;

	return ret;
}
DLLIMPORT int get6AH04(double *qps, double *kinArray, double *arrayH04) {
	int ret = -1;

	FwdK6A fk(kinArray2Struct6A(kinArray));
	Eigen::Matrix4d H04 = fk.qps2H04(qps);

	arrayH04[0] = H04(0, 0);
	arrayH04[1] = H04(1, 0);
	arrayH04[2] = H04(2, 0);
	arrayH04[3] = H04(0, 1);
	arrayH04[4] = H04(1, 1);
	arrayH04[5] = H04(2, 1);
	arrayH04[6] = H04(0, 2);
	arrayH04[7] = H04(1, 2);
	arrayH04[8] = H04(2, 2);
	arrayH04[9] = H04(0, 3);
	arrayH04[10] = H04(1, 3);
	arrayH04[11] = H04(2, 3);

	ret = 0;

	return ret;
}
DLLIMPORT int get6AH05(double *qps, double *kinArray, double *arrayH05) {
	FwdK6A fk(kinArray2Struct6A(kinArray));
	Eigen::Matrix4d H05 = fk.qps2H05(qps);
	arrayH05[0] = H05(0, 0);
	arrayH05[1] = H05(1, 0);
	arrayH05[2] = H05(2, 0);
	arrayH05[3] = H05(0, 1);
	arrayH05[4] = H05(1, 1);
	arrayH05[5] = H05(2, 1);
	arrayH05[6] = H05(0, 2);
	arrayH05[7] = H05(1, 2);
	arrayH05[8] = H05(2, 2);
	arrayH05[9] = H05(0, 3);
	arrayH05[10] = H05(1, 3);
	arrayH05[11] = H05(2, 3);
	return 0;
}

DLLIMPORT int get11AH01(double *qps, double *kinArray, double *arrayH) {
	FwdK11A fk(kinArray2Struct11A(kinArray));
	Eigen::Matrix4d H = fk.qps2H01(qps);

	arrayH[0] = H(0, 0);
	arrayH[1] = H(1, 0);
	arrayH[2] = H(2, 0);
	arrayH[3] = H(0, 1);
	arrayH[4] = H(1, 1);
	arrayH[5] = H(2, 1);
	arrayH[6] = H(0, 2);
	arrayH[7] = H(1, 2);
	arrayH[8] = H(2, 2);
	arrayH[9] = H(0, 3);
	arrayH[10] = H(1, 3);
	arrayH[11] = H(2, 3);
	return 0;
}
DLLIMPORT int get11AH02(double *qps, double *kinArray, double *arrayH) {
	FwdK11A fk(kinArray2Struct11A(kinArray));
	Eigen::Matrix4d H = fk.qps2H02(qps);

	arrayH[0] = H(0, 0);
	arrayH[1] = H(1, 0);
	arrayH[2] = H(2, 0);
	arrayH[3] = H(0, 1);
	arrayH[4] = H(1, 1);
	arrayH[5] = H(2, 1);
	arrayH[6] = H(0, 2);
	arrayH[7] = H(1, 2);
	arrayH[8] = H(2, 2);
	arrayH[9] = H(0, 3);
	arrayH[10] = H(1, 3);
	arrayH[11] = H(2, 3);
	return 0;
}
DLLIMPORT int get11AH03(double *qps, double *kinArray, double *arrayH) {
	FwdK11A fk(kinArray2Struct11A(kinArray));
	Eigen::Matrix4d H = fk.qps2H03(qps);

	arrayH[0] = H(0, 0);
	arrayH[1] = H(1, 0);
	arrayH[2] = H(2, 0);
	arrayH[3] = H(0, 1);
	arrayH[4] = H(1, 1);
	arrayH[5] = H(2, 1);
	arrayH[6] = H(0, 2);
	arrayH[7] = H(1, 2);
	arrayH[8] = H(2, 2);
	arrayH[9] = H(0, 3);
	arrayH[10] = H(1, 3);
	arrayH[11] = H(2, 3);
	return 0;
}
DLLIMPORT int get11AH04(double *qps, double *kinArray, double *arrayH) {
	FwdK11A fk(kinArray2Struct11A(kinArray));
	Eigen::Matrix4d H = fk.qps2H04(qps);

	arrayH[0] = H(0, 0);
	arrayH[1] = H(1, 0);
	arrayH[2] = H(2, 0);
	arrayH[3] = H(0, 1);
	arrayH[4] = H(1, 1);
	arrayH[5] = H(2, 1);
	arrayH[6] = H(0, 2);
	arrayH[7] = H(1, 2);
	arrayH[8] = H(2, 2);
	arrayH[9] = H(0, 3);
	arrayH[10] = H(1, 3);
	arrayH[11] = H(2, 3);
	return 0;
}
DLLIMPORT int get11AH05(double *qps, double *kinArray, double *arrayH) {
	FwdK11A fk(kinArray2Struct11A(kinArray));
	Eigen::Matrix4d H = fk.qps2H05(qps);

	arrayH[0] = H(0, 0);
	arrayH[1] = H(1, 0);
	arrayH[2] = H(2, 0);
	arrayH[3] = H(0, 1);
	arrayH[4] = H(1, 1);
	arrayH[5] = H(2, 1);
	arrayH[6] = H(0, 2);
	arrayH[7] = H(1, 2);
	arrayH[8] = H(2, 2);
	arrayH[9] = H(0, 3);
	arrayH[10] = H(1, 3);
	arrayH[11] = H(2, 3);
	return 0;
}

// wrap taskDefinitions
DLLIMPORT int getTask5A_xyz(double *qps, double *kinArray, double *xyz) {
	FwdK5A fk5a(kinArray2Struct5A(kinArray));
	TaskXYZ<FwdK5A> task(fk5a);
	task.qps2task(qps, xyz);
	return 0;
}
DLLIMPORT int getTask6A_xyz(double *qps, double *kinArray, double *xyz) {
	FwdK6A fk6a(kinArray2Struct6A(kinArray));
	TaskXYZ<FwdK6A> task(fk6a);
	task.qps2task(qps, xyz);
	return 0;
}
DLLIMPORT int getTask11A_xyz(double *qps, double *kinArray, double *xyz) {
	FwdK11A fk(kinArray2Struct11A(kinArray));
	TaskXYZ<FwdK11A> task(fk);
	task.qps2task(qps, xyz);
	return 0;
}
DLLIMPORT int getTask6A_phiPsi(double *qps, double *kinArray, double *pp) {
	FwdK6A fk6a(kinArray2Struct6A(kinArray));
	TaskPhiPsi<FwdK6A> task(fk6a);
	task.qps2task(qps, pp);
	return 0;
}
DLLIMPORT int getTask11A_phiPsi(double *qps, double *kinArray, double *pp) {
	FwdK11A fk11a(kinArray2Struct11A(kinArray));
	TaskPhiPsi<FwdK11A> task(fk11a);
	task.qps2task(qps, pp);
	return 0;
}
DLLIMPORT int getTask5A_xyzuxuyuz(double *qps, double *kinArray, double *xyz, double *uxyz) {
	FwdK5A fk5a(kinArray2Struct5A(kinArray));
	TaskXYZUxUyUz<FwdK5A> task(fk5a);
	double x[6] = { xyz[0],xyz[1],xyz[2], uxyz[0],uxyz[1],uxyz[2] };
	task.qps2task(qps, x);
	xyz[0] = x[0]; xyz[1] = x[1]; xyz[2] = x[2];
	uxyz[0] = x[3]; uxyz[1] = x[4]; uxyz[2] = x[5];
	return 0;
}
DLLIMPORT int getTask6A_xyzuxuyuz(double *qps, double *kinArray, double *xyz, double *uxyz) {
	FwdK6A fk6a(kinArray2Struct6A(kinArray));
	TaskXYZUxUyUz<FwdK6A> task(fk6a);
	double x[6] = { xyz[0],xyz[1],xyz[2], uxyz[0],uxyz[1],uxyz[2] };
	task.qps2task(qps, x);
	xyz[0] = x[0]; xyz[1] = x[1]; xyz[2] = x[2];
	uxyz[0] = x[3]; uxyz[1] = x[4]; uxyz[2] = x[5];
	return 0;
}
DLLIMPORT int getTask11A_xyzuxuyuz(double *qps, double *kinArray, double *xyz, double *uxyz) {
	FwdK11A fk11a(kinArray2Struct11A(kinArray));
	TaskXYZUxUyUz<FwdK11A> task(fk11a);
	double x[6] = { xyz[0],xyz[1],xyz[2], uxyz[0],uxyz[1],uxyz[2] };
	task.qps2task(qps, x);
	xyz[0] = x[0]; xyz[1] = x[1]; xyz[2] = x[2];
	uxyz[0] = x[3]; uxyz[1] = x[4]; uxyz[2] = x[5];
	return 0;
}

// wrap funs
DLLIMPORT int fun_qp0_xyzuxuyuz5A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin) {
	InvP_nlopt<KINEMATICPARAMS5A, FwdK5A> ip;
	ip.funQp(nSamps, stackedQ, stackedX, stackedU, qp0, kn0, fmin);
	return 0;
}
DLLIMPORT int fun_qp0_xyzuxuyuz6A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin) {
	InvP_nlopt<KINEMATICPARAMS6A, FwdK6A> ip;
	ip.funQp(nSamps, stackedQ, stackedX, stackedU, qp0, kn0, fmin);
	return 0;
}
DLLIMPORT int fun_qp0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin) {
	InvP_nlopt<KINEMATICPARAMS11A, FwdK11A> ip;
	ip.funQp(nSamps, stackedQ, stackedX, stackedU, qp0, kn0, fmin);
	return 0;
}

DLLIMPORT int fun_kn0_xyzuxuyuz5A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *kns0, double *fmin) {
	InvP_nlopt<KINEMATICPARAMS5A, FwdK5A> ip;
	ip.funKn(nSamps, stackedQ, stackedX, stackedU, kns0, fmin);
	return 0;
}
DLLIMPORT int fun_kn0_xyzuxuyuz6A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *kns0, double *fmin) {
	InvP_nlopt<KINEMATICPARAMS6A, FwdK6A> ip;
	ip.funKn(nSamps, stackedQ, stackedX, stackedU, kns0, fmin);
	return 0;
}
DLLIMPORT int fun_kn0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *kns0, double *fmin) {
	InvP_nlopt<KINEMATICPARAMS11A, FwdK11A> ip;
	ip.funKn(nSamps, stackedQ, stackedX, stackedU, kns0, fmin);
	return 0;
}

DLLIMPORT int fun_qp0kn0_xyzuxuyuz5A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *kns0, double *fmin) {
	InvP_nlopt<KINEMATICPARAMS5A, FwdK5A> ip;
	ip.funQpKn(nSamps, stackedQ, stackedX, stackedU, qps0, kns0, fmin);
	return 0;
}
DLLIMPORT int fun_qp0kn0_xyzuxuyuz6A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *kns0, double *fmin) {
	InvP_nlopt<KINEMATICPARAMS6A, FwdK6A> ip;
	ip.funQpKn(nSamps, stackedQ, stackedX, stackedU, qps0, kns0, fmin);
	return 0;
}
DLLIMPORT int fun_qp0kn0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *kns0, double *fmin) {
	InvP_nlopt<KINEMATICPARAMS11A, FwdK11A> ip;
	ip.funQpKn(nSamps, stackedQ, stackedX, stackedU, qps0, kns0, fmin);
	return 0;
}


// wrap inverse kinematic solvers
DLLIMPORT int estimate_qps_xyz5A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *fmin) {
	int ret = -1;

	FwdK5A fk5a(kinArray2Struct5A(kinArray));
	InvK_nlopt<TaskXYZ<FwdK5A>, FwdK5A> ik(fk5a, jntArray2Struct(jntArray), nlArray2Struct(nlArray));
	ret = ik.solve(qps, xyz, fmin);

	return ret;
}
DLLIMPORT int estimate_qps_xyz6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *fmin) {
	int ret = -1;
	
	FwdK6A fk6a(kinArray2Struct6A(kinArray));
	InvK_nlopt<TaskXYZ<FwdK6A>, FwdK6A> ik(fk6a, jntArray2Struct(jntArray), nlArray2Struct(nlArray));
	ret = ik.solve(qps, xyz, fmin);
	
	return ret;
}
DLLIMPORT int estimate_qps_xyz11A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *fmin) {
	int ret = -1;

	FwdK11A fk11a(kinArray2Struct11A(kinArray));
	InvK_nlopt<TaskXYZ<FwdK11A>, FwdK11A> ik(fk11a, jntArray2Struct(jntArray), nlArray2Struct(nlArray));
	ret = ik.solve(qps, xyz, fmin);

	return ret;
}
DLLIMPORT int estimate_qps_xyzuxuyuz5A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *uxyz, double *fmin) {
	int ret = -1;

	FwdK5A fk5a(kinArray2Struct5A(kinArray));
	InvK_nlopt<TaskXYZUxUyUz<FwdK5A>, FwdK5A> ik(fk5a, jntArray2Struct(jntArray), nlArray2Struct(nlArray));
	double x[6] = { xyz[0],xyz[1],xyz[2], uxyz[0],uxyz[1],uxyz[2] };
	ret = ik.solve(qps, x, fmin);
	xyz[0] = x[0]; xyz[1] = x[1]; xyz[2] = x[2];
	uxyz[0] = x[3]; uxyz[1] = x[4]; uxyz[2] = x[5];

	return ret;
}
DLLIMPORT int estimate_qps_xyzuxuyuz6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *uxyz, double *fmin) {
	int ret = -1;

	FwdK6A fk6a(kinArray2Struct6A(kinArray));
	InvK_nlopt<TaskXYZUxUyUz<FwdK6A>, FwdK6A> ik(fk6a, jntArray2Struct(jntArray), nlArray2Struct(nlArray));
	double x[6] = { xyz[0],xyz[1],xyz[2], uxyz[0],uxyz[1],uxyz[2] };
	ret = ik.solve(qps, x, fmin);
	xyz[0] = x[0]; xyz[1] = x[1]; xyz[2] = x[2];
	uxyz[0] = x[3]; uxyz[1] = x[4]; uxyz[2] = x[5];

	return ret;
}
DLLIMPORT int estimate_qps_xyzuxuyuz11A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *uxyz, double *fmin) {
	int ret = -1;

	FwdK11A fk11a(kinArray2Struct11A(kinArray));
	InvK_nlopt<TaskXYZUxUyUz<FwdK11A>, FwdK11A> ik(fk11a, jntArray2Struct(jntArray), nlArray2Struct(nlArray));
	double x[6] = { xyz[0],xyz[1],xyz[2], uxyz[0],uxyz[1],uxyz[2] };
	ret = ik.solve(qps, x, fmin);
	xyz[0] = x[0]; xyz[1] = x[1]; xyz[2] = x[2];
	uxyz[0] = x[3]; uxyz[1] = x[4]; uxyz[2] = x[5];

	return ret;
}

// wrap initial joint angle estimators
DLLIMPORT int estimate_qp0_xyzuxuyuz5A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *kn0, double *nlArray, double *fmin) {
	KINEMATICPARAMS5A knup; //default ok, not used
	KINEMATICPARAMS5A kndn;
	JOINTLIMITS qp0Lims = jntArray2Struct(q0Lims);
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS5A, FwdK5A> ip(kinArray2Struct5A(kn0),kinArray2Struct5A(kn0), qp0Lims, nlParams);//knup & dn aren't used
	return ip.estimateQp(nSamps, stackedQ, stackedX, stackedU, qp0, kn0, fmin);
}
DLLIMPORT int estimate_qp0_xyzuxuyuz6A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *kn0, double *nlArray, double *fmin) {
	KINEMATICPARAMS6A knup; //default ok, not used
	KINEMATICPARAMS6A kndn;
	JOINTLIMITS qp0Lims = jntArray2Struct(q0Lims);
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS6A, FwdK6A> ip(knup, kndn, qp0Lims, nlParams);
	return ip.estimateQp(nSamps, stackedQ, stackedX, stackedU, qp0, kn0, fmin);
}
DLLIMPORT int estimate_qp0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *kn0, double *nlArray, double *fmin) {
	KINEMATICPARAMS11A knup;
	KINEMATICPARAMS11A kndn;
	JOINTLIMITS qp0Lims = jntArray2Struct(q0Lims);
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS11A, FwdK11A> ip(knup, kndn, qp0Lims, nlParams);
	return ip.estimateQp(nSamps, stackedQ, stackedX, stackedU, qp0, kn0, fmin);
}

// wrap inverse parameter solvers
DLLIMPORT int estimate_kn0_xyzuxuyuz5A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *knup, double *kndn, double *nlArray, double *fmin) {
	KINEMATICPARAMS5A kup = kinArray2Struct5A(knup);
	KINEMATICPARAMS5A kdn = kinArray2Struct5A(kndn);
	JOINTLIMITS qp0Lims; //default ok, not used
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS5A, FwdK5A> ip(kup, kdn, qp0Lims, nlParams);
	return ip.estimateKn(nSamps, stackedQ, stackedX, stackedU, kn0, fmin);
}
DLLIMPORT int estimate_kn0_xyzuxuyuz6A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *knup, double *kndn, double *nlArray, double *fmin) {
	KINEMATICPARAMS6A kup = kinArray2Struct6A(knup);
	KINEMATICPARAMS6A kdn = kinArray2Struct6A(kndn);
	JOINTLIMITS qp0Lims; //default ok, not used
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS6A, FwdK6A> ip(kup, kdn, qp0Lims, nlParams);
	return ip.estimateKn(nSamps, stackedQ, stackedX, stackedU, kn0, fmin);
}
DLLIMPORT int estimate_kn0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *knup, double *kndn, double *nlArray, double *fmin) {
	KINEMATICPARAMS11A kup = kinArray2Struct11A(knup);
	KINEMATICPARAMS11A kdn = kinArray2Struct11A(kndn);
	JOINTLIMITS qp0Lims; //default ok, not used
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS11A, FwdK11A> ip(kup, kdn, qp0Lims, nlParams);
	return ip.estimateKn(nSamps, stackedQ, stackedX, stackedU, kn0, fmin);
}

DLLIMPORT int estimate_kn0_xyzuxuyuz5A_subset(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *knup, double *kndn, bool *knSub, double *nlArray, double *fmin) {
	KINEMATICPARAMS5A kup = kinArray2Struct5A(knup);
	KINEMATICPARAMS5A kdn = kinArray2Struct5A(kndn);
	JOINTLIMITS qp0Lims; //default ok, not used
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS5A, FwdK5A> ip(kup, kdn, knSub, qp0Lims, nlParams);
	return ip.estimateKn(nSamps, stackedQ, stackedX, stackedU, kn0, fmin);
}
DLLIMPORT int estimate_kn0_xyzuxuyuz6A_subset(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *knup, double *kndn, bool *knSub, double *nlArray, double *fmin) {
	KINEMATICPARAMS6A kup = kinArray2Struct6A(knup);
	KINEMATICPARAMS6A kdn = kinArray2Struct6A(kndn);
	JOINTLIMITS qp0Lims; //default ok, not used
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS6A, FwdK6A> ip(kup, kdn, knSub, qp0Lims, nlParams);
	return ip.estimateKn(nSamps, stackedQ, stackedX, stackedU, kn0, fmin);
}
DLLIMPORT int estimate_kn0_xyzuxuyuz11A_subset(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *knup, double *kndn, bool *knSub, double *nlArray, double *fmin) {
	KINEMATICPARAMS11A kup = kinArray2Struct11A(knup);
	KINEMATICPARAMS11A kdn = kinArray2Struct11A(kndn);
	JOINTLIMITS qp0Lims; //default ok, not used
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS11A, FwdK11A> ip(kup, kdn, knSub, qp0Lims, nlParams);
	return ip.estimateKn(nSamps, stackedQ, stackedX, stackedU, kn0, fmin);
}

// wrap joint qp0 & kn0
DLLIMPORT int estimate_qp0kn0_xyzuxuyuz5A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *kn0, double *knup, double *kndn, double *nlArray, double *fmin) {
	
	KINEMATICPARAMS5A kup = kinArray2Struct5A(knup);
	KINEMATICPARAMS5A kdn = kinArray2Struct5A(kndn);
	JOINTLIMITS qp0Lims = jntArray2Struct(q0Lims);
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);
	
	InvP_nlopt<KINEMATICPARAMS5A, FwdK5A> ip(kup, kdn, qp0Lims, nlParams);
	return ip.estimateQpKn(nSamps, stackedQ, stackedX, stackedU, qp0, kn0, fmin);
}
DLLIMPORT int estimate_qp0kn0_xyzuxuyuz6A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *kn0, double *knup, double *kndn, double *nlArray, double *fmin) {
	KINEMATICPARAMS6A kup = kinArray2Struct6A(knup);
	KINEMATICPARAMS6A kdn = kinArray2Struct6A(kndn);
	JOINTLIMITS qp0Lims = jntArray2Struct(q0Lims);
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS6A, FwdK6A> ip(kup, kdn, qp0Lims, nlParams);
	return ip.estimateQpKn(nSamps, stackedQ, stackedX, stackedU, qp0, kn0, fmin);
}
DLLIMPORT int estimate_qp0kn0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *kn0, double *knup, double *kndn, double *nlArray, double *fmin) {
	KINEMATICPARAMS11A kup = kinArray2Struct11A(knup);
	KINEMATICPARAMS11A kdn = kinArray2Struct11A(kndn);
	JOINTLIMITS qp0Lims = jntArray2Struct(q0Lims);
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS11A, FwdK11A> ip(kup, kdn, qp0Lims, nlParams);
	return ip.estimateQpKn(nSamps, stackedQ, stackedX, stackedU, qp0, kn0, fmin);
}
DLLIMPORT int estimate_qp0kn0_xyzuxuyuz5A_subset(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *kn0, double *knup, double *kndn, bool *knSub, double *nlArray, double *fmin) {
	KINEMATICPARAMS5A kup = kinArray2Struct5A(knup);
	KINEMATICPARAMS5A kdn = kinArray2Struct5A(kndn);
	JOINTLIMITS qp0Lims = jntArray2Struct(q0Lims);
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS5A, FwdK5A> ip(kup, kdn, knSub, qp0Lims, nlParams);
	return ip.estimateQpKn(nSamps, stackedQ, stackedX, stackedU, qp0, kn0, fmin);
}
DLLIMPORT int estimate_qp0kn0_xyzuxuyuz6A_subset(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *kn0, double *knup, double *kndn, bool *knSub, double *nlArray, double *fmin) {
	KINEMATICPARAMS6A kup = kinArray2Struct6A(knup);
	KINEMATICPARAMS6A kdn = kinArray2Struct6A(kndn);
	JOINTLIMITS qp0Lims = jntArray2Struct(q0Lims);
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS6A, FwdK6A> ip(kup, kdn, knSub, qp0Lims, nlParams);
	return ip.estimateQpKn(nSamps, stackedQ, stackedX, stackedU, qp0, kn0, fmin);
}
DLLIMPORT int estimate_qp0kn0_xyzuxuyuz11A_subset(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *kn0, double *knup, double *kndn, bool *knSub, double *nlArray, double *fmin) {
	KINEMATICPARAMS11A kup = kinArray2Struct11A(knup);
	KINEMATICPARAMS11A kdn = kinArray2Struct11A(kndn);
	JOINTLIMITS qp0Lims = jntArray2Struct(q0Lims);
	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);

	InvP_nlopt<KINEMATICPARAMS11A, FwdK11A> ip(kup, kdn, knSub, qp0Lims, nlParams);
	return ip.estimateQpKn(nSamps, stackedQ, stackedX, stackedU, qp0, kn0, fmin);
}

KINEMATICPARAMS5A kinArray2Struct5A(double *kinArray) {
	KINEMATICPARAMS5A params;
	params.tx01 = kinArray[0];
	params.ty01 = kinArray[1];
	params.tz01 = kinArray[2];
	params.rz01 = kinArray[3];
	params.lCath = kinArray[4];
	return params;
}
KINEMATICPARAMS6A kinArray2Struct6A(double *kinArray) {
	KINEMATICPARAMS6A params;
	params.tx01 = kinArray[0];
	params.ty01 = kinArray[1];
	params.tz01 = kinArray[2];
	params.rz01 = kinArray[3];
	params.ry34 = kinArray[4];
	params.lCath = kinArray[5];
	return params;
}
KINEMATICPARAMS11A kinArray2Struct11A(double *kinArray) {
	KINEMATICPARAMS11A params;
	//tx01,ty01,tz01,ry01,rz01, ry34,rz34,kAlpha,eAlpha,cathL, ry45
	params.tx01 = kinArray[0];
	params.ty01 = kinArray[1];
	params.tz01 = kinArray[2];
	params.ry01 = kinArray[3];
	params.rz01 = kinArray[4];
	params.ry34 = kinArray[5];
	params.rz34 = kinArray[6];
	params.kAlpha = kinArray[7];
	params.eAlpha = kinArray[8];
	params.lCath = kinArray[9];
	params.ry45 = kinArray[10];
	return params;
}

NLOPTPARAMS nlArray2Struct(double *nlArray) {
	NLOPTPARAMS params;
	params.maxIts = (int)nlArray[0];
	params.maxTimeSec = nlArray[1];
	params.method = (nlMethod)(int)nlArray[2];
	params.minFunVal = nlArray[3];
	params.tolFunAbs = nlArray[4];
	params.tolXAbs = nlArray[5];
	return params;
}

//[down; up]
JOINTLIMITS jntArray2Struct(double *jntArray) {
	JOINTLIMITS jntLims;
	for (int i = 0; i < 5; i++) {
		jntLims.dn[i] = jntArray[i];
		jntLims.up[i] = jntArray[i + 5];
	}
	return jntLims;
}