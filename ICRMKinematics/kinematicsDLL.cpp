#include "kinematicsDLL.h"
#include "taskDefinitions.h"
#include <iostream>


// wrap forwardK
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
	int ret = -1;

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

	ret = 0;

	return ret;
}

DLLIMPORT int get11AH05(double *qps, double *kinArray, double *arrayH05) {
	int ret = -1;

	FwdK11A fk(kinArray2Struct11A(kinArray));
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

// wrap taskDefinitions
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

// wrap inverse kinematic solvers
DLLIMPORT int getQps_IKnlopt_xyz6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz) {
	int ret = -1;
	
	FwdK6A fk6a(kinArray2Struct6A(kinArray));
	InvK_nlopt<TaskXYZ<FwdK6A>, FwdK6A> ik(fk6a, jntArray2Struct(jntArray), nlArray2Struct(nlArray));
	ret = ik.solve(qps, xyz);
	
	return ret;
}
DLLIMPORT int getQps_IKnlopt_xyz11A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz) {
	int ret = -1;

	FwdK11A fk11a(kinArray2Struct11A(kinArray));
	InvK_nlopt<TaskXYZ<FwdK11A>, FwdK11A> ik(fk11a, jntArray2Struct(jntArray), nlArray2Struct(nlArray));
	ret = ik.solve(qps, xyz);

	return ret;
}
DLLIMPORT int getQps_IKnlopt_xyzuxuyuz6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *uxyz) {
	int ret = -1;

	FwdK6A fk6a(kinArray2Struct6A(kinArray));
	InvK_nlopt<TaskXYZUxUyUz<FwdK6A>, FwdK6A> ik(fk6a, jntArray2Struct(jntArray), nlArray2Struct(nlArray));
	double x[6] = { xyz[0],xyz[1],xyz[2], uxyz[0],uxyz[1],uxyz[2] };
	ret = ik.solve(qps, x);
	xyz[0] = x[0]; xyz[1] = x[1]; xyz[2] = x[2];
	uxyz[0] = x[3]; uxyz[4] = x[1]; uxyz[5] = x[2];

	return ret;
}
DLLIMPORT int getQps_IKnlopt_xyzuxuyuz11A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *uxyz) {
	int ret = -1;

	FwdK11A fk11a(kinArray2Struct11A(kinArray));
	InvK_nlopt<TaskXYZUxUyUz<FwdK11A>, FwdK11A> ik(fk11a, jntArray2Struct(jntArray), nlArray2Struct(nlArray));
	double x[6] = { xyz[0],xyz[1],xyz[2], uxyz[0],uxyz[1],uxyz[2] };
	ret = ik.solve(qps, x);
	xyz[0] = x[0]; xyz[1] = x[1]; xyz[2] = x[2];
	uxyz[0] = x[3]; uxyz[4] = x[1]; uxyz[5] = x[2];

	return ret;
}

//DLLIMPORT int getQps_IKNLOpt_xyzpp6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *phipsi) {
//	int ret = -1;
//
//	InvKNLOpt_xyzpp6A ik(kinArray2Struct6A(kinArray), nlArray2Struct(nlArray), jntArray2Struct(jntArray));
//	ret = ik.solve(qps, xyz, phipsi);
//
//	return ret;
//}
//DLLIMPORT int getQps_IKNLOpt_xyzCtrSum6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz) {
//	int ret = -1;
//
//	InvKNLOpt_xyzCenterSum6A ik(kinArray2Struct6A(kinArray), nlArray2Struct(nlArray), jntArray2Struct(jntArray));
//	ret = ik.solve(qps, xyz);
//
//	return ret;
//}
//DLLIMPORT int getQps_IKNLOpt_xyzppCtrSum6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *pp) {
//	int ret = -1;
//
//	InvKNLOpt_xyzppCenterSum6A ik(kinArray2Struct6A(kinArray), nlArray2Struct(nlArray), jntArray2Struct(jntArray));
//	ret = ik.solve(qps, xyz, pp);
//
//	return ret;
//}
//DLLIMPORT int getQps_IKNLOpt_xyzMinChange6A(double *qps, double *xyz, double *qpsLast, double *kinArray, double *nlArray, double *jntArray) {
//	int ret = -1;
//
//	InvKNLOpt_xyzMinChange6A ik(kinArray2Struct6A(kinArray), nlArray2Struct(nlArray), jntArray2Struct(jntArray));
//	ret = ik.solve(qps, xyz, qpsLast);
//
//	return ret;
//}
//DLLIMPORT int getQps_IKNLOpt_xyzKalX6A(double *qps, double *xyz, double *qpsLast, double *qdsLast, double *kinArray, double *nlArray, double *jntArray, double *K11s, double *K21s, double tsSec) {
//	int ret = -1;
//
//	InvKNLOpt_xyzKalX6A ik(kinArray2Struct6A(kinArray), nlArray2Struct(nlArray), jntArray2Struct(jntArray), tsSec);
//	ret = ik.solve(qps, xyz, qpsLast, qdsLast, K11s, K21s);
//
//	return ret;
//}
//DLLIMPORT int getQps_IKNLOpt_xyzKalBounds6A(double *qps, double *xyz, double *qpsLast, double *qdsLast, double *kinArray, double *nlArray, double *jntArray, double *K11s, double *K21s, double *stdDevs, double tsSec) {
//	int ret = -1;
//
//	InvKNLOpt_xyzKalBounds6A ik(kinArray2Struct6A(kinArray), nlArray2Struct(nlArray), jntArray2Struct(jntArray), tsSec);
//	ret = ik.solve(qps, xyz, qpsLast, qdsLast, K11s, K21s, stdDevs);
//
//	return ret;
//}
//DLLIMPORT int getQps_IKNLOpt_xyzpp11A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *phipsi) {
//	int ret = -1;
//
//	InvKNLOpt_xyzpp11A ik(kinArray2Struct11A(kinArray), nlArray2Struct(nlArray), jntArray2Struct(jntArray));
//	ret = ik.solve(qps, xyz, phipsi);
//
//	return ret;
//}
//DLLIMPORT int getFunVal_xyzpp11A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *phipsi, double *fval) {
//	int ret = -1;
//	InvKNLOpt_xyzpp11A ik(kinArray2Struct11A(kinArray), nlArray2Struct(nlArray), jntArray2Struct(jntArray));
//	ret = ik.getFval(qps, xyz, phipsi, fval);
//
//	return ret;
//}

// wrap inverse parameter solvers
//DLLIMPORT int funIP_xyzdotu11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *pms0, double *fmin) {
//	InvPNLOpt_xyzdotu11A ip;
//	ip.funIP_UX11A(nSamps, stackedQ, stackedU, stackedX, qps0, pms0, fmin);
//	return 0;
//}
//DLLIMPORT int estimatePmsQ_IPNLOpt_xyzdotu11A_assumeX0(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *k11up, double *k11dn, double *q0Lims, double *nlArray, double *fmin) {
//	int ret = -99;
//
//	KINEMATICPARAMS11A kp11up = kinArray2Struct11A(k11up);
//	KINEMATICPARAMS11A kp11dn = kinArray2Struct11A(k11dn);
//	JOINTLIMITS q0pLims = jntArray2Struct(q0Lims);
//	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);
//
//	InvPNLOpt_xyzdotu11A ip(kp11up, kp11dn, q0pLims, nlParams);
//	ret = ip.estimate(nSamps, stackedQ, stackedU, stackedX, fmin);
//	
//	return ret;
//}
//DLLIMPORT int estimatePmsQ_IPNLOpt_xyzdotu11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *k11up, double *k11dn, double *q0Lims, double *nlArray, double *qps0, double *kps0, double *fmin) {
//	int ret = -99;
//
//	KINEMATICPARAMS11A kp11up = kinArray2Struct11A(k11up);
//	KINEMATICPARAMS11A kp11dn = kinArray2Struct11A(k11dn);
//	JOINTLIMITS q0pLims = jntArray2Struct(q0Lims);
//	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);
//
//	InvPNLOpt_xyzdotu11A ip(kp11up, kp11dn, q0pLims, nlParams);
//	double x0[5 + 11];
//	for (int i = 0; i < 5; i++) { x0[i] = qps0[i]; }
//	for (int i = 0; i < 11; i++) { x0[i+5] = kps0[i]; }
//
//	ret = ip.estimate(nSamps, stackedQ, stackedU, stackedX, x0, fmin);
//
//	for (int i = 0; i < 5; i++) { qps0[i] = x0[i]; }
//	for (int i = 0; i < 11; i++) { kps0[i] = x0[i + 5]; }
//	return ret;
//}
//
//DLLIMPORT int funIP_xyzpp11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *pms0, double *fmin) {
//	InvPNLOpt_xyzpp11A ip;
//	ip.funIP_xyzpp11A(nSamps, stackedQ, stackedU, stackedX, qps0, pms0, fmin);
//	return 0;
//}
//DLLIMPORT int estimatePmsQ_IPNLOpt_xyzpp11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *k11up, double *k11dn, double *q0Lims, double *nlArray, double *qps0, double *kps0, double *fmin) {
//	int ret = -99;
//
//	KINEMATICPARAMS11A kp11up = kinArray2Struct11A(k11up);
//	KINEMATICPARAMS11A kp11dn = kinArray2Struct11A(k11dn);
//	JOINTLIMITS q0pLims = jntArray2Struct(q0Lims);
//	NLOPTPARAMS nlParams = nlArray2Struct(nlArray);
//
//	InvPNLOpt_xyzpp11A ip(kp11up, kp11dn, q0pLims, nlParams);
//	double x0[5 + 11];
//	for (int i = 0; i < 5; i++) { x0[i] = qps0[i]; }
//	for (int i = 0; i < 11; i++) { x0[i + 5] = kps0[i]; }
//
//	ret = ip.estimate(nSamps, stackedQ, stackedU, stackedX, x0, fmin);
//
//	for (int i = 0; i < 5; i++) { qps0[i] = x0[i]; }
//	for (int i = 0; i < 11; i++) { kps0[i] = x0[i + 5]; }
//	return ret;
//}


KINEMATICPARAMS6A kinArray2Struct6A(double *kinArray) {
	KINEMATICPARAMS6A params;
	params.cathL = kinArray[0];
	params.rz01 = kinArray[1];
	params.tx01 = kinArray[2];
	params.ty01 = kinArray[3];
	params.tz01 = kinArray[4];
	params.tx23 = kinArray[5];
	return params;
}
KINEMATICPARAMS11A kinArray2Struct11A(double *kinArray) {
	KINEMATICPARAMS11A params;
	//tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45
	params.tx01 = kinArray[0];
	params.ty01 = kinArray[1];
	params.tz01 = kinArray[2];
	params.ry01 = kinArray[3];
	params.rz01 = kinArray[4];
	params.tx23 = kinArray[5];
	params.ry34 = kinArray[6];
	params.rz34 = kinArray[7];
	params.cathL = kinArray[8];
	params.ry45 = kinArray[9];
	params.rz45 = kinArray[10];
	return params;
}

//NEWTONPARAMS nrArray2Struct(double *nrArray) {
//	NEWTONPARAMS params;
//	params.epslon = nrArray[0];
//	params.errTol = nrArray[1];
//	params.maxIts = nrArray[2];
//	params.fact01 = nrArray[3];
//	params.fact02 = nrArray[4];
//	return params;
//}

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

JOINTLIMITS jntArray2Struct(double *jntArray) {
	JOINTLIMITS jntLims;
	for (int i = 0; i < 5; i++) {
		jntLims.dn[i] = jntArray[i];
		jntLims.up[i] = jntArray[i + 5];
	}
	return jntLims;
}