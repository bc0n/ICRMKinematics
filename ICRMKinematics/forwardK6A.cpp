#include "forwardKinematics.h"
#include <iostream>

FwdK6A::FwdK6A() {
}
FwdK6A::FwdK6A(KINEMATICPARAMS6A inParams) {
	kinParams = inParams;
}
Eigen::Matrix4d FwdK6A::qps2H01(double *qps) {
	//printf("kinParams = [%f %f %f %f %f %f]\n", kinParams.cathL, kinParams.rz01, kinParams.tx01, kinParams.ty01, kinParams.tz01, kinParams.tx23);

	Eigen::Matrix4d Txyz;
	Txyz.setIdentity();
	Txyz(0, 3) = kinParams.tx01;
	Txyz(1, 3) = kinParams.ty01;
	Txyz(2, 3) = kinParams.tz01;

	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;
	R = Eigen::AngleAxisd(kinParams.rz01, Eigen::Vector3d::UnitZ())
		* Eigen::AngleAxisd(qps[0], Eigen::Vector3d::UnitX());
	Rxyz.setIdentity();
	Rxyz.block<3, 3>(0, 0) = R;

	Eigen::Matrix4d H01;
	H01.setIdentity();
	H01 = H01 * Txyz * Rxyz;

	return H01;
}
Eigen::Matrix4d FwdK6A::qps2H02(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;

	R = Eigen::AngleAxisd(qps[1], Eigen::Vector3d::UnitZ());
	Rxyz.setIdentity();
	Rxyz.block<3, 3>(0, 0) = R;

	return FwdK6A::qps2H01(qps) * Rxyz;
}
Eigen::Matrix4d FwdK6A::qps2H03(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;

	R = Eigen::AngleAxisd(qps[2], Eigen::Vector3d::UnitX());
	Rxyz.setIdentity();
	Rxyz(0, 3) = kinParams.tx23; //Tx
	Rxyz.block<3, 3>(0, 0) = R;

	return FwdK6A::qps2H02(qps) * Rxyz;
}
Eigen::Matrix4d FwdK6A::qps2H04(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;

	if (qps[3] < 1e-3) {
		qps[3] = 1e-3;
	}

	double r = kinParams.cathL / qps[3]; //radius of catheter arc
	double ty = r * (1 - cos(qps[3]));
	double tx = r * sin(qps[3]);

	R = Eigen::AngleAxisd(qps[3], Eigen::Vector3d::UnitZ());
	Rxyz.setIdentity();
	Rxyz(0, 3) = tx;
	Rxyz(1, 3) = ty;
	Rxyz.block<3, 3>(0, 0) = R;

	return FwdK6A::qps2H03(qps) * Rxyz;
}
Eigen::Matrix4d FwdK6A::qps2H05(double *qps) {
	Eigen::Matrix4d Rxyz;
	Eigen::Matrix4d H05;

	Rxyz.setIdentity();
	Rxyz(0, 3) = qps[4]; //Tx

	return FwdK6A::qps2H04(qps) * Rxyz;
}

