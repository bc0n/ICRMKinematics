#include "forwardKinematics.h"
#include <iostream>

FwdK6A::FwdK6A() {
}
FwdK6A::FwdK6A(double *inArray) {
	kinParams.tx01 = inArray[0];
	kinParams.ty01 = inArray[1];
	kinParams.tz01 = inArray[2];
	kinParams.rz01 = inArray[3];
	kinParams.ry34 = inArray[4];
	kinParams.lCath = inArray[5];
}
FwdK6A::FwdK6A(KINEMATICPARAMS6A inParams) {
	kinParams = inParams;
}
Eigen::Matrix4d FwdK6A::qps2H01(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;
	R = Eigen::AngleAxisd(kinParams.rz01, Eigen::Vector3d::UnitZ())	* Eigen::AngleAxisd(qps[0], Eigen::Vector3d::UnitX());
	Rxyz.setIdentity();
	Rxyz.block<3, 3>(0, 0) = R;
	Rxyz(0, 3) = kinParams.tx01;
	Rxyz(1, 3) = kinParams.ty01;
	Rxyz(2, 3) = kinParams.tz01;

	Eigen::Matrix4d H01;
	H01.setIdentity();
	H01 = H01 * Rxyz;

	return H01;
}
Eigen::Matrix4d FwdK6A::qps2H02(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;

	R = Eigen::AngleAxisd(qps[1], Eigen::Vector3d::UnitZ());
	Rxyz.setIdentity();
	Rxyz.block<3, 3>(0, 0) = R;
	Rxyz(0, 3) = physicalDims.lProx; // translate from the coaxial input to the pitch axis

	return FwdK6A::qps2H01(qps) * Rxyz;
}
Eigen::Matrix4d FwdK6A::qps2H03(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;

	R = Eigen::AngleAxisd(qps[2], Eigen::Vector3d::UnitX());
	Rxyz.setIdentity();
	Rxyz(0, 3) = physicalDims.lPtch; // translate from the ptich axis to the roll face
	Rxyz.block<3, 3>(0, 0) = R;

	return FwdK6A::qps2H02(qps) * Rxyz;
}
Eigen::Matrix4d FwdK6A::qps2H04(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz, H;

	Rxyz.setIdentity();
	H.setIdentity();

	//translate from pitch to catheter base
	H(0, 3) = physicalDims.lRoll;
	//rotate to catheter plane
	R = Eigen::AngleAxisd( kinParams.ry34, Eigen::Vector3d::UnitY());
	H.block<3, 3>(0, 0) = R;
	//translate to the tip
	if (qps[3] < 1e-3) {
		qps[3] = 1e-3;
	}
	double r = kinParams.lCath / qps[3]; //radius of catheter arc
	Rxyz(0, 3) = r * sin(qps[3]);
	Rxyz(1, 3) = r * (1 - cos(qps[3]));
	H = H*Rxyz;
	//rotate to tip
	R = Eigen::AngleAxisd(qps[3], Eigen::Vector3d::UnitZ());
	Rxyz.setIdentity();
	Rxyz.block<3, 3>(0, 0) = R;
	H = H*Rxyz;

	return FwdK6A::qps2H03(qps) * H;
}
Eigen::Matrix4d FwdK6A::qps2H05(double *qps) {
	Eigen::Matrix4d Rxyz;

	Rxyz.setIdentity();
	Rxyz(0, 3) = qps[4]; //Tx

	return FwdK6A::qps2H04(qps) * Rxyz;
}

