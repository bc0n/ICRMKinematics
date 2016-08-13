#include "forwardKinematics.h"
#include <iostream>

FwdK11A::FwdK11A() {
}
FwdK11A::FwdK11A(double *inArray) {
	kinParams.tx01 = inArray[0];
	kinParams.ty01 = inArray[1];
	kinParams.tz01 = inArray[2];
	kinParams.ry01 = inArray[3];
	kinParams.rz01 = inArray[4];
	kinParams.ry34 = inArray[5];
	kinParams.rz34 = inArray[6];
	kinParams.kAlpha = inArray[7];
	kinParams.eAlpha = inArray[8];
	kinParams.lCath = inArray[9];
	kinParams.ry45 = inArray[10];
}
FwdK11A::FwdK11A(KINEMATICPARAMS11A inParams) {
	kinParams = inParams;
}
Eigen::Matrix4d FwdK11A::qps2H01(double *qps) {
	Eigen::Matrix4d Rxyz;
	Rxyz.setIdentity();
	Rxyz(0, 3) = kinParams.tx01;
	Rxyz(1, 3) = kinParams.ty01;
	Rxyz(2, 3) = kinParams.tz01;
	Eigen::Matrix3d R;
	R = Eigen::AngleAxisd(kinParams.ry01, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(kinParams.rz01, Eigen::Vector3d::UnitZ())
		* Eigen::AngleAxisd(qps[0], Eigen::Vector3d::UnitX());
	Rxyz.block<3, 3>(0, 0) = R;

	Eigen::Matrix4d H01;
	H01.setIdentity();
	H01 = H01 * Rxyz;

	return H01;
}
Eigen::Matrix4d FwdK11A::qps2H02(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;

	R = Eigen::AngleAxisd(qps[1], Eigen::Vector3d::UnitZ());
	Rxyz.setIdentity();
	Rxyz.block<3, 3>(0, 0) = R;
	Rxyz(0, 3) = physicalDims.lProx; // translate from the coaxial input to the pitch axis

	return FwdK11A::qps2H01(qps) * Rxyz;
}
Eigen::Matrix4d FwdK11A::qps2H03(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;

	R = Eigen::AngleAxisd(qps[2], Eigen::Vector3d::UnitX());
	Rxyz.setIdentity();
	Rxyz(0, 3) = physicalDims.lPtch; // translate from the ptich axis to the roll face
	Rxyz.block<3, 3>(0, 0) = R;

	return FwdK11A::qps2H02(qps) * Rxyz;
}
Eigen::Matrix4d FwdK11A::qps2H04(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz, H;

	Rxyz.setIdentity();
	H.setIdentity();

	//translate from pitch to catheter base
	H(0, 3) = physicalDims.lRoll;
	//rotate to catheter plane
	R = Eigen::AngleAxisd(kinParams.ry34, Eigen::Vector3d::UnitY()) * Eigen::AngleAxisd(kinParams.rz34, Eigen::Vector3d::UnitZ());
	H.block<3, 3>(0, 0) = R;

	//translate to the tip
	double al = pow(qps[3] * kinParams.kAlpha, kinParams.eAlpha); //parabolic catheter
	if (al < 1e-3) {
		al = 1e-3;
	}
	double r = kinParams.lCath / al; //'radius' of catheter arc
	Rxyz(0, 3) = r * sin(al);
	Rxyz(1, 3) = r * (1 - cos(al));
	H = H*Rxyz;
	//rotate to tip
	R = Eigen::AngleAxisd(al, Eigen::Vector3d::UnitZ());
	Rxyz.setIdentity();
	Rxyz.block<3, 3>(0, 0) = R;
	H = H*Rxyz;

	return FwdK11A::qps2H03(qps) * H;
}
Eigen::Matrix4d FwdK11A::qps2H05(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;
	Eigen::Matrix4d H05;

	R = Eigen::AngleAxisd(kinParams.ry45, Eigen::Vector3d::UnitY());

	Rxyz.setIdentity();
	Rxyz.block<3, 3>(0, 0) = R;
	Rxyz(0, 3) = qps[4]; //Tx

	return FwdK11A::qps2H04(qps) * Rxyz;
}
