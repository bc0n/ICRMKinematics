#include "forwardKinematics.h"
#include <iostream>

FwdK5A::FwdK5A() {
}
FwdK5A::FwdK5A(double *inArray) {
	kinParams.tx01 = inArray[0];
	kinParams.ty01 = inArray[1];
	kinParams.tz01 = inArray[2];
	kinParams.rz01 = inArray[3];
	kinParams.lCath = inArray[4];
}
FwdK5A::FwdK5A(KINEMATICPARAMS5A inParams) {
	kinParams = inParams;
}
Eigen::Matrix4d FwdK5A::qps2H01(double *qps) {
	//printf("kinParams = [%f %f %f %f %f]\n", kinParams.tx01, kinParams.ty01, kinParams.tz01, kinParams.rz01, kinParams.lCath);

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
Eigen::Matrix4d FwdK5A::qps2H02(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;
	
	R = Eigen::AngleAxisd(qps[1], Eigen::Vector3d::UnitZ());
	Rxyz.setIdentity();
	Rxyz.block<3, 3>(0, 0) = R;
	Rxyz(0, 3) = physicalDims.lProx; // translate from the coaxial input to the pitch axis

	return FwdK5A::qps2H01(qps) * Rxyz;
}
Eigen::Matrix4d FwdK5A::qps2H03(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;

	R = Eigen::AngleAxisd(qps[2], Eigen::Vector3d::UnitX());
	Rxyz.setIdentity();
	Rxyz(0, 3) = physicalDims.lPtch; // translate from the ptich axis to the roll face
	Rxyz.block<3, 3>(0, 0) = R;

	return FwdK5A::qps2H02(qps) * Rxyz;
}
Eigen::Matrix4d FwdK5A::qps2H04(double *qps) {
	Eigen::Matrix3d R;
	Eigen::Matrix4d Rxyz;

	if (qps[3] < 1e-3) {
		qps[3] = 1e-3;
	}

	double r = kinParams.lCath / qps[3]; //radius of catheter arc
	double ty = r * (1 - cos(qps[3]));
	double tx = r * sin(qps[3]);

	R = Eigen::AngleAxisd(qps[3], Eigen::Vector3d::UnitZ());
	Rxyz.setIdentity();
	Rxyz(0, 3) = tx + physicalDims.lRoll; // translate the to the catheter base, then x,y,rotate
	Rxyz(1, 3) = ty;
	Rxyz.block<3, 3>(0, 0) = R;

	return FwdK5A::qps2H03(qps) * Rxyz;
}
Eigen::Matrix4d FwdK5A::qps2H05(double *qps) {
	Eigen::Matrix4d Rxyz;
	Eigen::Matrix4d H05;

	Rxyz.setIdentity();
	Rxyz(0, 3) = qps[4]; // translate the virtual distance

	return FwdK5A::qps2H04(qps) * Rxyz;
}

