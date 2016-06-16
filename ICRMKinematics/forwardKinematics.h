#pragma once


#include "kinematics_structs.h"
#include <vector>

class FwdK6A {
public:
	KINEMATICPARAMS6A kinParams;

	FwdK6A(); //use default kinematic params

	FwdK6A(KINEMATICPARAMS6A inParams); //specified

	// forwardK methods
	Eigen::Matrix4d qps2H01(double *qps);
	Eigen::Matrix4d qps2H02(double *qps);
	Eigen::Matrix4d qps2H03(double *qps);
	Eigen::Matrix4d qps2H04(double *qps);
	Eigen::Matrix4d qps2H05(double *qps);
};

class FwdK11A {
public:
	KINEMATICPARAMS11A kinParams;

	FwdK11A(); // default params

	FwdK11A(KINEMATICPARAMS11A inParams); //specified

	// forwardK methods
	Eigen::Matrix4d qps2H01(double *qps); //adds Ry01
	Eigen::Matrix4d qps2H02(double *qps);
	Eigen::Matrix4d qps2H03(double *qps);
	Eigen::Matrix4d qps2H04(double *qps); //adds Ry34 and Rz34
	Eigen::Matrix4d qps2H05(double *qps); //adds Ry45 and Rz45
};
