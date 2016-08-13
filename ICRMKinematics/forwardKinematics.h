#pragma once


#include "kinematics_structs.h"
#include <vector>


class FwdK5A {
private:
	PHYSICALDIMENSIONS physicalDims;
public:
	KINEMATICPARAMS5A kinParams;
	int nParams = 5;

	FwdK5A(); //use default kinematic params
	FwdK5A(double *inArray); // inArray = [tx01,ty01,tz01,rz01, cathL]
	FwdK5A(KINEMATICPARAMS5A inParams); //specified

	// forwardK methods
	Eigen::Matrix4d qps2H01(double *qps);
	Eigen::Matrix4d qps2H02(double *qps);
	Eigen::Matrix4d qps2H03(double *qps);
	Eigen::Matrix4d qps2H04(double *qps);
	Eigen::Matrix4d qps2H05(double *qps);

};

class FwdK6A {
private:
	PHYSICALDIMENSIONS physicalDims;
public:
	KINEMATICPARAMS6A kinParams;
	int nParams = 6;

	FwdK6A(); //use default kinematic params
	FwdK6A(double *inArray); // inArray = [tx01,ty01,tz01,rz01, ry34, cathL]
	FwdK6A(KINEMATICPARAMS6A inParams);

	// forwardK methods
	Eigen::Matrix4d qps2H01(double *qps);
	Eigen::Matrix4d qps2H02(double *qps);
	Eigen::Matrix4d qps2H03(double *qps);
	Eigen::Matrix4d qps2H04(double *qps);
	Eigen::Matrix4d qps2H05(double *qps);
};

class FwdK11A {
private:
	PHYSICALDIMENSIONS physicalDims;
public:
	KINEMATICPARAMS11A kinParams;
	int nParams = 11;

	FwdK11A(); // default params
	FwdK11A(double *inArray); // inArray = [tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45]
	FwdK11A(KINEMATICPARAMS11A inParams); 

	// forwardK methods
	Eigen::Matrix4d qps2H01(double *qps); //adds Ry01
	Eigen::Matrix4d qps2H02(double *qps);
	Eigen::Matrix4d qps2H03(double *qps);
	Eigen::Matrix4d qps2H04(double *qps); //adds Ry34 and Rz34
	Eigen::Matrix4d qps2H05(double *qps); //adds Ry45 and Rz45
};
