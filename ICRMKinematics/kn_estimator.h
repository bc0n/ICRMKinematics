#pragma once

#include "kinematics_structs.h"
#include "taskDefinitions.h"
#include "nlopt.hpp"

typedef struct tagFKNDATA {
	double *cmdQps;
	double *cmdQpss;
	double *meaXsYsZs;
	int nSets;
	Eigen::Vector3d meaVec;
} FKNDATA;



class KNEstimate {
private:
	NLOPTPARAMS nlParams;
public:
	KNEstimate(NLOPTPARAMS nParams);
	int estimate_AllSingleXYZ(double *cmdQps, double *meaXYZ, double *kinArray);
	//int estimateAll_costH05(double *cmdQps, double *meaXYZ, double *kinArray); ??

	int estimate_AllMultipleXYZ(int nSets, double *cmdQpss, double *meaXsYsZs, double *kinArray, double *kinLim, double *fmin);
};


nlopt::algorithm translateNLAlgForKN(nlMethod method);