#pragma once

#ifndef CICRM
#define CICRM

#include "kinematics_structs.h"
#include "forwardKinematics.h"
#include "taskDefinitions.h"
#include "ik_newton.h"
#include "ik_nlopt.h"
#include "ip_nlopt.h"
#include "kn_estimator.h"

// building a DLL
#define DLLIMPORT __declspec (dllexport)

#ifdef __cplusplus
extern "C" { // using a C++ compiler
#endif

	typedef struct kinematics kinematics; //make class opaque to the wrapper---why do we want this?

	DLLIMPORT kinematics* createKinematics(void);

	DLLIMPORT int get6AH01(double *qps, double *kinArray, double *arrayH01);
	DLLIMPORT int get6AH02(double *qps, double *kinArray, double *arrayH02);
	DLLIMPORT int get6AH03(double *qps, double *kinArray, double *arrayH03);
	DLLIMPORT int get6AH04(double *qps, double *kinArray, double *arrayH04);
	DLLIMPORT int get6AH05(double *qps, double *kinArray, double *arrayH05);
	DLLIMPORT int get11AH05(double *qps, double *kinArray, double *arrayH05);

	//taskDefinitions
	DLLIMPORT int get6APoint(double *qps, double *kinArray, double *xyz);
	DLLIMPORT int get11APoint(double *qps, double *kinArray, double *xyz);
	DLLIMPORT int get6APhiPsi(double *qps, double *kinArray, double *pp);
	DLLIMPORT int get11APhiPsi(double *qps, double *kinArray, double *pp);
	DLLIMPORT int getCenterSum(double *qps, double *jntArray, double *cs);

	//inverse kinematic solvers
	DLLIMPORT int getQps_IKNewtonSinglePoint(double *qps, double *kinArray, double *nrArray, double *jntArray, double *xyz);
	DLLIMPORT int getQps_IKNewtonPriorityPointPhiPsi(double *qps, double *kinArray, double *nrArray, double *jntArray, double *xyz, double *phipsi);
	DLLIMPORT int getQps_IKNewtonPriorityPointCenterSum(double *qps, double *kinArray, double *nrArray, double *jntArray, double *xyz);
	DLLIMPORT int getQps_IKNewtonPriorityPointCenterSumStops(double *qps, double *kinArray, double *nrArray, double *jntArray, double *xyz);
	DLLIMPORT int getQps_IKNewtonAugmentedPointCenterSumStops(double *qps, double *kinArray, double *nrArray, double *jntArray, double *xyz);
	DLLIMPORT int getQps_IKNewtonAugmentedPointPhiPsiStops(double *qps, double *kinArray, double *nrArray, double *jntArray, double *xyz, double *phipsi);

	// search begins from qps|xyz, result returned there
	// returns -1 = failure
	//         -2 = invalid args
	//         -3 = out of memory
	//         -4 = roundoff limited
	//         -5 = forced stop
	//          1 = success
	//          2 = fun stopval reached
	//          3 = fun tol reached
	//          4 = x tol reached
	//          5 = max evals
	//          6 = max time reached
	DLLIMPORT int getQps_IKNLOpt_xyz6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz);
	DLLIMPORT int getQps_IKNLOpt_xyzpp6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *phipsi);
	DLLIMPORT int getQps_IKNLOpt_xyzpp11A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *phipsi);
	DLLIMPORT int getFunVal_xyzpp11A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *phipsi, double *fval);
	DLLIMPORT int getQps_IKNLOpt_xyzCtrSum6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz);
	DLLIMPORT int getQps_IKNLOpt_xyzppCtrSum6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *pp);
	DLLIMPORT int getQps_IKNLOpt_xyzMinChange6A(double *qps, double *xyz, double *qpsLast, double *kinArray, double *nlArray, double *jntArray);
	DLLIMPORT int getQps_IKNLOpt_xyzKalX6A(double *qps, double *xyz, double *qpsLast, double *qdsLast, double *kinArray, double *nlArray, double *jntArray, double *K11s, double *K21s, double tsSec);
	DLLIMPORT int getQps_IKNLOpt_xyzKalBounds6A(double *qps, double *xyz, double *qpsLast, double *qdsLast, double *kinArray, double *nlArray, double *jntArray, double *K11s, double *K21s, double *stdDevs, double tsSec);
	
	//inverse paramameter estimators
	DLLIMPORT int funIP_xyzdotu11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *pms0, double *fmin);
	//assumes qps0 = 0 and kpms0 are centered between up&dn
	DLLIMPORT int estimatePmsQ_IPNLOpt_xyzdotu11A_assumeX0(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *k11up, double *k11dn, double *q0Lims, double *nlArray, double *fmin);
	//uses the given qps0 and kpms0
	DLLIMPORT int estimatePmsQ_IPNLOpt_xyzdotu11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *k11up, double *k11dn, double *q0Lims, double *nlArray, double *qps0, double *kps0, double *fmin);

	DLLIMPORT int funIP_xyzpp11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *pms0, double *fmin);
	DLLIMPORT int estimatePmsQ_IPNLOpt_xyzpp11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *k11up, double *k11dn, double *q0Lims, double *nlArray, double *qps0, double *kps0, double *fmin);

#ifdef __cplusplus
}
#endif

KINEMATICPARAMS6A kinArray2Struct6A(double *kinArray);
KINEMATICPARAMS11A kinArray2Struct11A(double *kinArray);
NEWTONPARAMS nrArray2Struct(double *nrArray);
NLOPTPARAMS nlArray2Struct(double *nlArray);
JOINTLIMITS jntArray2Struct(double *jntArray);

#endif //CICRM