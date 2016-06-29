#pragma once

#ifndef CICRM
#define CICRM

#include "kinematics_structs.h"
#include "forwardKinematics.h"
#include "taskDefinitions.h"
#include "ik_nlopt.h"
#include "ip_nlopt.h"

// building a DLL
#define DLLIMPORT __declspec (dllexport)

#ifdef __cplusplus
extern "C" { // using a C++ compiler
#endif

	typedef struct kinematics kinematics; //make class opaque to the wrapper---why do we want this?
	DLLIMPORT kinematics* createKinematics(void);


	DLLIMPORT int get5AH01(double *qps, double *kinArray, double *arrayH01);
	DLLIMPORT int get5AH02(double *qps, double *kinArray, double *arrayH02);
	DLLIMPORT int get5AH03(double *qps, double *kinArray, double *arrayH03);
	DLLIMPORT int get5AH04(double *qps, double *kinArray, double *arrayH04);
	DLLIMPORT int get5AH05(double *qps, double *kinArray, double *arrayH05);

	DLLIMPORT int get6AH01(double *qps, double *kinArray, double *arrayH01);
	DLLIMPORT int get6AH02(double *qps, double *kinArray, double *arrayH02);
	DLLIMPORT int get6AH03(double *qps, double *kinArray, double *arrayH03);
	DLLIMPORT int get6AH04(double *qps, double *kinArray, double *arrayH04);
	DLLIMPORT int get6AH05(double *qps, double *kinArray, double *arrayH05);
	DLLIMPORT int get11AH05(double *qps, double *kinArray, double *arrayH05);

	//taskDefinitions --- expanded to encompass each FwdK and Task template variant
	DLLIMPORT int getTask5A_xyz(double *qps, double *kinArray, double *xyz);
	DLLIMPORT int getTask6A_xyz(double *qps, double *kinArray, double *xyz);
	DLLIMPORT int getTask11A_xyz(double *qps, double *kinArray, double *xyz);
	DLLIMPORT int getTask6A_phiPsi(double *qps, double *kinArray, double *pp);
	DLLIMPORT int getTask11A_phiPsi(double *qps, double *kinArray, double *pp);
	DLLIMPORT int getTask6A_xyzuxuyuz(double *qps, double *kinArray, double *xyz, double *uxyz);
	DLLIMPORT int getTask11A_xyzuxuyuz(double *qps, double *kinArray, double *xyz, double *uxyz);
	

	// optimizers
	// search begins from kns|qps|xyz, result returned there
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

	//inverse kinematics
	DLLIMPORT int getQps_IKnlopt_xyz5A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *fmin);
	DLLIMPORT int getQps_IKnlopt_xyz6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *fmin);
	DLLIMPORT int getQps_IKnlopt_xyz11A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *fmin);
	DLLIMPORT int getQps_IKnlopt_xyzuxuyuz6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *uxyz, double *fmin);
	DLLIMPORT int getQps_IKnlopt_xyzuxuyuz11A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *uxyz, double *fmin);

	//initial joint angle estimation
	DLLIMPORT int funIP_qp0_xyz5A(int nSamps, double *stackedQ, double *stackedX, double *kn5a, double *qp0, double *fmin);
	DLLIMPORT int estimate_qp0_xyz5A(int nSamps, double *stackedQ, double *stackedX, double *kn5a, double *qp0, double *qpup, double *qpdn, double *nlArray, double *fmin);
	DLLIMPORT int funIP_qp0_xyzuxuyuz5A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn5a, double *qp0, double *fmin);
	DLLIMPORT int estimate_qp0_xyzuxuyuz5A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn5a, double *qp0, double *qpup, double *qpdn, double *nlArray, double *fmin);

	//inverse paramameter estimation
	DLLIMPORT int estimate_kn0_xyz5A(int nSamps, double *stackedQ, double *stackedX, double *kn0, double *k5up, double *k5dn, double *nlArray, double *fmin);

	DLLIMPORT int fun_kn0_xyz5A(int nSamps, double *stackedQ, double *stackedX, double *kn0, double *fmin);
	DLLIMPORT int fun_qp0kn0_xyzdotu11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *pms0, double *fmin);
	DLLIMPORT int fun_qp0kn0_xyzpp11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *pms0, double *fmin);

	//simultaneous inverse parameter and initial joint estimation
	DLLIMPORT int estimate_qp0kn0_xyz5A(int nSamps, double *stackedQ, double *stackedX, double *kn0, double *k5up, double *k5dn, double *qp0, double *qpupdn, double *nlArray, double *fmin);
	DLLIMPORT int estimate_qp0kn0_xyzdotu11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *k110, double *k11up, double *k11dn, double *qp0, double *q0Lims, double *nlArray, double *fmin);
	DLLIMPORT int estimate_qp0kn0_xyzpp11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *k110, double *k11up, double *k11dn, double *qp0, double *q0Lims, double *nlArray, double *fmin);
	DLLIMPORT int estimate_qp0kn0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *k110, double *k11up, double *k11dn, double *qp0, double *q0Lims, double *nlArray, double *fmin);


#ifdef __cplusplus
}
#endif

KINEMATICPARAMS5A kinArray2Struct5A(double *kinArray);
KINEMATICPARAMS6A kinArray2Struct6A(double *kinArray);
KINEMATICPARAMS11A kinArray2Struct11A(double *kinArray);

NLOPTPARAMS nlArray2Struct(double *nlArray);
JOINTLIMITS jntArray2Struct(double *jntArray);

#endif //CICRM