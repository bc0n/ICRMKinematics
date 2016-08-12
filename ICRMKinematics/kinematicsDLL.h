#pragma once

#ifndef CICRM
#define CICRM

// building a DLL
#define DLLIMPORT __declspec (dllexport)

#ifdef __cplusplus
extern "C" { // using a C++ compiler
#endif

	typedef struct kinematics kinematics; //make class opaque to the wrapper---why do we want this?
	//DLLIMPORT kinematics* createKinematics(void); apparently not needed

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

	DLLIMPORT int get11AH01(double *qps, double *kinArray, double *arrayH01);
	DLLIMPORT int get11AH02(double *qps, double *kinArray, double *arrayH02);
	DLLIMPORT int get11AH03(double *qps, double *kinArray, double *arrayH03);
	DLLIMPORT int get11AH04(double *qps, double *kinArray, double *arrayH04);
	DLLIMPORT int get11AH05(double *qps, double *kinArray, double *arrayH05);

	//taskDefinitions --- expanded to encompass each FwdK and Task template variant
	DLLIMPORT int getTask5A_xyz(double *qps, double *kinArray, double *xyz);
	DLLIMPORT int getTask6A_xyz(double *qps, double *kinArray, double *xyz);
	DLLIMPORT int getTask11A_xyz(double *qps, double *kinArray, double *xyz);
	DLLIMPORT int getTask6A_phiPsi(double *qps, double *kinArray, double *pp);
	DLLIMPORT int getTask11A_phiPsi(double *qps, double *kinArray, double *pp);
	DLLIMPORT int getTask5A_xyzuxuyuz(double *qps, double *kinArray, double *xyz, double *uxyz);
	DLLIMPORT int getTask6A_xyzuxuyuz(double *qps, double *kinArray, double *xyz, double *uxyz);
	DLLIMPORT int getTask11A_xyzuxuyuz(double *qps, double *kinArray, double *xyz, double *uxyz);
	
	//funs
	DLLIMPORT int fun_qp0_xyzuxuyuz5A( int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin);
	DLLIMPORT int fun_qp0_xyzuxuyuz6A( int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin);
	DLLIMPORT int fun_qp0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin);
	DLLIMPORT int fun_kn0_xyzuxuyuz5A( int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kns0, double *fmin);
	DLLIMPORT int fun_kn0_xyzuxuyuz6A( int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kns0, double *fmin);
	DLLIMPORT int fun_kn0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kns0, double *fmin);
	DLLIMPORT int fun_qp0kn0_xyzuxuyuz5A( int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qps0, double *kns0, double *fmin);
	DLLIMPORT int fun_qp0kn0_xyzuxuyuz6A( int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qps0, double *kns0, double *fmin);
	DLLIMPORT int fun_qp0kn0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qps0, double *kns0, double *fmin);

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
	DLLIMPORT int estimate_qps_xyz5A( double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *fmin);
	DLLIMPORT int estimate_qps_xyz6A( double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *fmin);
	DLLIMPORT int estimate_qps_xyz11A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *fmin);
	DLLIMPORT int estimate_qps_xyzuxuyuz5A( double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *uxyz, double *fmin);
	DLLIMPORT int estimate_qps_xyzuxuyuz6A( double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *uxyz, double *fmin);
	DLLIMPORT int estimate_qps_xyzuxuyuz11A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz, double *uxyz, double *fmin);

	//initial joint angle estimation
	DLLIMPORT int estimate_qp0_xyzuxuyuz5A( int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *kn0, double *nlArray, double *fmin);
	DLLIMPORT int estimate_qp0_xyzuxuyuz6A( int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *kn0, double *nlArray, double *fmin);
	DLLIMPORT int estimate_qp0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *kn0, double *nlArray, double *fmin);

	//inverse paramameter estimation
	DLLIMPORT int estimate_kn0_xyzuxuyuz5A( int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *kup, double *kdn, double *nlArray, double *fmin);
	DLLIMPORT int estimate_kn0_xyzuxuyuz6A( int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *kup, double *kdn, double *nlArray, double *fmin);
	DLLIMPORT int estimate_kn0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *kup, double *kdn, double *nlArray, double *fmin);

	//simultaneous inverse parameter and initial joint estimation
	DLLIMPORT int estimate_qp0kn0_xyzuxuyuz5A( int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *k50, double *k5up, double *k5dn,  double *nlArray, double *fmin);
	DLLIMPORT int estimate_qp0kn0_xyzuxuyuz6A( int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *k60, double *k6up, double *k6dn, double *nlArray, double *fmin);
	DLLIMPORT int estimate_qp0kn0_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *q0Lims, double *k110, double *k11up, double *k11dn, double *nlArray, double *fmin);


#ifdef __cplusplus
}
#endif

#endif //CICRM