#include "ip_nlopt.h"
#include <iostream>


template <class KPMS, class TFK>
InvP_nlopt<KPMS, TFK>::InvP_nlopt() {
	KPMS kup;
	KPMS kdn;
	JOINTLIMITS q0Lims;
	NLOPTPARAMS nlParams;
}
template <class KPMS, class TFK>
InvP_nlopt<KPMS, TFK>::InvP_nlopt(KPMS kupArg, KPMS kdnArg, JOINTLIMITS q0LimsArg, NLOPTPARAMS nlArg) {
	kup = kupArg;
	kdn = kdnArg;
	q0Lims = q0LimsArg;
	nlParams = nlArg;
}

// qp0
template <class KPMS, class TFK>
double err_Qp_xyzuxuyuz(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData;
	TFK *pfk;
	pfk = (TFK*)(pfip->fk);

	double qps[5];
	double accum = 0.0;
	
	//calculate error over all samps
	Eigen::Matrix4d Hnew;
	Eigen::Vector3d uxyz;
	Eigen::Vector3d unew;
	Eigen::Vector3d xyz;
	Eigen::VectorXd evec(6);
	for (int isamp = 0; isamp < pfip->nSamps; isamp++) {
		//update q
		for (int j = 0; j < 5; j++) {
			qps[j] = pfip->stackedQ[isamp * 5 + j] + x[j];
		}

		//new position
		Hnew = (*pfk).qps2H05(qps);

		//error
		evec(0) = pfip->stackedX[isamp * 3 + 0] - Hnew(0, 3);
		evec(1) = pfip->stackedX[isamp * 3 + 1] - Hnew(1, 3);
		evec(2) = pfip->stackedX[isamp * 3 + 2] - Hnew(2, 3);

		evec(3) = pfip->stackedU[isamp * 3 + 0] - Hnew(0, 0);
		evec(4) = pfip->stackedU[isamp * 3 + 1] - Hnew(1, 0);
		evec(5) = pfip->stackedU[isamp * 3 + 2] - Hnew(2, 0);

		accum += evec.norm();
	}
	return accum;
}

template <class KPMS, class TFK>
void InvP_nlopt<KPMS, TFK>::funQp(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qp0, double *kn0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;
	KPMS kp;
	kinematicArray2Struct(kn0, kp);
	TFK fk(kp);
	fip.fk = &fk;
	
	std::vector<double> x(kup.nParams);
	for (int i = 0; i < kup.nParams; i++) { x[i] = kn0[i]; }

	*fmin = err_Qp_xyzuxuyuz<KPMS, TFK>(x, x, &fip);
}

template <class KPMS, class TFK>
int InvP_nlopt<KPMS, TFK>::estimateQp(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;
	KPMS kp;
	kinematicArray2Struct(kn0, kp);
	TFK fk(kp);
	fip.fk = &fk;

	// get an alg object
	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), nQps); //there are 5 qps + 11 kinematic params in knParams11A

	// set objective
	void *pfip;
	pfip = &fip;
	alg.set_min_objective(err_Qp_xyzuxuyuz<KPMS, TFK>, pfip);

	//FIPDATA *ppfip;
	//ppfip = (FIPDATA*)pfip;
	////FwdK5A *pfk5;
	////pfk5 = (FwdK5A*)(ppfip->fk);
	//TFK *pfk;
	//pfk = (TFK*)(ppfip->fk);
	//std::cout << (*pfk).nParams << std::endl;

	// set initial step size (only used by nongradient methods)
	std::vector<double> iniStep(nQps);
	for (int i = 0; i < nQps; i++) { iniStep[i] = .001; }
	alg.set_initial_step(iniStep);
	alg.set_default_initial_step(iniStep);

	//set boundary constraints based on qpLast and stdev
	std::vector<double> limup(nQps);
	std::vector<double> limdn(nQps);
	for (int i = 0; i < nQps; i++) {
		limup[i] = q0Lims.up[i];
		limdn[i] = q0Lims.dn[i];
	}

	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to start within bounds
	std::vector<double> x(nQps);
	for (int i = 0; i < nQps; i++) { x[i] = qp0[i]; }
	for (int i = 0; i < nQps; i++) {
		if (x[i] <= limdn[i]) {
			x[i] = limdn[i];
		}
		if (limup[i] <= x[i]) {
			x[i] = limup[i];
		}
	}

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	//may need to set the local optimizer, eg MLSL http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
	if (ipTranslateNLOptAlg(nlParams.method) == GN_MLSL || ipTranslateNLOptAlg(nlParams.method) == GN_MLSL_LDS) {
		nlopt::opt lAlg(nlopt::LN_NELDERMEAD, alg.get_dimension());

		lAlg.set_min_objective(err_Qp_xyzuxuyuz<KPMS, TFK>, pfip);
		lAlg.set_upper_bounds(limup);
		lAlg.set_lower_bounds(limdn);
		lAlg.set_stopval(nlParams.minFunVal);
		lAlg.set_ftol_abs(nlParams.tolFunAbs);
		lAlg.set_xtol_abs(nlParams.tolXAbs);
		lAlg.set_maxeval(nlParams.maxIts);
		lAlg.set_maxtime(10.0);//lAlg.set_maxtime(nlParams.maxTimeSec);

		alg.set_local_optimizer(lAlg);
	}

	// solve
	nlopt::result res;
	try {
		res = alg.optimize(x, *fmin);
	}
	catch (nlopt::roundoff_limited e) {
		res = nlopt::FAILURE;
		std::cout << "Caught RoundoffLimited: " << e.what() << std::endl;
	}
	catch (nlopt::forced_stop e) {
		res = nlopt::FAILURE;
		std::cout << "Caught ForcedStop: " << e.what() << std::endl;
	}
	catch (std::runtime_error e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	catch (std::invalid_argument e) {
		res = nlopt::FAILURE;
		printf("Caught: %s for x=[ %4.4f %4.4f %4.4f %4.4f %4.4f...]\n", e.what(), x[0], x[1], x[2], x[3], x[4]);
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}

	// unpack from solver
	for (int i = 0; i < nQps; i++) { qp0[i] = x[i]; }
	*fmin /= nSamps;

	return res;
}


// kn0
template <class KPMS, class TFK>
double err_Kn_xyzuxuyuz(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData;

	double qps[5];
	double accum = 0.0;

	//        0    1    2    3    4     5     6    7     8     9   10
	//x = [tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45
	KPMS kp;
	double kna[kp.nParams];
	for (int i = 0; i < kp.nParams; i++) {
		kna[i] = x[i]; //kn params follow the q0s
	}
	kinematicArray2Struct(kna, kp);

	// FK with the new params
	TFK fk(kp);

	//calculate error over all samps
	Eigen::Matrix4d Hnew;
	Eigen::Vector3d uxyz;
	Eigen::Vector3d unew;
	Eigen::Vector3d xyz;
	Eigen::VectorXd evec(6);
	for (int isamp = 0; isamp < pfip->nSamps; isamp++) {
		//update q
		for (int j = 0; j < 5; j++) {
			qps[j] = pfip->stackedQ[isamp * 5 + j];
		}

		//new position
		Hnew = fk.qps2H05(qps);

		//error
		evec(0) = pfip->stackedX[isamp * 3 + 0] - Hnew(0, 3);
		evec(1) = pfip->stackedX[isamp * 3 + 1] - Hnew(1, 3);
		evec(2) = pfip->stackedX[isamp * 3 + 2] - Hnew(2, 3);

		evec(3) = pfip->stackedU[isamp * 3 + 0] - Hnew(0, 0);
		evec(4) = pfip->stackedU[isamp * 3 + 1] - Hnew(1, 0);
		evec(5) = pfip->stackedU[isamp * 3 + 2] - Hnew(2, 0);

		accum += evec.norm();
	}
	return accum;
}

template <class KPMS, class TFK>
void InvP_nlopt<KPMS, TFK>::funKn(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *kn0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	std::vector<double> x(kup.nParams);
	for (int i = 0; i < kup.nParams; i++) { x[i] = kn0[i]; }

	*fmin = err_QpKn_xyzuxuyuz<KPMS, TFK>(x, x, &fip);
}

template <class KPMS, class TFK>
int InvP_nlopt<KPMS, TFK>::estimateKn(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	// get an alg object
	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), kup.nParams); //there are 5 qps + 11 kinematic params in knParams11A
	
	// set objective
	void *pfip;
	pfip = &fip;
	alg.set_min_objective(err_Kn_xyzuxuyuz<KPMS, TFK>, pfip);

	// set initial step size (only used by nongradient methods)
	std::vector<double> iniStep(kup.nParams);
	for (int i = 0; i < kup.nParams; i++) { iniStep[i] = .001; }
	alg.set_initial_step(iniStep);
	alg.set_default_initial_step(iniStep);

	//set boundary constraints based on qpLast and stdev
	//      0  1  2  3  4     5    6    7    8    9    10
	//x = [tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45
	std::vector<double> limup(kup.nParams);
	std::vector<double> limdn(kup.nParams);
	double knaup[kup.nParams], knadn[kup.nParams];
	kinematicStruct2Array(kup, knaup);
	kinematicStruct2Array(kdn, knadn);
	for (int i = 0; i < kup.nParams; i++) {
		limup[i] = knaup[i];
		limdn[i] = knadn[i];
	}

	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to start within bounds
	std::vector<double> x(kup.nParams);
	for (int i = 0; i < kup.nParams; i++) { x[i] = kn0[i]; }
	for (int i = 0; i < kup.nParams; i++) {
		if (x[i] <= limdn[i]) {
			x[i] = limdn[i];
		}
		if (limup[i] <= x[i]) {
			x[i] = limup[i];
		}
	}

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	//may need to set the local optimizer, eg MLSL http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
	if (ipTranslateNLOptAlg(nlParams.method) == GN_MLSL || ipTranslateNLOptAlg(nlParams.method) == GN_MLSL_LDS) {
		nlopt::opt lAlg(nlopt::LN_NELDERMEAD, alg.get_dimension());

		lAlg.set_min_objective(err_Kn_xyzuxuyuz<KPMS, TFK>, pfip);
		lAlg.set_upper_bounds(limup);
		lAlg.set_lower_bounds(limdn); 
		lAlg.set_stopval(nlParams.minFunVal);
		lAlg.set_ftol_abs(nlParams.tolFunAbs);
		lAlg.set_xtol_abs(nlParams.tolXAbs);
		lAlg.set_maxeval(nlParams.maxIts);
		lAlg.set_maxtime(10.0);//lAlg.set_maxtime(nlParams.maxTimeSec);

		alg.set_local_optimizer(lAlg);
	}

	// solve
	nlopt::result res;
	try {
		res = alg.optimize(x, *fmin);
	}
	catch (nlopt::roundoff_limited e) {
		res = nlopt::FAILURE;
		std::cout << "Caught RoundoffLimited: " << e.what() << std::endl;
	}
	catch (nlopt::forced_stop e) {
		res = nlopt::FAILURE;
		std::cout << "Caught ForcedStop: " << e.what() << std::endl;
	}
	catch (std::runtime_error e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	catch (std::invalid_argument e) {
		res = nlopt::FAILURE;
		printf("Caught: %s for x=[ %4.4f %4.4f %4.4f %4.4f %4.4f...]\n", e.what(), x[0], x[1], x[2], x[3], x[4]);
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}

	// unpack from solver
	for (int i = 0; i < kup.nParams; i++) { kn0[i] = x[i]; }
	*fmin /= nSamps;

	return res;
}


//simultaneous qp0 & kn0
template <class KPMS, class TFK>
double err_QpKn_xyzuxuyuz(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData;

	double qps[5];
	double accum = 0.0;

	//      0  1  2  3  4     5    6    7    8    9    10    11   12    13    14   15
	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45
	KPMS kp;
	double kna[kp.nParams];
	for (int i = 0; i < kp.nParams; i++) {
		kna[i] = x[5 + i]; //kn params follow the q0s
	}
	kinematicArray2Struct(kna, kp);
	
	// FK with the new params
	TFK fk(kp);

	//calculate error over all samps
	Eigen::Matrix4d Hnew;
	Eigen::Vector3d uxyz;
	Eigen::Vector3d unew;
	Eigen::Vector3d xyz;
	Eigen::VectorXd evec(6);
	for (int isamp = 0; isamp < pfip->nSamps; isamp++) {
		//update q
		for (int j = 0; j < 5; j++) {
			qps[j] = x[j] + pfip->stackedQ[isamp * 5 + j];
		}

		//new position
		Hnew = fk.qps2H05(qps);

		//error
		evec(0) = pfip->stackedX[isamp * 3 + 0] - Hnew(0, 3);
		evec(1) = pfip->stackedX[isamp * 3 + 1] - Hnew(1, 3);
		evec(2) = pfip->stackedX[isamp * 3 + 2] - Hnew(2, 3);

		evec(3) = pfip->stackedU[isamp * 3 + 0] - Hnew(0, 0);
		evec(4) = pfip->stackedU[isamp * 3 + 1] - Hnew(1, 0);
		evec(5) = pfip->stackedU[isamp * 3 + 2] - Hnew(2, 0);

		accum += evec.norm();
	}
	return accum;
}

template <class KPMS, class TFK>
void InvP_nlopt<KPMS, TFK>::funQpKn(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qp0, double *kn0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	std::vector<double> x(nQps + kup.nParams);
	for (int i = 0; i < nQps; i++) { x[i] = qp0[i]; }
	for (int i = 0; i < kup.nParams; i++) { x[i + nQps] = kn0[i]; }

	*fmin = err_QpKn_xyzuxuyuz<KPMS,TFK>(x, x, &fip);
}

template <class KPMS, class TFK>
int InvP_nlopt<KPMS, TFK>::estimateQpKn(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;
	
	// get an alg object
	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), nQps + kup.nParams); //there are 5 qps + n( kinematic params in knParams11A )

	// set objective
	void *pfip;
	pfip = &fip;
	alg.set_min_objective(err_QpKn_xyzuxuyuz<KPMS,TFK>, pfip);

	// set initial step size (only used by nongradient methods)
	std::vector<double> iniStep(nQps + kup.nParams);
	for (int i = 0; i < nQps + kup.nParams; i++) { iniStep[i] = .001; }
	alg.set_initial_step(iniStep);
	alg.set_default_initial_step(iniStep);

	//set boundary constraints based on qpLast and stdev
	//      0  1  2  3  4     5    6    7    8    9    10    11   12    13    14   15
	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45
	std::vector<double> limup(nQps + kup.nParams);
	std::vector<double> limdn(nQps + kup.nParams);
	for (int i = 0; i < nQps; i++) {
		limup[i] = q0Lims.up[i];
		limdn[i] = q0Lims.dn[i];
	}
	double knaup[kup.nParams], knadn[kup.nParams];
	kinematicStruct2Array(kup, knaup);
	kinematicStruct2Array(kdn, knadn);
	for (int i = 0; i < kup.nParams; i++) {
		limup[nQps + i] = knaup[i];
		limdn[nQps + i] = knadn[i];
	}

	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to be within bounds
	std::vector<double> x(nQps + kup.nParams);
	for (int i = 0; i < nQps; i++) { x[i] = qp0[i]; }
	for (int i = 0; i < kup.nParams; i++) { x[nQps + i] = kn0[i]; }
	for (int i = 0; i < nQps + kup.nParams; i++) {
		if (x[i] <= limdn[i]) {
			x[i] = limdn[i];
		}
		if (limup[i] <= x[i]) {
			x[i] = limup[i];
		}
	}

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	//may need to set the local optimizer, eg MLSL http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
	if (ipTranslateNLOptAlg(nlParams.method) == GN_MLSL || ipTranslateNLOptAlg(nlParams.method) == GN_MLSL_LDS) {
		nlopt::opt lAlg(nlopt::LN_NELDERMEAD, alg.get_dimension());

		lAlg.set_min_objective(err_QpKn_xyzuxuyuz<KPMS, TFK>, pfip);
		lAlg.set_upper_bounds(limup);
		lAlg.set_lower_bounds(limdn);
		lAlg.set_stopval(nlParams.minFunVal);
		lAlg.set_ftol_abs(nlParams.tolFunAbs);
		lAlg.set_xtol_abs(nlParams.tolXAbs);
		lAlg.set_maxeval(nlParams.maxIts);
		lAlg.set_maxtime(10.0);//lAlg.set_maxtime(nlParams.maxTimeSec);

		alg.set_local_optimizer(lAlg);
	}

	// solve
	nlopt::result res;
	try {
		res = alg.optimize(x, *fmin);
	}
	catch (nlopt::roundoff_limited e) {
		res = nlopt::FAILURE;
		std::cout << "Caught RoundoffLimited: " << e.what() << std::endl;
	}
	catch (nlopt::forced_stop e) {
		res = nlopt::FAILURE;
		std::cout << "Caught ForcedStop: " << e.what() << std::endl;
	}
	catch (std::runtime_error e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	catch (std::invalid_argument e) {
		res = nlopt::FAILURE;
		printf("Caught: %s for x=[", e.what());
		for (int i = 0; i < nQps + kup.nParams; i++) { printf("%f ", x[i]); }
		printf("]\n");
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	
	// unpack from solver
	for (int i = 0; i < nQps; i++) { qp0[i] = x[i]; }
	for (int i = 0; i < kup.nParams; i++) { kn0[i] = x[nQps + i]; }
	*fmin /= nSamps;

	return res;
}



//Struct to Array
void kinematicStruct2Array(KINEMATICPARAMS5A kns, double *kna) {
	kna[0] = kns.tx01;
	kna[1] = kns.ty01;
	kna[2] = kns.tz01;
	kna[3] = kns.rz01;
	kna[4] = kns.lCath;
}
void kinematicStruct2Array(KINEMATICPARAMS6A kns, double *kna) {
	kna[0] = kns.tx01;
	kna[1] = kns.ty01;
	kna[2] = kns.tz01;
	kna[3] = kns.rz01;
	kna[4] = kns.ry34;
	kna[5] = kns.lCath;
}
void kinematicStruct2Array(KINEMATICPARAMS11A kns, double *kna) {
	kna[0] = kns.tx01;
	kna[1] = kns.ty01;
	kna[2] = kns.tz01;
	kna[3] = kns.ry01;
	kna[4] = kns.rz01;
	kna[5] = kns.ry34;
	kna[6] = kns.rz34;
	kna[7] = kns.kAlpha;
	kna[8] = kns.eAlpha;
	kna[9] = kns.lCath;
	kna[10] = kns.ry45;
}

//Array to Struct
void kinematicArray2Struct(double *kna, KINEMATICPARAMS5A kns) {
	kns.tx01 = kna[0];
	kns.ty01 = kna[1];
	kns.tz01 = kna[2];
	kns.rz01 = kna[3];
	kns.lCath = kna[4];
}
void kinematicArray2Struct(double *kna, KINEMATICPARAMS6A kns) {
	kns.tx01 = kna[0];
	kns.ty01 = kna[1];
	kns.tz01 = kna[2];
	kns.rz01 = kna[3];
	kns.ry34 = kna[4];
	kns.lCath = kna[5];
}
void kinematicArray2Struct(double *kna, KINEMATICPARAMS11A kns) {
	kns.tx01 = kna[0];
	kns.ty01 = kna[1];
	kns.tz01 = kna[2];
	kns.ry01 = kna[3];
	kns.rz01 = kna[4];
	kns.ry34 = kna[5];
	kns.rz34 = kna[6];
	kns.kAlpha = kna[7];
	kns.eAlpha = kna[8];
	kns.lCath = kna[9];
	kns.ry45 = kna[10];
}

nlopt::algorithm ipTranslateNLOptAlg(nlMethod method) {
	//std::cout << "recv " << method << std::endl;
	switch (method) {
	case nlMethod::GN_DIRECT: return nlopt::GN_DIRECT; //0
	case nlMethod::GN_DIRECT_L: return nlopt::GN_DIRECT_L;//1
	case nlMethod::GN_DIRECT_L_NOSCAL: return nlopt::GN_DIRECT_L_NOSCAL;//2
	case nlMethod::GN_DIRECT_L_RAND: return nlopt::GN_DIRECT_L_RAND;//3
	case nlMethod::GN_DIRECT_L_RAND_NOSCAL: return nlopt::GN_DIRECT_L_RAND_NOSCAL;//4
	case nlMethod::GN_DIRECT_NOSCAL: return nlopt::GN_DIRECT_NOSCAL;//5
	case nlMethod::GN_ESCH: return nlopt::GN_ESCH;//6
	case nlMethod::GN_ISRES: return nlopt::GN_ISRES;//7
	case nlMethod::GN_MLSL: return nlopt::GN_MLSL;//8
	case nlMethod::GN_MLSL_LDS: return nlopt::GN_MLSL_LDS;//9
	case nlMethod::GN_ORIG_DIRECT: return nlopt::GN_ORIG_DIRECT;//10
	case nlMethod::GN_ORIG_DIRECT_L: return nlopt::GN_ORIG_DIRECT_L;//11
	case nlMethod::LN_BOBYQA: return nlopt::LN_BOBYQA;//12
	case nlMethod::LN_COBYLA: return nlopt::LN_COBYLA;//13
	case nlMethod::LN_NELDERMEAD: return nlopt::LN_NELDERMEAD;//14
	case nlMethod::LN_NEWUOA: return nlopt::LN_NEWUOA;//15
	case nlMethod::LN_NEWUOA_BOUND: return nlopt::LN_NEWUOA_BOUND;//16
	case nlMethod::LN_PRAXIS: return nlopt::LN_PRAXIS;//17
	case nlMethod::LN_SBPLX: return nlopt::LN_SBPLX;//18
	default: return nlopt::LN_BOBYQA;
	}
}

//useless definition to help linker, one for every possible template value; absence leads to LNK2001s
template class InvP_nlopt<KINEMATICPARAMS5A, FwdK5A>;
template class InvP_nlopt<KINEMATICPARAMS6A, FwdK6A>;
template class InvP_nlopt<KINEMATICPARAMS11A, FwdK11A>;

