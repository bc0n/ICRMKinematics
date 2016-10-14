#include "ip_nlopt.h"
#include <iostream>


template <class KPMS, class TFK>
InvP_nlopt<KPMS, TFK>::InvP_nlopt() {
	KPMS kup;
	KPMS kdn;
	JOINTLIMITS q0Lims;
	NLOPTPARAMS nlParams;

	static bool knSub[kup.nParams]; //must be static
	for (int i = 0; i < kup.nParams; i++) {
		knSub[i] = true;
	}
	knSubset = knSub;

	double knDft[kup.nParams];
	kinematicStruct2Array(&kup, knDft);
	knDefault = knDft;
}
template <class KPMS, class TFK>
InvP_nlopt<KPMS, TFK>::InvP_nlopt(KPMS kupArg, KPMS kdnArg, JOINTLIMITS q0LimsArg, NLOPTPARAMS nlArg) {
	kup = kupArg;
	kdn = kdnArg;
	q0Lims = q0LimsArg;
	nlParams = nlArg;

	static bool knSub[kup.nParams];
	for (int i = 0; i < kup.nParams; i++) {
		knSub[i] = true;
	}
	knSubset = knSub;

	double knDft[kup.nParams];
	kinematicStruct2Array(&kup, knDft);
	knDefault = knDft;
}
template <class KPMS, class TFK>
InvP_nlopt<KPMS, TFK>::InvP_nlopt(KPMS kupArg, KPMS kdnArg, JOINTLIMITS q0LimsArg, NLOPTPARAMS nlArg, NLOPTPARAMS lnlArg) {
	kup = kupArg;
	kdn = kdnArg;
	q0Lims = q0LimsArg;
	nlParams = nlArg;
	localNLParams = lnlArg;

	static bool knSub[kup.nParams];
	for (int i = 0; i < kup.nParams; i++) {
		knSub[i] = true;
	}
	knSubset = knSub;

	double knDft[kup.nParams];
	kinematicStruct2Array(&kup, knDft);
	knDefault = knDft;
}
template <class KPMS, class TFK>
InvP_nlopt<KPMS, TFK>::InvP_nlopt(KPMS kupArg, KPMS kdnArg, bool *knSubsetArg, JOINTLIMITS q0LimsArg, NLOPTPARAMS nlArg) {
	kup = kupArg;
	kdn = kdnArg;
	q0Lims = q0LimsArg;
	nlParams = nlArg;
	knSubset = knSubsetArg;

	double knDft[kup.nParams];
	kinematicStruct2Array(&kup, knDft);//the non-estimated have kup=kdn
	knDefault = knDft;	
}
template <class KPMS, class TFK>
InvP_nlopt<KPMS, TFK>::InvP_nlopt(KPMS kupArg, KPMS kdnArg, bool *knSubsetArg, JOINTLIMITS q0LimsArg, NLOPTPARAMS nlArg, NLOPTPARAMS localNLArg) {
	kup = kupArg;
	kdn = kdnArg;
	q0Lims = q0LimsArg;
	nlParams = nlArg;
	localNLParams = localNLArg;
	knSubset = knSubsetArg;

	double knDft[kup.nParams];
	kinematicStruct2Array(&kup, knDft);//the non-estimated have kup=kdn
	knDefault = knDft;
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

	TFK fk(kn0);
	fip.fk = &fk;
	
	std::vector<double> x(kup.nParams);
	for (int i = 0; i < kup.nParams; i++) { x[i] = kn0[i]; }

	*fmin = err_Qp_xyzuxuyuz<KPMS, TFK>(x, x, &fip) / nSamps;
}

template <class KPMS, class TFK>
int InvP_nlopt<KPMS, TFK>::estimateQp(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;
	KPMS kp;
	kinematicArray2Struct(kn0, &kp);
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
	if (alg.get_algorithm() == nlopt::GN_MLSL_LDS) {
		nlopt::opt lAlg(ipTranslateNLOptAlg(localNLParams.method), alg.get_dimension());
std::cout << "Local optimization " << lAlg.get_algorithm_name() << std::endl;
		lAlg.set_min_objective(err_QpKn_xyzuxuyuz<KPMS, TFK>, pfip);
		lAlg.set_upper_bounds(limup);
		lAlg.set_lower_bounds(limdn);
		lAlg.set_stopval(localNLParams.minFunVal);
		lAlg.set_ftol_abs(localNLParams.tolFunAbs);
		lAlg.set_xtol_abs(localNLParams.tolXAbs);
		lAlg.set_maxeval(localNLParams.maxIts);
		lAlg.set_maxtime(localNLParams.maxTimeSec);

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

	//x = [tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45
	KPMS kp;
	double kna[kp.nParams];
	int ix = 0; // since nX <= nParams, need a separate index
	for (int i = 0; i < kp.nParams; i++) {
		if (pfip->knSubset[i]) {
			kna[i] = x[ix]; //kn params follow the q0s
			ix++;
		}else {
			kna[i] = pfip->knDefault[i];
		}
	}
	
	// FK with the new params
	TFK fk(kna);

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
		//printf("%d: qps[%f %f %f %f %f] - Hnew[%f %f %f %f %f %f] evec[%f %f %f %f %f %f] = %f = %f\n", isamp, qps[0],qps[1],qps[2],qps[3],qps[4], Hnew(0,3),Hnew(1,3),Hnew(2,3),Hnew(0,0),Hnew(1,0),Hnew(2,0), evec(0), evec(1), evec(2), evec(3), evec(4), evec(5), evec.norm(), accum);
	}
	return accum;
}

template <class KPMS, class TFK>
void InvP_nlopt<KPMS, TFK>::funKn(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;
	fip.knSubset = knSubset; //knSubset is set by the constructor above
	fip.knDefault = kn0;

	// find objective length
	int numX = 0;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i])
			numX++;
	}
	std::vector<double> x(numX);
	int iX = 0;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i]) {
			x[iX] = kn0[i];
			iX++;
		}
	}
	*fmin = err_Kn_xyzuxuyuz<KPMS, TFK>(x, x, &fip) / nSamps;
}

template <class KPMS, class TFK>
int InvP_nlopt<KPMS, TFK>::estimateKn(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;
	fip.knSubset = knSubset;
	fip.knDefault = kn0;

	// find objective length
	int iX = 0;
	int numX = 0;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i])
			numX++;
	}

	// get an alg object
	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), numX);
	
	// set objective
	void *pfip;
	pfip = &fip;
	alg.set_min_objective(err_Kn_xyzuxuyuz<KPMS, TFK>, pfip);
	
	// set initial step size (only used by nongradient methods)
	std::vector<double> iniStep(numX);
	for (int i = 0; i < numX; i++) { iniStep[i] = .001; }
	alg.set_initial_step(iniStep);
	alg.set_default_initial_step(iniStep);
	
	//set boundary constraints
	std::vector<double> limup(numX);
	std::vector<double> limdn(numX);
	double knaup[kup.nParams], knadn[kup.nParams];
	kinematicStruct2Array(&kup, knaup);
	kinematicStruct2Array(&kdn, knadn);
	iX = 0;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i]) {
			limup[iX] = knaup[i];
			limdn[iX] = knadn[i];
			iX++;
		}
	}
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to be within bounds
	std::vector<double> x(numX);
	iX = 0;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i]) {
			x[iX] = kn0[i];
			iX++;
		}
	}
	for (int i = 0; i < numX; i++) {
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
	if (alg.get_algorithm() == nlopt::GN_MLSL_LDS) {
		nlopt::opt lAlg(ipTranslateNLOptAlg(localNLParams.method), alg.get_dimension());
		std::cout << "Local optimization " << lAlg.get_algorithm_name() << std::endl;
		lAlg.set_min_objective(err_QpKn_xyzuxuyuz<KPMS, TFK>, pfip);
		lAlg.set_upper_bounds(limup);
		lAlg.set_lower_bounds(limdn);
		lAlg.set_stopval(localNLParams.minFunVal);
		lAlg.set_ftol_abs(localNLParams.tolFunAbs);
		lAlg.set_xtol_abs(localNLParams.tolXAbs);
		lAlg.set_maxeval(localNLParams.maxIts);
		lAlg.set_maxtime(localNLParams.maxTimeSec);

		alg.set_local_optimizer(lAlg);
	}
	
	//FILE *fid;
	//fopen_s(&fid, "estimate_kn0.txt", "w+");
	//fprintf(fid, "estimateKn with nSamps[%d] numX[%d]\n", nSamps, numX);
	//fprintf(fid, "knSubset(%d)[", kup.nParams);
	//for (int i = 0; i < kup.nParams; i++) {
	//	fprintf(fid, "%d ", knSubset[i]);
	//}
	//fprintf(fid, "]\nlimup[");
	//for (int i = 0; i < numX; i++) {
	//	fprintf(fid, "%f ", limup[i]);
	//}
	//fprintf(fid, "]\n    x[");
	//for (int i = 0; i < numX; i++) {
	//	fprintf(fid, "%f ", x[i]);
	//}
	//fprintf(fid, "]\nlimdn[");
	//for (int i = 0; i < numX; i++) {
	//	fprintf(fid, "%f ", limdn[i]);
	//}
	


	// solve
	nlopt::result res;
	try {
		res = alg.optimize(x, *fmin);
	}
	catch (nlopt::roundoff_limited e) {
		res = nlopt::FAILURE;
		std::cout << "Caught RoundoffLimited: " << e.what() << std::endl;
		//fclose(fid);
	}
	catch (nlopt::forced_stop e) {
		res = nlopt::FAILURE;
		std::cout << "Caught ForcedStop: " << e.what() << std::endl;
		//fclose(fid);
	}
	catch (std::runtime_error e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
		//fclose(fid);
	}
	catch (std::invalid_argument e) {
		res = nlopt::FAILURE;
		printf("Caught: %s for x=[", e.what());
		for (int i = 0; i < numX; i++) { printf("%f ", x[i]); }
		printf("]\n");
		//fclose(fid);
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
		//fclose(fid);
	}
		


	// unpack from solver
	iX = 0;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i]) {
			kn0[i] = x[iX];
			iX++;
		}
	}
	*fmin /= nSamps;

	//fprintf(fid, "]\n  kn1[");
	//for (int i = 0; i < kup.nParams; i++) {
	//	fprintf(fid, "%f ", kn0[i]);
	//}
	//fprintf(fid, "]\n fmin[%f]", *fmin);
	//fclose(fid);

	return res;
}


//simultaneous qp0 & kn0
//xyz only
template <class KPMS, class TFK>
double err_QpKn_xyz(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData;

	double qps[5];
	double accum = 0.0;

	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45
	KPMS kp;
	double kna[kp.nParams];
	int ix = 0; // since nX <= nParams, need a separate index
	for (int i = 0; i < kp.nParams; i++) {
		if (pfip->knSubset[i]) {
			kna[i] = x[5 + ix]; //kn params follow the q0s
			ix++;
		}
		else {
			kna[i] = pfip->knDefault[i];
		}
	}

	// FK with the new params
	TFK fk(kna);

	//calculate error over all samps
	Eigen::Matrix4d Hnew;
	Eigen::VectorXd evec(3);
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

		accum += evec.norm();
	}
	return accum;
}

template <class KPMS, class TFK>
void InvP_nlopt<KPMS, TFK>::funQpKn(int nSamps, double *stackedQ, double *stackedX, double *qp0, double *kn0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedX = stackedX;
	fip.knSubset = knSubset;
	fip.knDefault = kn0;

	std::vector<double> x(nQps + kup.nParams);
	for (int i = 0; i < nQps; i++) { x[i] = qp0[i]; }
	for (int i = 0; i < kup.nParams; i++) { x[i + nQps] = kn0[i]; }
	*fmin = err_QpKn_xyz<KPMS, TFK>(x, x, &fip) / nSamps;
}

template <class KPMS, class TFK>
int InvP_nlopt<KPMS, TFK>::estimateQpKn(int nSamps, double *stackedQ, double *stackedX, double *qp0, double *kn0, double *fmin) {
	//fill optimization data struct
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedX = stackedX;
	fip.knSubset = knSubset;
	fip.knDefault = kn0;

	// find objective length
	int iX = 0;
	int numX = nQps;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i])
			numX++;
	}

	// get an alg object
	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), numX); //there are 5 qps + n( kinematic params in knParams11A )

	// set objective
	void *pfip;
	pfip = &fip;
	alg.set_min_objective(err_QpKn_xyz<KPMS, TFK>, pfip);

	// set initial step size (only used by nongradient methods)
	std::vector<double> iniStep(numX);
	for (int i = 0; i < numX; i++) { iniStep[i] = .001; }
	alg.set_initial_step(iniStep);
	alg.set_default_initial_step(iniStep);

	//set boundary constraints
	std::vector<double> limup(numX);
	std::vector<double> limdn(numX);
	for (int i = 0; i < nQps; i++) {
		limup[i] = q0Lims.up[i];
		limdn[i] = q0Lims.dn[i];
	}
	double knaup[kup.nParams], knadn[kup.nParams];
	kinematicStruct2Array(&kup, knaup);
	kinematicStruct2Array(&kdn, knadn);
	iX = 0;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i]) {
			limup[nQps + iX] = knaup[i];
			limdn[nQps + iX] = knadn[i];
			iX++;
		}
	}
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to be within bounds
	std::vector<double> x(numX);
	for (int i = 0; i < nQps; i++) { x[i] = qp0[i]; }
	iX = 0;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i]) {
			x[nQps + iX] = kn0[i];
			iX++;
		}
	}
	for (int i = 0; i < numX; i++) {
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
	if (alg.get_algorithm() == nlopt::GN_MLSL_LDS) {
		nlopt::opt lAlg(ipTranslateNLOptAlg(localNLParams.method), alg.get_dimension());
std::cout << "Local optimization " << lAlg.get_algorithm_name() << std::endl;
		lAlg.set_min_objective(err_QpKn_xyzuxuyuz<KPMS, TFK>, pfip);
		lAlg.set_upper_bounds(limup);
		lAlg.set_lower_bounds(limdn);
		lAlg.set_stopval(localNLParams.minFunVal);
		lAlg.set_ftol_abs(localNLParams.tolFunAbs);
		lAlg.set_xtol_abs(localNLParams.tolXAbs);
		lAlg.set_maxeval(localNLParams.maxIts);
		lAlg.set_maxtime(localNLParams.maxTimeSec);

		alg.set_local_optimizer(lAlg);
	}

	//FILE *fid;
	//fopen_s(&fid, "estimate_qp0kn0_xyz.txt", "w+");
	//fprintf(fid, "estimateKn with nSamps[%d] numX[%d]\n", nSamps, numX);
	//fprintf(fid, "knSubset(%d)[", kup.nParams);
	//for (int i = 0; i < kup.nParams; i++) {
	//	fprintf(fid, "%d ", knSubset[i]);
	//}
	//fprintf(fid, "]\nlimup[");
	//for (int i = 0; i < numX; i++) {
	//	fprintf(fid, "%f ", limup[i]);
	//}
	//fprintf(fid, "]\n    x[");
	//for (int i = 0; i < numX; i++) {
	//	fprintf(fid, "%f ", x[i]);
	//}
	//fprintf(fid, "]\nlimdn[");
	//for (int i = 0; i < numX; i++) {
	//	fprintf(fid, "%f ", limdn[i]);
	//}


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
		for (int i = 0; i < numX; i++) { printf("%f ", x[i]); }
		printf("]\n");
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}

	// unpack from solver
	for (int i = 0; i < nQps; i++) { qp0[i] = x[i]; }
	iX = 0;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i]) {
			kn0[i] = x[nQps + iX];
			iX++;
		}
	}
	*fmin /= nSamps;

	//fprintf(fid, "]\n  kn1[");
	//for (int i = 0; i < kup.nParams; i++) {
	//	fprintf(fid, "%f ", kn0[i]);
	//}
	//fprintf(fid, "]\n fmin[%f]", *fmin);
	//fclose(fid);

	return res;
}
//xyzuxuyuz
template <class KPMS, class TFK>
double err_QpKn_xyzuxuyuz(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData;

	double qps[5];
	double accum = 0.0;

	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45
	KPMS kp;
	double kna[kp.nParams];
	int ix = 0; // since nX <= nParams, need a separate index
	for (int i = 0; i < kp.nParams; i++) {
		if ( pfip->knSubset[i]) {
			kna[i] = x[5 + ix]; //kn params follow the q0s
			ix++;
		}else {
			kna[i] = pfip->knDefault[i];
		}
	}

	// FK with the new params
	TFK fk(kna);

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
	fip.knSubset = knSubset;
	fip.knDefault = kn0;

	std::vector<double> x(nQps + kup.nParams);
	for (int i = 0; i < nQps; i++) { x[i] = qp0[i]; }
	for (int i = 0; i < kup.nParams; i++) { x[i + nQps] = kn0[i]; }

	*fmin = err_QpKn_xyzuxuyuz<KPMS,TFK>(x, x, &fip) / nSamps;
}

template <class KPMS, class TFK>
int InvP_nlopt<KPMS, TFK>::estimateQpKn(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *kn0, double *fmin) {
	//fill optimization data struct
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;
	fip.knSubset = knSubset;
	fip.knDefault = kn0;
	
	// find objective length
	int iX = 0;
	int numX = nQps;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i])
			numX++;
	}

	// get an alg object
	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), numX);

	// set objective
	void *pfip;
	pfip = &fip;
	alg.set_min_objective(err_QpKn_xyzuxuyuz<KPMS,TFK>, pfip);

	// set initial step size (only used by nongradient methods)
	std::vector<double> iniStep(numX);
	for (int i = 0; i < numX; i++) { iniStep[i] = .001; }
	alg.set_initial_step(iniStep);
	alg.set_default_initial_step(iniStep);

	//set boundary constraints
	std::vector<double> limup(numX);
	std::vector<double> limdn(numX);
	for (int i = 0; i < nQps; i++) {
		limup[i] = q0Lims.up[i];
		limdn[i] = q0Lims.dn[i];
	}
	double knaup[kup.nParams], knadn[kup.nParams];
	kinematicStruct2Array(&kup, knaup);
	kinematicStruct2Array(&kdn, knadn);
	iX = 0;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i]) {
			limup[nQps + iX] = knaup[i];
			limdn[nQps + iX] = knadn[i];
			iX++;
		}
	}
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);
	
	// set start x to be within bounds
	std::vector<double> x(numX);
	for (int i = 0; i < nQps; i++) { x[i] = qp0[i]; }
	iX = 0;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i]) {
			x[nQps + iX] = kn0[i];
			iX++;
		}
	}
	for (int i = 0; i < numX; i++) {
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
	if (alg.get_algorithm() == nlopt::GN_MLSL_LDS) {
		nlopt::opt lAlg(ipTranslateNLOptAlg(localNLParams.method), alg.get_dimension());
std::cout << "Local optimization " << lAlg.get_algorithm_name() << std::endl;
		lAlg.set_min_objective(err_QpKn_xyzuxuyuz<KPMS, TFK>, pfip);
		lAlg.set_upper_bounds(limup);
		lAlg.set_lower_bounds(limdn);
		lAlg.set_stopval(localNLParams.minFunVal);
		lAlg.set_ftol_abs(localNLParams.tolFunAbs);
		lAlg.set_xtol_abs(localNLParams.tolXAbs);
		lAlg.set_maxeval(localNLParams.maxIts);
		lAlg.set_maxtime(localNLParams.maxTimeSec);

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
		for (int i = 0; i < numX; i++) { printf("%f ", x[i]); }
		printf("]\n");
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	
	// unpack from solver
	for (int i = 0; i < nQps; i++) { qp0[i] = x[i]; }
	iX = 0;
	for (int i = 0; i < kup.nParams; i++) {
		if (knSubset[i]) {
			kn0[i] = x[nQps + iX];
			iX++;
		}
	}
	*fmin /= nSamps;

	return res;
}


//Kinematic parameter struct to array
void kinematicStruct2Array(KINEMATICPARAMS5A *kns, double *kna) {
	kna[0] = (*kns).tx01;
	kna[1] = (*kns).ty01;
	kna[2] = (*kns).tz01;
	kna[3] = (*kns).rz01;
	kna[4] = (*kns).lCath;
}
void kinematicStruct2Array(KINEMATICPARAMS6A *kns, double *kna) {
	kna[0] = (*kns).tx01;
	kna[1] = (*kns).ty01;
	kna[2] = (*kns).tz01;
	kna[3] = (*kns).rz01;
	kna[4] = (*kns).ry34;
	kna[5] = (*kns).lCath;
}
void kinematicStruct2Array(KINEMATICPARAMS11A *kns, double *kna) {
	kna[0] = (*kns).tx01;
	kna[1] = (*kns).ty01;
	kna[2] = (*kns).tz01;
	kna[3] = (*kns).ry01;
	kna[4] = (*kns).rz01;
	kna[5] = (*kns).ry34;
	kna[6] = (*kns).rz34;
	kna[7] = (*kns).kAlpha;
	kna[8] = (*kns).eAlpha;
	kna[9] = (*kns).lCath;
	kna[10] = (*kns).ry45;
}

//Kinematic parameter array to struct
void kinematicArray2Struct(double *kna, KINEMATICPARAMS5A *kns) {
	(*kns).tx01 = kna[0];
	(*kns).ty01 = kna[1];
	(*kns).tz01 = kna[2];
	(*kns).rz01 = kna[3];
	(*kns).lCath = kna[4];
}
void kinematicArray2Struct(double *kna, KINEMATICPARAMS6A *kns) {
	(*kns).tx01 = kna[0];
	(*kns).ty01 = kna[1];
	(*kns).tz01 = kna[2];
	(*kns).rz01 = kna[3];
	(*kns).ry34 = kna[4];
	(*kns).lCath = kna[5];
}
void kinematicArray2Struct(double *kna, KINEMATICPARAMS11A *kns) {
	(*kns).tx01 = kna[0];
	(*kns).ty01 = kna[1];
	(*kns).tz01 = kna[2];
	(*kns).ry01 = kna[3];
	(*kns).rz01 = kna[4];
	(*kns).ry34 = kna[5];
	(*kns).rz34 = kna[6];
	(*kns).kAlpha = kna[7];
	(*kns).eAlpha = kna[8];
	(*kns).lCath = kna[9];
	(*kns).ry45 = kna[10];
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

