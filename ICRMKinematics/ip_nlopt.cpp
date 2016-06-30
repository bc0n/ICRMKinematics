#include "ip_nlopt.h"
#include <iostream>


// initial joint angles
double err_qp0_xyz5A(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData; // cast from void

	double qps[5];
	double accum = 0.0;

	//calculate error over all samps
	Eigen::Matrix4d Hnew;
	Eigen::Vector3d evec;
	FwdK5A *fk = (FwdK5A*)(pfip->fk);
	for (int isamp = 0; isamp < pfip->nSamps; isamp++) {
		//update q
		for (int j = 0; j < 5; j++) {
			qps[j] = pfip->stackedQ[isamp * 5 + j] + x[j]; //try out the candidate qp0
		}
		//new position
		Hnew = fk->qps2H05(qps);

		//error
		evec(0) = pfip->stackedX[isamp * 3 + 0] - Hnew(0, 3);
		evec(1) = pfip->stackedX[isamp * 3 + 1] - Hnew(1, 3);
		evec(2) = pfip->stackedX[isamp * 3 + 2] - Hnew(2, 3);
		accum += evec.norm();
	}
	return accum;
}
IPnlopt_qp0_xyz5A::IPnlopt_qp0_xyz5A() {
	KINEMATICPARAMS5A kn;
	NLOPTPARAMS nlParams;
}
IPnlopt_qp0_xyz5A::IPnlopt_qp0_xyz5A(KINEMATICPARAMS5A kn5a) {
	kn = kn5a;
}
IPnlopt_qp0_xyz5A::IPnlopt_qp0_xyz5A(KINEMATICPARAMS5A kn5a, NLOPTPARAMS nl) {
	kn = kn5a;
	nlParams = nl;
}
void IPnlopt_qp0_xyz5A::funIP_qp0_xyz5A(int nSamps, double *stackedQ, double *stackedX, double *qp0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedX = stackedX;

	FwdK5A fk(kn);
	fip.fk = (void*)(&fk);

	std::vector<double> x(5), g(5);
	for (int i = 0; i < 5; i++) { x[i] = qp0[i]; }

	*fmin = err_qp0_xyz5A(x, g, &fip);
}
int IPnlopt_qp0_xyz5A::estimate(int nSamps, double *stackedQ, double *stackedX, double *qp0, double *qpup, double *qpdn, double *fmin) {
	int ret = -99;
	
	// set method
	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), 5); //there are 5 kinematic params in knParams5A

	// set objective
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedX = stackedX;
	FwdK5A fk(kn);
	fip.fk = (void*)&fk;
	void *vfip = &fip;
	alg.set_min_objective(err_qp0_xyz5A, vfip);

	// set initial step size (only used by nongradient methods)
	//

	//set boundary constraints
	std::vector<double> limup(5), limdn(5);
	for (int i = 0; i < 5; i++) {
		limup[i] = qpup[i];
		limdn[i] = qpdn[i];
	}
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to be within bounds
	std::vector<double> x(5);
	for (int i = 0; i < 5; i++) {
		x[i] = qp0[i];
		if (x[i] <= limdn[i]) {
			x[i] = limdn[i];
		}
		if (limup[i] <= x[i]) {
			x[i] = limup[i];
		}
	}
	//printf("Starting Initial Joint Angle Search from qp0[%f,%f,%f,%f,%f]\n", x[0], x[1], x[2], x[3], x[4]);

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	//may need to set the local optimizer, eg MLSL http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
	if (ipTranslateNLOptAlg(nlParams.method) == GN_MLSL || ipTranslateNLOptAlg(nlParams.method) == GN_MLSL_LDS) {
		nlopt::opt lAlg(nlopt::LN_NELDERMEAD, alg.get_dimension());

		lAlg.set_min_objective(err_qp0_xyz5A, vfip);
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
	*fmin = 1e3;

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
		printf("Caught: %s for x=[ %4.4f %4.4f %4.4f %4.4f %4.4f]\n", e.what(), x[0], x[1], x[2], x[3], x[4]);
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	catch (...) {
		std::cout << "Caught ellipsis" << std::endl;
	}
	*fmin /= nSamps;
	//unpack from solver
	for (int i = 0; i < 5; i++) {
		qp0[i] = x[i];
	}
	//printf("Found err %f and qp0[%f,%f,%f,%f,%f]\n", *fmin, qp0[0], qp0[1], qp0[2], qp0[3], qp0[4]);

	return res;
}

double err_qp0_xyzuxuyuz5A(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData; // cast from void

	double qps[5];
	double accum = 0.0;

	//calculate error over all samps
	Eigen::Matrix4d Hnew;
	Eigen::VectorXd evec(6);
	FwdK5A *fk = (FwdK5A*)(pfip->fk);
	for (int isamp = 0; isamp < pfip->nSamps; isamp++) {
		//update q
		for (int j = 0; j < 5; j++) {
			qps[j] = pfip->stackedQ[isamp * 5 + j] + x[j]; //try out the candidate qp0
		}
		//new position
		Hnew = fk->qps2H05(qps);

		//error
		evec(0) = pfip->stackedX[isamp * 3 + 0] - Hnew(0, 3);//x
		evec(1) = pfip->stackedX[isamp * 3 + 1] - Hnew(1, 3);//y
		evec(2) = pfip->stackedX[isamp * 3 + 2] - Hnew(2, 3);//z
		evec(3) = pfip->stackedU[isamp * 3 + 0] - Hnew(0, 0);//ux
		evec(4) = pfip->stackedU[isamp * 3 + 1] - Hnew(1, 0);//uy
		evec(5) = pfip->stackedU[isamp * 3 + 2] - Hnew(2, 0);//uz
		accum += evec.norm();
	}
	return accum;
}
IPnlopt_qp0_xyzuxuyuz5A::IPnlopt_qp0_xyzuxuyuz5A() {
	KINEMATICPARAMS5A kn;
	NLOPTPARAMS nlParams;
}
IPnlopt_qp0_xyzuxuyuz5A::IPnlopt_qp0_xyzuxuyuz5A(KINEMATICPARAMS5A kn5a) {
	kn = kn5a;
}
IPnlopt_qp0_xyzuxuyuz5A::IPnlopt_qp0_xyzuxuyuz5A(KINEMATICPARAMS5A kn5a, NLOPTPARAMS nl) {
	kn = kn5a;
	nlParams = nl;
}
void IPnlopt_qp0_xyzuxuyuz5A::funIP_qp0_xyzuxuyuz5A(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedX = stackedX;
	fip.stackedU = stackedU;

	FwdK5A fk(kn);
	fip.fk = (void*)(&fk);

	std::vector<double> x(5), g(5);
	for (int i = 0; i < 5; i++) { x[i] = qp0[i]; }

	*fmin = err_qp0_xyzuxuyuz5A(x, g, &fip);
}
int IPnlopt_qp0_xyzuxuyuz5A::estimate(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *qp0, double *qpup, double *qpdn, double *fmin) {
	int ret = -99;

	// set method
	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), 5); //there are 5 kinematic params in knParams5A

															 // set objective
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedX = stackedX;
	fip.stackedU = stackedU;
	FwdK5A fk(kn);
	fip.fk = (void*)&fk;
	void *vfip = &fip;
	alg.set_min_objective(err_qp0_xyzuxuyuz5A, vfip);

	// set initial step size (only used by nongradient methods)
	//

	//set boundary constraints
	std::vector<double> limup(5), limdn(5);
	for (int i = 0; i < 5; i++) {
		limup[i] = qpup[i];
		limdn[i] = qpdn[i];
	}
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to be within bounds
	std::vector<double> x(5);
	for (int i = 0; i < 5; i++) {
		x[i] = qp0[i];
		if (x[i] <= limdn[i]) {
			x[i] = limdn[i];
		}
		if (limup[i] <= x[i]) {
			x[i] = limup[i];
		}
	}
	//printf("Starting Initial Joint Angle Search from qp0[%f,%f,%f,%f,%f]\n", x[0], x[1], x[2], x[3], x[4]);

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	//may need to set the local optimizer, eg MLSL http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
	if (ipTranslateNLOptAlg(nlParams.method) == GN_MLSL || ipTranslateNLOptAlg(nlParams.method) == GN_MLSL_LDS) {
		nlopt::opt lAlg(nlopt::LN_NELDERMEAD, alg.get_dimension());

		lAlg.set_min_objective(err_qp0_xyzuxuyuz5A, vfip);
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
	*fmin = 1e3;

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
		printf("Caught: %s for x=[ %4.4f %4.4f %4.4f %4.4f %4.4f]\n", e.what(), x[0], x[1], x[2], x[3], x[4]);
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	catch (...) {
		std::cout << "Caught ellipsis" << std::endl;
	}
	*fmin /= nSamps;
	//unpack from solver
	for (int i = 0; i < 5; i++) {
		qp0[i] = x[i];
	}
	//printf("Found err %f and qp0[%f,%f,%f,%f,%f]\n", *fmin, qp0[0], qp0[1], qp0[2], qp0[3], qp0[4]);

	return res;
}

// inverse parameter
double err_kn0_xyz5A(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData;

	double qps[5];
	double accum = 0.0;

	//        0    1    2    3    4     5    6    7      8     9   10
	//x = [tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,lCath, ry45,rz45

	//new kinematic params
	KINEMATICPARAMS5A kp;
	kp.tx01 = x[0];
	kp.ty01 = x[1];
	kp.tz01 = x[2];
	kp.rz01 = x[3];
	kp.lCath = x[4];
	
	// FK with the new params
	FwdK5A fk(kp);

	//calculate error over all samps
	Eigen::Matrix4d Hnew;
	Eigen::Vector3d xyz;
	Eigen::Vector3d evec;
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
		accum += evec.norm();

		//for verifiying against matlab
		//printf("%d: q[%+5.2f,%+5.2f,%+5.2f,%+5.2f,%+5.2f]", isamp, qps[0], qps[1], qps[2], qps[3], qps[4]);
		//printf(" u[%5.2f,%5.2f,%5.2f] x[%5.2f,%5.2f,%5.2f]", uxyz(0),uxyz(1),uxyz(2), xyz(0),xyz(1),xyz(2));
		//printf("\n");
		//printf(" u[%5.2f,%5.2f,%5.2f] x[%5.2f,%5.2f,%5.2f] = %f\n", unew(0), unew(1), unew(2), Hnew(0, 3), Hnew(1, 3), Hnew(2, 3), accum);

	}
	return accum;
}
IPnlopt_kn0_xyz5A::IPnlopt_kn0_xyz5A() {
	KINEMATICPARAMS5A k;
	k5up = k;
	k5dn = k;
	NLOPTPARAMS n;
	nlParams = n;
}
IPnlopt_kn0_xyz5A::IPnlopt_kn0_xyz5A(KINEMATICPARAMS5A kup, KINEMATICPARAMS5A kdn, NLOPTPARAMS nl) {
	k5up = kup;
	k5dn = kdn;
	nlParams = nl;
}
void IPnlopt_kn0_xyz5A::funIP_kn0_xyz5A(int nSamps, double *stackedQ, double *stackedX, double *kn0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedX = stackedX;

	std::vector<double> x(5), g(5);
	for (int i = 0; i < 5; i++) { x[i] = kn0[i]; }

	*fmin = err_kn0_xyz5A(x, g, &fip);
}
int IPnlopt_kn0_xyz5A::estimate(int nSamps, double *stackedQ, double *stackedX, double *kn0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedX = stackedX;

	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), 5); //there are 5 kinematic params in knParams5A

	// set objective
	void *pfip;
	pfip = &fip;
	alg.set_min_objective(err_kn0_xyz5A, pfip);

	// set initial step size (only used by nongradient methods)

	//set boundary constraints
	//      0  1  2  3  4
	//x = [tx01,ty01,tz01,rz01,lCath
	std::vector<double> limup(5);
	std::vector<double> limdn(5);
	limup[0] = k5up.tx01;
	limup[1] = k5up.ty01;
	limup[2] = k5up.tz01;
	limup[3] = k5up.rz01;
	limup[4] = k5up.lCath;
	limdn[0] = k5dn.tx01;
	limdn[1] = k5dn.ty01;
	limdn[2] = k5dn.tz01;
	limdn[3] = k5dn.rz01;
	limdn[4] = k5dn.lCath;
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	//printf("knU[%f,%f,%f,%f,%f]\n", limup[0], limup[1], limup[2], limup[3], limup[4]);
	//printf("kn0[%f,%f,%f,%f,%f]\n", kn0[0], kn0[1], kn0[2], kn0[3], kn0[4]);
	//printf("knD[%f,%f,%f,%f,%f]\n", limdn[0], limdn[1], limdn[2], limdn[3], limdn[4]);

	// set start x to be within bounds
	std::vector<double> x(5);
	for (int i = 0; i < 5; i++) {
		x[i] = kn0[i];
		if (x[i] <= limdn[i]) {
			x[i] = limdn[i];
		}
		if (limup[i] <= x[i]) {
			x[i] = limup[i];
		}
	}
	//printf("Starting Inverse Parameter Search from kn0[%f,%f,%f,%f,%f]\n", x[0], x[1], x[2], x[3], x[4]);

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	//may need to set the local optimizer, eg MLSL http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
	if (ipTranslateNLOptAlg(nlParams.method) == GN_MLSL || ipTranslateNLOptAlg(nlParams.method) == GN_MLSL_LDS) {
		nlopt::opt lAlg(nlopt::LN_NELDERMEAD, alg.get_dimension());

		lAlg.set_min_objective(err_kn0_xyz5A, pfip);
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
		printf("Caught: %s for x=[ %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f]\n", e.what(), x[0], x[1], x[2], x[3], x[4], x[5]);
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	*fmin /= nSamps;
	// unpack from solver
	for (int i = 0; i < 5; i++) {
		kn0[i] = x[i];
	}
	//printf("Found err %f and kn0[%f,%f,%f,%f,%f]\n", fmn, kn0[0], kn0[1], kn0[2], kn0[3], kn0[4]);

	return res;
}


// simultaneous qp0 kn0
double err_qp0kn0_xyz5A(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData;

	double qps[5];
	double accum = 0.0;

	//      0  1  2  3  4     5    6    7    8     9
	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,rz01,lCath

	//new kinematic params
	KINEMATICPARAMS5A kp;
	kp.tx01 =  x[5];
	kp.ty01 =  x[6];
	kp.tz01 =  x[7];
	kp.rz01 =  x[8];
	kp.lCath = x[9];

	// FK with the new params
	FwdK5A fk(kp);

	//calculate error over all samps
	Eigen::Matrix4d Hnew;
	Eigen::Vector3d xyz;
	Eigen::Vector3d evec;
	for (int isamp = 0; isamp < pfip->nSamps; isamp++) {
		//update q
		for (int j = 0; j < 5; j++) {
			qps[j] = pfip->stackedQ[isamp * 5 + j] + x[j]; //new qp0
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
IPnlopt_qp0kn0_xyz5A::IPnlopt_qp0kn0_xyz5A() {
	KINEMATICPARAMS5A k;
	k5up = k;
	k5dn = k;
	JOINTLIMITS qp0Lim;
	NLOPTPARAMS nlParams;
}
IPnlopt_qp0kn0_xyz5A::IPnlopt_qp0kn0_xyz5A(KINEMATICPARAMS5A kup, KINEMATICPARAMS5A kdn, JOINTLIMITS qp0Limits, NLOPTPARAMS nl) {
	k5up = kup;
	k5dn = kdn;
	qp0Lim = qp0Limits;
	nlParams = nl;
}
void IPnlopt_qp0kn0_xyz5A::funIP_qp0kn0_xyz5A(int nSamps, double *stackedQ, double *stackedX, double *kn0, double *qp0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedX = stackedX;

	std::vector<double> x(10), g(10);
	for (int i = 0; i < 5; i++) {
		x[i] = qp0[i];
		x[5 + i] = kn0[i];
	}

	*fmin = err_qp0kn0_xyz5A(x, g, &fip);
}
int IPnlopt_qp0kn0_xyz5A::estimate(int nSamps, double *stackedQ, double *stackedX, double *kn0, double *qp0, double *fmin) {

	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedX = stackedX;

	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), 10); //there are 5 qp0 + 5 kinematic params

	// set objective
	void *pfip;
	pfip = &fip;
	alg.set_min_objective(err_qp0kn0_xyz5A, pfip);

	// set initial step size (only used by nongradient methods)

	//set boundary constraints
	//      0  1  2  3  4     5    6    7    8     9
	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,rz01,lCath
	std::vector<double> limup(10);
	std::vector<double> limdn(10);
	limup[0] = qp0Lim.up[0];
	limup[1] = qp0Lim.up[1];
	limup[2] = qp0Lim.up[2];
	limup[3] = qp0Lim.up[3];
	limup[4] = qp0Lim.up[4];
	limup[5] = k5up.tx01;
	limup[6] = k5up.ty01;
	limup[7] = k5up.tz01;
	limup[8] = k5up.rz01;
	limup[9] = k5up.lCath;
	limdn[0] = qp0Lim.dn[0];
	limdn[1] = qp0Lim.dn[1];
	limdn[2] = qp0Lim.dn[2];
	limdn[3] = qp0Lim.dn[3];
	limdn[4] = qp0Lim.dn[4];
	limdn[5] = k5dn.tx01;
	limdn[6] = k5dn.ty01;
	limdn[7] = k5dn.tz01;
	limdn[8] = k5dn.rz01;
	limdn[9] = k5dn.lCath;
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to be within bounds
	std::vector<double> x(10);
	for (int i = 0; i < 5; i++) {
		x[i] = qp0[i];
		x[5 + i] = kn0[i];
		if (x[i] <= limdn[i]) {
			x[i] = limdn[i];
		}
		if (limup[i] <= x[i]) {
			x[i] = limup[i];
		}
	}
	//printf("Searching for qp0 & kn0 from in\n");
	//printf("up["); for (int i = 0; i < 10; i++) { printf("%8.3f ", limup[i]); } printf("]\n");
	//printf("x0["); for (int i = 0; i < 10; i++) { printf("%8.3f ", x[i]); } printf("]\n");
	//printf("dn["); for (int i = 0; i < 10; i++) { printf("%8.3f ", limdn[i]); } printf("]\n");

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	//may need to set the local optimizer, eg MLSL http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
	if (ipTranslateNLOptAlg(nlParams.method) == GN_MLSL || ipTranslateNLOptAlg(nlParams.method) == GN_MLSL_LDS) {
		nlopt::opt lAlg(nlopt::LN_NELDERMEAD, alg.get_dimension());

		lAlg.set_min_objective(err_qp0kn0_xyz5A, pfip);
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
		printf("Caught: %s for x=[ %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f]\n", e.what(), x[0], x[1], x[2], x[3], x[4], x[5]);
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	*fmin /= nSamps;
	// unpack from solver
	for (int i = 0; i < 5; i++) {
		qp0[i] = x[i];
		kn0[i] = x[5+i];
	}
	//printf("Found err %f and qp0[", *fmin); for (int i = 0; i < 5; i++) { printf("%8.3f ", qp0[i]); }
	//printf("] kn0["); for (int i = 0; i < 5; i++) { printf("%8.3f ", kn0[i]); } printf("]\n");

	return res;
}



//nlopt complains if a member of InvP11
double funIPUX11A(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {

	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData;

	double qps[5];
	double accum = 0.0;

	//      0  1  2  3  4     5    6    7    8    9    10    11   12    13    14   15
	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45

	//new kinematic params
	KINEMATICPARAMS11A kp;
	kp.tx01 = x[5];
	kp.ty01 = x[6];
	kp.tz01 = x[7];
	kp.ry01 = x[8];
	kp.rz01 = x[9];
	kp.ry34 = x[10];
	kp.rz34 = x[11];
	kp.kAlpha = x[12];
	kp.eAlpha = x[13];
	kp.lCath = x[14];
	kp.ry45 = x[15];

	// FK with the new params
	FwdK11A fk(kp);

	//calculate error over all samps
	Eigen::Matrix4d Hnew;
	Eigen::Vector3d uxyz;
	Eigen::Vector3d unew;
	Eigen::Vector3d xyz;
	for (int isamp = 0; isamp < pfip->nSamps; isamp++) {
		//update q
		for (int j = 0; j < 5; j++) {
			qps[j] = x[j] + pfip->stackedQ[isamp * 5 + j];
		}
		//printf("%d: %f %f %f %f %f\n", isamp, qps[0], qps[1], qps[2], qps[3], qps[4]);

		//new position
		Hnew = fk.qps2H05(qps);

		//dot of U and Unew
		uxyz << pfip->stackedU[isamp * 3 + 0], pfip->stackedU[isamp * 3 + 1], pfip->stackedU[isamp * 3 + 2];
		unew << Hnew(0, 0), Hnew(1, 0), Hnew(2, 0);
		accum += 1 - uxyz.dot(unew); //dot == 1 if aligned, so 1-()

		// norm distance
		xyz << pfip->stackedX[isamp * 3 + 0], pfip->stackedX[isamp * 3 + 1], pfip->stackedX[isamp * 3 + 2];
		accum += (xyz - Hnew.block(0, 3, 3, 1)).norm();

		//for verifiying against matlab
		//printf("%d: q[%+5.2f,%+5.2f,%+5.2f,%+5.2f,%+5.2f]", isamp, qps[0], qps[1], qps[2], qps[3], qps[4]);
		//printf(" u[%5.2f,%5.2f,%5.2f] x[%5.2f,%5.2f,%5.2f]", uxyz(0),uxyz(1),uxyz(2), xyz(0),xyz(1),xyz(2));
		//printf("\n");
		//printf(" u[%5.2f,%5.2f,%5.2f] x[%5.2f,%5.2f,%5.2f] = %f,%f = %f\n", unew(0), unew(1), unew(2), Hnew(0, 3), Hnew(1, 3), Hnew(2, 3), 1 - uxyz.dot(unew), (xyz - Hnew.block(0, 3, 3, 1)).norm(), accum);
	}
	return accum;
}
InvPNLOpt_xyzdotu11A::InvPNLOpt_xyzdotu11A() {
	KINEMATICPARAMS11A k;
	k11up = k;
	k11dn = k;
	JOINTLIMITS j;
	q0Lims = j;
	NLOPTPARAMS n;
	nlParams = n;
}
InvPNLOpt_xyzdotu11A::InvPNLOpt_xyzdotu11A(KINEMATICPARAMS11A kup, KINEMATICPARAMS11A kdn, JOINTLIMITS q0, NLOPTPARAMS nl) {
	k11up = kup;
	k11dn = kdn;
	q0Lims = q0;
	nlParams = nl;
}
void InvPNLOpt_xyzdotu11A::funIP_xyzdotu11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *kn0, double *qp0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	std::vector<double> x(5 + 11);
	for (int i = 0; i < 5; i++) { x[i] = qp0[i]; }
	for (int i = 0; i < 11; i++) { x[i + 5] = kn0[i]; }

	*fmin = funIPUX11A(x, x, &fip);
}
int InvPNLOpt_xyzdotu11A::estimate(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *qp0, double *fmin) {
	//printf("Starting Inverse Parameter Search from q0[");
	//for (int i = 0; i < 5; i++) { printf("%8.3f ", qp0[i]); } printf("] kn[");
	//for (int i = 0; i < 11; i++) { printf("%8.3f ", kn0[i]); } printf("]\n");

	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), 16); //there are 5 qps + 11 kinematic params in knParams11A

															  // set objective
	void *pfip;
	pfip = &fip;
	alg.set_min_objective(funIPUX11A, pfip);

	// set initial step size (only used by nongradient methods)

	//set boundary constraints based on qpLast and stdev
	//      0  1  2  3  4     5    6    7    8    9    10    11   12    13    14   15
	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45
	std::vector<double> limup(5 + 11);
	std::vector<double> limdn(5 + 11);
	for (int i = 0; i < 5; i++) {
		limup[i] = q0Lims.up[i];
		limdn[i] = q0Lims.dn[i];
	}
	limup[5] = k11up.tx01;
	limup[6] = k11up.ty01;
	limup[7] = k11up.tz01;
	limup[8] = k11up.ry01;
	limup[9] = k11up.rz01;
	limup[10] = k11up.ry34;
	limup[11] = k11up.rz34;
	limup[12] = k11up.kAlpha;
	limup[13] = k11up.eAlpha;
	limup[14] = k11up.lCath;
	limup[15] = k11up.ry45;
	limdn[5] = k11dn.tx01;
	limdn[6] = k11dn.ty01;
	limdn[7] = k11dn.tz01;
	limdn[8] = k11dn.ry01;
	limdn[9] = k11dn.rz01;
	limdn[10] = k11dn.ry34;
	limdn[11] = k11dn.rz34;
	limdn[12] = k11dn.kAlpha;
	limdn[13] = k11dn.eAlpha;
	limdn[14] = k11dn.lCath;
	limdn[15] = k11dn.ry45;
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to be within bounds
	std::vector<double> x(5 + 11);
	for (int i = 0; i < 5; i++) { x[i] = qp0[i]; }
	for (int i = 0; i < 11; i++) { x[5+i] = kn0[i]; }
	for (int i = 0; i < 5 + 11; i++) {
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

		lAlg.set_min_objective(funIPUX11A, pfip);
		lAlg.set_upper_bounds(limup);
		lAlg.set_lower_bounds(limdn);
		lAlg.set_stopval(nlParams.minFunVal);
		lAlg.set_ftol_abs(nlParams.tolFunAbs);
		lAlg.set_xtol_abs(nlParams.tolXAbs);
		lAlg.set_maxeval(nlParams.maxIts);
		lAlg.set_maxtime(nlParams.maxTimeSec);

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
		for (int i = 0; i < 5 + 11; i++) { printf("%5.3f, ", x[i]); } printf("\n");
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	*fmin /= nSamps;
	//unpack from solver
	for (int i = 0; i < 5; i++) { qp0[i] = x[i]; }
	for (int i = 0; i < 11; i++) { kn0[i] = x[5+i]; }
	//printf("Found err %8.3fq0[", *fmin);
	//for (int i = 0; i < 5; i++) { printf("%8.3f ", qp0[i]); } printf("] kn[");
	//for (int i = 0; i < 11; i++) { printf("%8.3f ", kn0[i]); } printf("]\n");

	return res;
}

double funIPxyzpp11A(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData;

	double qps[5];
	double accum = 0.0;

	//      0  1  2  3  4     5    6    7    8    9    10    11   12    13    14   15
	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45

	//new kinematic params
	KINEMATICPARAMS11A kp;
	kp.tx01 = x[5];
	kp.ty01 = x[6];
	kp.tz01 = x[7];
	kp.ry01 = x[8];
	kp.rz01 = x[9];
	kp.ry34 = x[10];
	kp.rz34 = x[11];
	kp.kAlpha = x[12];
	kp.eAlpha = x[13];
	kp.lCath = x[14];
	kp.ry45 = x[15];

	// FK with the new params
	FwdK11A fk(kp);

	//calculate error over all samps
	Eigen::Matrix4d Hnew;
	Eigen::Vector3d uxyz;
	Eigen::Vector3d unew;
	Eigen::Vector3d xyz;
	Eigen::VectorXd evec(5);
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
		evec(3) = atan2(pfip->stackedU[isamp * 3 + 2], sqrt(pow(pfip->stackedU[isamp * 3 + 0], 2) + pow(pfip->stackedU[isamp * 3 + 1], 2)));
		evec(3) -= atan2(Hnew(2, 0), sqrt(pow(Hnew(0, 0), 2) + pow(Hnew(1, 0), 2))); //z-xy plane
		evec(4) = atan2(pfip->stackedU[isamp * 3 + 1], pfip->stackedU[isamp * 3 + 0]);
		evec(4) -= atan2(Hnew(1, 0), Hnew(0, 0)); //x-y plane

		accum += evec.norm();

		//for verifiying against matlab
		//printf("%d: q[%+5.2f,%+5.2f,%+5.2f,%+5.2f,%+5.2f]", isamp, qps[0], qps[1], qps[2], qps[3], qps[4]);
		//printf(" u[%5.2f,%5.2f,%5.2f] x[%5.2f,%5.2f,%5.2f]", uxyz(0),uxyz(1),uxyz(2), xyz(0),xyz(1),xyz(2));
		//printf("\n");
		//printf(" u[%5.2f,%5.2f,%5.2f] x[%5.2f,%5.2f,%5.2f] = %f\n", unew(0), unew(1), unew(2), Hnew(0, 3), Hnew(1, 3), Hnew(2, 3), accum);

	}
	return accum;
}
InvPNLOpt_xyzpp11A::InvPNLOpt_xyzpp11A() {
	KINEMATICPARAMS11A k;
	k11up = k;
	k11dn = k;
	JOINTLIMITS j;
	q0Lims = j;
	NLOPTPARAMS n;
	nlParams = n;
}
InvPNLOpt_xyzpp11A::InvPNLOpt_xyzpp11A(KINEMATICPARAMS11A kup, KINEMATICPARAMS11A kdn, JOINTLIMITS q0, NLOPTPARAMS nl) {
	k11up = kup;
	k11dn = kdn;
	q0Lims = q0;
	nlParams = nl;
}
void InvPNLOpt_xyzpp11A::funIP_xyzpp11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *kn0, double *qp0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	std::vector<double> x(5 + 11);
	for (int i = 0; i < 5; i++) { x[i] = qp0[i]; }
	for (int i = 0; i < 11; i++) { x[i + 5] = kn0[i]; }

	*fmin = funIPxyzpp11A(x, x, &fip);
}
int InvPNLOpt_xyzpp11A::estimate(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *qp0, double *fmin) {

	//printf("Starting Inverse Parameter Search from q0[");
	//for (int i = 0; i < 5; i++) { printf("%8.3f ", qp0[i]); } printf("] kn[");
	//for (int i = 0; i < 11; i++) { printf("%8.3f ", kn0[i]); } printf("]\n");

	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), 16); //there are 5 qps + 11 kinematic params in knParams11A

															  // set objective
	void *pfip;
	pfip = &fip;
	alg.set_min_objective(funIPxyzpp11A, pfip);

	// set initial step size (only used by nongradient methods)

	//set boundary constraints based on qpLast and stdev
	//      0  1  2  3  4     5    6    7    8    9    10    11   12    13    14   15
	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45
	std::vector<double> limup(5 + 11);
	std::vector<double> limdn(5 + 11);
	for (int i = 0; i < 5; i++) {
		limup[i] = q0Lims.up[i];
		limdn[i] = q0Lims.dn[i];
	}
	limup[5] = k11up.tx01;
	limup[6] = k11up.ty01;
	limup[7] = k11up.tz01;
	limup[8] = k11up.ry01;
	limup[9] = k11up.rz01;
	limup[10] = k11up.ry34;
	limup[11] = k11up.rz34;
	limup[12] = k11up.kAlpha;
	limup[13] = k11up.eAlpha;
	limup[14] = k11up.lCath;
	limup[15] = k11up.ry45;
	limdn[5] = k11dn.tx01;
	limdn[6] = k11dn.ty01;
	limdn[7] = k11dn.tz01;
	limdn[8] = k11dn.ry01;
	limdn[9] = k11dn.rz01;
	limdn[10] = k11dn.ry34;
	limdn[11] = k11dn.rz34;
	limdn[12] = k11dn.kAlpha;
	limdn[13] = k11dn.eAlpha;
	limdn[14] = k11dn.lCath;
	limdn[15] = k11dn.ry45;
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to be within bounds
	std::vector<double> x(5 + 11);
	for (int i = 0; i < 5; i++) { x[i] = qp0[i]; }
	for (int i = 0; i < 11; i++) { x[5 + i] = kn0[i]; }
	for (int i = 0; i < 5 + 11; i++) {
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

		lAlg.set_min_objective(funIPUX11A, pfip);
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
		printf("Caught: %s for x=[ %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f]\n", e.what(), x[0], x[1], x[2], x[3], x[4], x[5]);
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	*fmin /= nSamps;
	// unpack from solver
	for (int i = 0; i < 5; i++) { qp0[i] = x[i]; }
	for (int i = 0; i < 11; i++) { kn0[i] = x[5 + i]; }
	//printf("Found err %8.3f q0[", *fmin);
	//for (int i = 0; i < 5; i++) { printf("%8.3f ", qp0[i]); } printf("] kn[");
	//for (int i = 0; i < 11; i++) { printf("%8.3f ", kn0[i]); } printf("]\n");

	return res;
}

double err_qp0kn0_xyzuxuyuz5A(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData;

	double qps[5];
	double accum = 0.0;

	//      0  1  2  3  4     5    6    7    8      9
	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,rz01, cathL

	//new kinematic params
	KINEMATICPARAMS5A kp;
	kp.tx01 = x[5];
	kp.ty01 = x[6];
	kp.tz01 = x[7];
	kp.rz01 = x[8];
	kp.lCath = x[9];

	// FK with the new params
	FwdK5A fk(kp);

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
InvPNLOpt_xyzuxuyuz5A::InvPNLOpt_xyzuxuyuz5A() {
	KINEMATICPARAMS5A k;
	k5up = k;
	k5dn = k;
	JOINTLIMITS j;
	q0Lims = j;
	NLOPTPARAMS n;
	nlParams = n;
}
InvPNLOpt_xyzuxuyuz5A::InvPNLOpt_xyzuxuyuz5A(KINEMATICPARAMS5A kup, KINEMATICPARAMS5A kdn, JOINTLIMITS q0, NLOPTPARAMS nl) {
	k5up = kup;
	k5dn = kdn;
	q0Lims = q0;
	nlParams = nl;
}
void InvPNLOpt_xyzuxuyuz5A::funIP_xyzuxuyuz5A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *kn0, double *qp0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	std::vector<double> x(5 + 5);
	for (int i = 0; i < 5; i++) { x[i] = qp0[i]; }
	for (int i = 0; i < 5; i++) { x[i + 5] = kn0[i]; }

	*fmin = err_qp0kn0_xyzuxuyuz5A(x, x, &fip);
}
int InvPNLOpt_xyzuxuyuz5A::estimate(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *qp0, double *fmin) {

	//printf("Starting Inverse Parameter Search from q0[");
	//for (int i = 0; i < 5; i++) { printf("%8.3f ", qp0[i]); } printf("] kn[");
	//for (int i = 0; i < 5; i++) { printf("%8.3f ", kn0[i]); } printf("]\n");

	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), 5+5); //there are 5 qps + 5 kinematic params in knParams5a

															  // set objective
	void *pfip;
	pfip = &fip;
	alg.set_min_objective(err_qp0kn0_xyzuxuyuz5A, pfip);

	// set initial step size (only used by nongradient methods)

	//set boundary constraints based on qpLast and stdev
	//      0  1  2  3  4     5    6    7    8      9
	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,rz01, cathL
	std::vector<double> limup(5 + 5);
	std::vector<double> limdn(5 + 5);
	for (int i = 0; i < 5; i++) {
		limup[i] = q0Lims.up[i];
		limdn[i] = q0Lims.dn[i];
	}
	limup[5] = k5up.tx01;
	limup[6] = k5up.ty01;
	limup[7] = k5up.tz01;
	limup[8] = k5up.rz01;
	limup[9] = k5up.lCath;
	limdn[5] = k5dn.tx01;
	limdn[6] = k5dn.ty01;
	limdn[7] = k5dn.tz01;
	limdn[8] = k5dn.rz01;
	limdn[9] = k5dn.lCath;
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to be within bounds
	std::vector<double> x(5 + 5);
	for (int i = 0; i < 5; i++) {
		x[i] = qp0[i];
		x[5 + i] = kn0[i];
	}
	for (int i = 0; i < 5 + 5; i++) {
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

		lAlg.set_min_objective(err_qp0kn0_xyzuxuyuz5A, pfip);
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
		printf("Caught: %s for x=[ %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f]\n", e.what(), x[0], x[1], x[2], x[3], x[4], x[5]);
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	*fmin /= nSamps;

	// unpack from solver
	for (int i = 0; i < 5; i++) { 
		qp0[i] = x[i];
		kn0[i] = x[5 + i];
	}
	//printf("Found err %8.3f q0[", *fmin);
	//for (int i = 0; i < 5; i++) { printf("%8.3f ", qp0[i]); } printf("] kn[");
	//for (int i = 0; i < 11; i++) { printf("%8.3f ", kn0[i]); } printf("]\n");

	return res;
}


double err_xyzuxuyuz11A(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
	FIPDATA *pfip;
	pfip = (FIPDATA*)fipData;

	double qps[5];
	double accum = 0.0;

	//      0  1  2  3  4     5    6    7    8    9    10    11   12    13    14   15
	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45

	//new kinematic params
	KINEMATICPARAMS11A kp;
	kp.tx01 = x[5];
	kp.ty01 = x[6];
	kp.tz01 = x[7];
	kp.ry01 = x[8];
	kp.rz01 = x[9];
	kp.ry34 = x[10];
	kp.rz34 = x[11];
	kp.kAlpha = x[12];
	kp.eAlpha = x[13];
	kp.lCath = x[14];
	kp.ry45 = x[15];

	// FK with the new params
	FwdK11A fk(kp);

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
InvPNLOpt_xyzuxuyuz11A::InvPNLOpt_xyzuxuyuz11A() {
	KINEMATICPARAMS11A k;
	k11up = k;
	k11dn = k;
	JOINTLIMITS j;
	q0Lims = j;
	NLOPTPARAMS n;
	nlParams = n;
}
InvPNLOpt_xyzuxuyuz11A::InvPNLOpt_xyzuxuyuz11A(KINEMATICPARAMS11A kup, KINEMATICPARAMS11A kdn, JOINTLIMITS q0, NLOPTPARAMS nl) {
	k11up = kup;
	k11dn = kdn;
	q0Lims = q0;
	nlParams = nl;
}
void InvPNLOpt_xyzuxuyuz11A::funIP_xyzuxuyuz11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *kn0, double *qp0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	std::vector<double> x(5 + 11);
	for (int i = 0; i < 5; i++) { x[i] = qp0[i]; }
	for (int i = 0; i < 11; i++) { x[i + 5] = kn0[i]; }

	*fmin = err_xyzuxuyuz11A(x, x, &fip);
}
int InvPNLOpt_xyzuxuyuz11A::estimate(int nSamps, double *stackedQ, double *stackedX, double *stackedU, double *kn0, double *qp0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	nlopt::opt alg(ipTranslateNLOptAlg(nlParams.method), 16); //there are 5 qps + 11 kinematic params in knParams11A

															  // set objective
	void *pfip;
	pfip = &fip;
	alg.set_min_objective(err_xyzuxuyuz11A, pfip);

	// set initial step size (only used by nongradient methods)

	//set boundary constraints based on qpLast and stdev
	//      0  1  2  3  4     5    6    7    8    9    10    11   12    13    14   15
	//x = [q0,q1,q2,q3,q4, tx01,ty01,tz01,ry01,rz01, tx23, ry34,rz34,cathL, ry45,rz45
	std::vector<double> limup(5 + 11);
	std::vector<double> limdn(5 + 11);
	for (int i = 0; i < 5; i++) {
		limup[i] = q0Lims.up[i];
		limdn[i] = q0Lims.dn[i];
	}
	limup[5] = k11up.tx01;
	limup[6] = k11up.ty01;
	limup[7] = k11up.tz01;
	limup[8] = k11up.ry01;
	limup[9] = k11up.rz01;
	limup[10] = k11up.ry34;
	limup[11] = k11up.rz34;
	limup[12] = k11up.kAlpha;
	limup[13] = k11up.eAlpha;
	limup[14] = k11up.lCath;
	limup[15] = k11up.ry45;
	limdn[5] = k11dn.tx01;
	limdn[6] = k11dn.ty01;
	limdn[7] = k11dn.tz01;
	limdn[8] = k11dn.ry01;
	limdn[9] = k11dn.rz01;
	limdn[10] = k11dn.ry34;
	limdn[11] = k11dn.rz34;
	limdn[12] = k11dn.kAlpha;
	limdn[13] = k11dn.eAlpha;
	limdn[14] = k11dn.lCath;
	limdn[15] = k11dn.ry45;
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to be within bounds
	std::vector<double> x(5 + 11);
	for (int i = 0; i < 5; i++) { x[i] = qp0[i]; }
	for (int i = 0; i < 11; i++) { x[5 + i] = kn0[i]; }
	for (int i = 0; i < 5 + 11; i++) {
		if (x[i] <= limdn[i]) {
			x[i] = limdn[i];
		}
		if (limup[i] <= x[i]) {
			x[i] = limup[i];
		}
	}
	//printf("Starting Inverse Parameter Search from q0[");
	//printf(" up["); for (int i = 0; i < 5; i++) { printf("%8.3f ", q0Lims.up[i]); } printf("]\n");
	//printf("qp0["); for (int i = 0; i < 5; i++) { printf("%8.3f ", qp0[i]); } printf("]\n");
	//printf(" dn["); for (int i = 0; i < 5; i++) { printf("%8.3f ", q0Lims.dn[i]); } printf("]\n");

	//printf("kn0["); for (int i = 0; i < 11; i++) { printf("%8.3f ", kn0[i]); } printf("]\n");
	//
	//printf("lup["); for (int i = 0; i < 5 + 11; i++) { printf("%8.3f ", limup[i]); }printf("\n");
	//printf("  x["); for (int i = 0; i < 5 + 11; i++) { printf("%8.3f ", x[i]); }printf("\n");
	//printf("ldn["); for (int i = 0; i < 5 + 11; i++) { printf("%8.3f ", limdn[i]); }printf("\n");
	


	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	//may need to set the local optimizer, eg MLSL http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
	if (ipTranslateNLOptAlg(nlParams.method) == GN_MLSL || ipTranslateNLOptAlg(nlParams.method) == GN_MLSL_LDS) {
		nlopt::opt lAlg(nlopt::LN_NELDERMEAD, alg.get_dimension());

		lAlg.set_min_objective(err_xyzuxuyuz11A, pfip);
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
		printf("Caught: %s for x=[ %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f]\n", e.what(), x[0], x[1], x[2], x[3], x[4], x[5]);
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	*fmin /= nSamps;
	// unpack from solver
	for (int i = 0; i < 5; i++) { qp0[i] = x[i]; }
	for (int i = 0; i < 11; i++) { kn0[i] = x[5 + i]; }
	//printf("Found err %8.3f q0[", *fmin);
	//for (int i = 0; i < 5; i++) { printf("%8.3f ", qp0[i]); } printf("] kn[");
	//for (int i = 0; i < 11; i++) { printf("%8.3f ", kn0[i]); } printf("]\n");

	return res;
}



nlopt::algorithm ipTranslateNLOptAlg(nlMethod method) {
	//std::cout << "recv " << method << std::endl;
	switch (method) {
	case nlMethod::GN_DIRECT: return nlopt::GN_DIRECT; //0
	case nlMethod::GN_DIRECT_L: return nlopt::GN_DIRECT;//1
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