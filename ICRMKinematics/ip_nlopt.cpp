#include "ip_nlopt.h"
#include <iostream>

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
	kp.tx23 = x[10];
	kp.ry34 = x[11];
	kp.rz34 = x[12];
	kp.lCath = x[13];
	kp.ry45 = x[14];
	kp.rz45 = x[15];

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
void InvPNLOpt_xyzdotu11A::funIP_UX11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *pms0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	std::vector<double> x(5 + 11);
	for (int i = 0; i < 5; i++) { x[i] = qps0[i]; }
	for (int i = 0; i < 11; i++) { x[i+5] = pms0[i]; }

	*fmin = funIPUX11A(x, x, &fip);
}
int InvPNLOpt_xyzdotu11A::estimate(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *fmin) {
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
	limup[10] = k11up.tx23;
	limup[11] = k11up.ry34;
	limup[12] = k11up.rz34;
	limup[13] = k11up.lCath;
	limup[14] = k11up.ry45;
	limup[15] = k11up.rz45;
	limdn[5] = k11dn.tx01;
	limdn[6] = k11dn.ty01;
	limdn[7] = k11dn.tz01;
	limdn[8] = k11dn.ry01;
	limdn[9] = k11dn.rz01;
	limdn[10] = k11dn.tx23;
	limdn[11] = k11dn.ry34;
	limdn[12] = k11dn.rz34;
	limdn[13] = k11dn.lCath;
	limdn[14] = k11dn.ry45;
	limdn[15] = k11dn.rz45;

	// set start x to be within bounds
	double x0[5 + 11];
	for (int i = 0; i < 5 + 11; i++) {
		x0[i] = (limup[i] - limdn[i]) / 2 + limdn[i];
		if (x0[i] <= limdn[i]) {
			x0[i] = limdn[i];
		}
		if (limup[i] <= x0[i]) {
			x0[i] = limup[i];
		}
	}

	return estimate(nSamps, stackedQ, stackedU, stackedX, x0, fmin);
}
int InvPNLOpt_xyzdotu11A::estimate(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *x0, double *fmin) {
	int ret = -99;
	printf("Starting Inverse Parameter Search from\nq0[%f,%f,%f,%f,%f]\nk0[",x0[0],x0[1],x0[2],x0[3],x0[4]);
	for (int i = 0; i < 11; i++) { printf("%f,", x0[i + 5]); }; printf("]\n");

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
	std::vector<double> limup(5+11);
	std::vector<double> limdn(5+11);
	for (int i = 0; i < 5; i++) {
		limup[i] = q0Lims.up[i];
		limdn[i] = q0Lims.dn[i];
	}
	limup[5] = k11up.tx01;
	limup[6] = k11up.ty01;
	limup[7] = k11up.tz01;
	limup[8] = k11up.ry01;
	limup[9] = k11up.rz01;
	limup[10] = k11up.tx23;
	limup[11] = k11up.ry34;
	limup[12] = k11up.rz34;
	limup[13] = k11up.lCath;
	limup[14] = k11up.ry45;
	limup[15] = k11up.rz45;
	limdn[5] = k11dn.tx01;
	limdn[6] = k11dn.ty01;
	limdn[7] = k11dn.tz01;
	limdn[8] = k11dn.ry01;
	limdn[9] = k11dn.rz01;
	limdn[10] = k11dn.tx23;
	limdn[11] = k11dn.ry34;
	limdn[12] = k11dn.rz34;
	limdn[13] = k11dn.lCath;
	limdn[14] = k11dn.ry45;
	limdn[15] = k11dn.rz45;
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to be within bounds
	std::vector<double> x(5+11);
	for (int i = 0; i < 5+11; i++) {
		x[i] = x0[i]; //(limup[i] - limdn[i]) / 2 + limdn[i];
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
		nlopt::opt lAlg( nlopt::LN_NELDERMEAD, alg.get_dimension());

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
	double fmn = 1e9;
	try {
		res = alg.optimize(x, fmn);
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

	//unpack from solver
	*fmin = fmn;
	ret = res;

	for (int i = 0; i < 5 + 11; i++) {
		x0[i] = x[i];
	}
	printf("Found err %f\nq0[%f,%f,%f,%f,%f]\nk0[", fmn, x0[0], x0[1], x0[2], x0[3], x0[4]);
	for (int i = 0; i < 11; i++) { printf("%f,", x0[i + 5]); }; printf("]\n");

	return ret;
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
	kp.tx23 = x[10];
	kp.ry34 = x[11];
	kp.rz34 = x[12];
	kp.lCath = x[13];
	kp.ry45 = x[14];
	kp.rz45 = x[15];

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

		//new position
		Hnew = fk.qps2H05(qps);

		//error
		Eigen::VectorXd evec(5);
		evec(0) = pfip->stackedX[isamp * 3 + 0] - Hnew(0, 3);
		evec(1) = pfip->stackedX[isamp * 3 + 1] - Hnew(1, 3);
		evec(2) = pfip->stackedX[isamp * 3 + 2] - Hnew(2, 3);
		evec(3) = atan2( pfip->stackedU[isamp*3+2], sqrt(pow(pfip->stackedU[isamp*3+0],2) + pow(pfip->stackedU[isamp*3+1],2)) );
		evec(3) -= atan2( Hnew(2, 0), sqrt(pow(Hnew(0,0),2) + pow(Hnew(1,0),2)) ); //z-xy plane
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
void InvPNLOpt_xyzpp11A::funIP_xyzpp11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *qps0, double *pms0, double *fmin) {
	FIPDATA fip;
	fip.nSamps = (const int)nSamps;
	fip.stackedQ = stackedQ;
	fip.stackedU = stackedU;
	fip.stackedX = stackedX;

	std::vector<double> x(5 + 11);
	for (int i = 0; i < 5; i++) { x[i] = qps0[i]; }
	for (int i = 0; i < 11; i++) { x[i + 5] = pms0[i]; }

	*fmin = funIPxyzpp11A(x, x, &fip);
}
int InvPNLOpt_xyzpp11A::estimate(int nSamps, double *stackedQ, double *stackedU, double *stackedX, double *x0, double *fmin) {
	int ret = -99;
	printf("Starting Inverse Parameter Search from\nq0[%f,%f,%f,%f,%f]\nk0[", x0[0], x0[1], x0[2], x0[3], x0[4]);
	for (int i = 0; i < 11; i++) { printf("%f,", x0[i + 5]); }; printf("]\n");

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
	limup[10] = k11up.tx23;
	limup[11] = k11up.ry34;
	limup[12] = k11up.rz34;
	limup[13] = k11up.lCath;
	limup[14] = k11up.ry45;
	limup[15] = k11up.rz45;
	limdn[5] = k11dn.tx01;
	limdn[6] = k11dn.ty01;
	limdn[7] = k11dn.tz01;
	limdn[8] = k11dn.ry01;
	limdn[9] = k11dn.rz01;
	limdn[10] = k11dn.tx23;
	limdn[11] = k11dn.ry34;
	limdn[12] = k11dn.rz34;
	limdn[13] = k11dn.lCath;
	limdn[14] = k11dn.ry45;
	limdn[15] = k11dn.rz45;
	alg.set_upper_bounds(limup);
	alg.set_lower_bounds(limdn);

	// set start x to be within bounds
	std::vector<double> x(5 + 11);
	for (int i = 0; i < 5 + 11; i++) {
		x[i] = x0[i]; //(limup[i] - limdn[i]) / 2 + limdn[i];
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
	double fmn = 1e9;
	try {
		res = alg.optimize(x, fmn);
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

	// unpack from solver
	*fmin = fmn;
	ret = res;
	for (int i = 0; i < 5 + 11; i++) {
		x0[i] = x[i];
	}
	printf("Found err %f\nq0[%f,%f,%f,%f,%f]\nk0[", fmn, x0[0], x0[1], x0[2], x0[3], x0[4]);
	for (int i = 0; i < 11; i++) { printf("%f,", x0[i + 5]); }; printf("]\n");

	return ret;
}

double errIP_kn0_xyz5A(const std::vector<double> &x, std::vector<double> &grad, void *fipData) {
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
	for (int isamp = 0; isamp < pfip->nSamps; isamp++) {
		//update q
		for (int j = 0; j < 5; j++) {
			qps[j] = pfip->stackedQ[isamp * 5 + j];
		}

		//new position
		Hnew = fk.qps2H05(qps);

		//error
		Eigen::VectorXd evec(5);
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

	std::vector<double> x(5);
	for (int i = 0; i < 5; i++) { x[i] = kn0[i]; }

	*fmin = errIP_kn0_xyz5A(x, x, &fip);
}
int IPnlopt_kn0_xyz5A::estimate(int nSamps, double *stackedQ, double *stackedX, double *kn0, double *fmin) {
	int ret = -99;

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
	double fmn = 1e9;
	try {
		res = alg.optimize(x, fmn);
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

	// unpack from solver
	*fmin = fmn;
	ret = res;
	for (int i = 0; i < 5; i++) {
		kn0[i] = x[i];
	}
	printf("Found err %f and kn0[%f,%f,%f,%f,%f]\n", fmn, kn0[0], kn0[1], kn0[2], kn0[3], kn0[4]);

	return ret;
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