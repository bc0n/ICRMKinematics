
#include "ik_nlopt.h"
#include <iostream>

template <class TASK, class TFK>
double funIK_normTask(const std::vector<double> &x, std::vector<double> &grad, void *fData) {

	FDATA<TASK> *pfdata;
	pfdata = (FDATA<TASK>*)fData; // recover task & tgt

	double qps[5] = { x[0],x[1],x[2],x[3],x[4] };

	//find the current candidate's error
	return (*(pfdata->task)).qps2taskError(qps, pfdata->target);

}

// constructor receives the fk parameters through the templated fkArg
template <class TASK, class TFK>
InvK_nlopt<TASK,TFK>::InvK_nlopt(TFK fkArg, JOINTLIMITS jlArg, NLOPTPARAMS nlArg) {
	tfk = fkArg;
	jntLims = jlArg;
	nlParams = nlArg;
}

template <class TASK, class TFK>
int InvK_nlopt<TASK, TFK>::solve(double *qps, double *xyz){
	int ret = -1;

	//prep algorithm
	TASK task(tfk);
	FDATA<TASK> fdata = { &task, xyz};
	nlopt::opt alg(ikTranslateNLOptAlg(nlParams.method), 5);// 5 = size of x
	
	// set objective
	void *ptr;
	ptr = &fdata;
	alg.set_min_objective(funIK_normTask<TASK,TFK>, ptr);

	// set initial step size (only for nongrad)
	std::vector<double> dx0(5);
	for (int i = 0; i < 5; i++) { dx0[i] = .01; }
	alg.set_initial_step(dx0);

	// set boundary constraints on x
	std::vector<double> limUp(5);
	std::vector<double> limDn(5);
	for (int i = 0; i < 5; i++) {
		limUp[i] = jntLims.up[i];
		limDn[i] = jntLims.dn[i];
	}
	alg.set_upper_bounds(limUp);
	alg.set_lower_bounds(limDn);

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	//may need to set the local optimizer, eg MLSL http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
	if (ikTranslateNLOptAlg(nlParams.method) == GN_MLSL || ikTranslateNLOptAlg(nlParams.method) == GN_MLSL_LDS) {
		nlopt::opt lAlg(nlopt::LN_NELDERMEAD, alg.get_dimension());

		lAlg.set_min_objective(funIK_normTask<TASK,TFK>, ptr);
		lAlg.set_upper_bounds(limUp);
		lAlg.set_lower_bounds(limDn);
		lAlg.set_stopval(nlParams.minFunVal);
		lAlg.set_ftol_abs(nlParams.tolFunAbs);
		lAlg.set_xtol_abs(nlParams.tolXAbs);
		lAlg.set_maxeval(nlParams.maxIts);
		lAlg.set_maxtime(nlParams.maxTimeSec/100); //note same time used for overall and local searches, should probably reduce

		alg.set_local_optimizer(lAlg);
	}

	//set starting location
	std::vector<double> x(5);
	for (int i = 0; i < 5; i++) {
		if (qps[i] > limUp[i]) {
			x[i] = limUp[i];
		}else if (qps[i] < limDn[i]) { 
			x[i] = limDn[i];
		}else {
			x[i] = qps[i];
		}
	}
	//printf("starting at = [%f %f %f %f %f]\n", x[0], x[1], x[2], x[3], x[4]);

	// solve
	nlopt::result res;
	double fmin = 1e3;
	
	//std::cout << "funIK@x =" << funIK_xyz<TFK>(x, x, ptr) << std::endl;
	try {
		res = alg.optimize(x, fmin);
	}
	catch (nlopt::roundoff_limited e) {
		res = nlopt::FAILURE;
		std::cout << "Caught RoundoffLimited: " << e.what() << std::endl;
	}
	catch (nlopt::forced_stop e){
		res = nlopt::FAILURE;
		std::cout << "Caught ForcedStop: " << e.what() << std::endl;
	}
	catch (std::runtime_error e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	catch (std::invalid_argument e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	//printf("alg: res %d    fmin = %f    x = [%f %f %f %f %f]\n", res, fmin, x[0], x[1], x[2], x[3], x[4]);

	//assign outputs
	for (int i = 0; i < 5; i++) { qps[i] = x[i]; }
	task.qps2task(qps, xyz); //run the task on the found qps to capture the found xyz
	ret = res;

	return ret;
}

/*
// Point PhiPsi
double funIK_xyzpp6A(const std::vector<double> &x, std::vector<double> &grad, void *fData) {
	FDATA6A *ptr;
	ptr = (FDATA6A*)fData; // recover pt & tgt

	double qp[5] = { x[0],x[1],x[2],x[3],x[4] };
	Eigen::Vector3d xyzVec;
	ptr->pt.qps2point(qp, &xyzVec); // find current point

	Eigen::Vector2d ppVec;
	ptr->pp.qps2phipsi(qp, &ppVec);
	
	//double fval1 = (ptr->tgtXYZVec - xyzVec).norm();
	//double fval2 = (ptr->tgtPPVec - ppVec).norm();
	//return fval1 + fval2;
	//return fval1 + fval2 / fval1;
	//return fval1 + fval2 * fval1;
	//return fval1 / ptr->tgtXYZVec.norm() + fval2 / ptr->tgtPPVec.norm(); //-2  not necessary since seeking minimum, many close
	//if (fval1 < ptr->nlParams.minFunVal*2) {
	//	return fval1+fval2*800/fval1;
	//	//return fval1 + fval2 / ptr->nlParams.minFunVal;
	//}
	//return fval1;

	Eigen::VectorXd tgt(5);
	tgt << ptr->tgtXYZVec, ptr->tgtPPVec;
	Eigen::VectorXd cur(5);
	cur << xyzVec, ppVec;
	return (tgt - cur).norm(); //14, 18 work here
}
InvKNLOpt_xyzpp6A::InvKNLOpt_xyzpp6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nParams, JOINTLIMITS jLims) {
	kinParams = kParams;
	nlParams = nParams;
	jntLims = jLims;
}
int InvKNLOpt_xyzpp6A::solve(double *qps, double *xyzGoal, double *ppGoal) {
	int ret = -1;
	Point6A pt(kinParams);
	PhiPsi6A pp(kinParams);

	Eigen::VectorXd qpVec(5);
	Eigen::Vector3d curXYZVec;
	Eigen::Vector2d curPPVec;
	Eigen::Vector3d tgtXYZVec;
	Eigen::Vector2d tgtPPVec;
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];
	tgtXYZVec << xyzGoal[0], xyzGoal[1], xyzGoal[2];
	tgtPPVec << ppGoal[0], ppGoal[1];

	FDATA6A ptpptgt;
	ptpptgt.pt = pt;
	ptpptgt.tgtXYZVec = tgtXYZVec;
	ptpptgt.pp = pp;
	ptpptgt.tgtPPVec = tgtPPVec;
	ptpptgt.nlParams = nlParams;

	nlopt::opt alg(ikTranslateNLOptAlg(nlParams.method), 5);

	// set objective
	void *ptr;
	ptr = &ptpptgt;
	alg.set_min_objective(funIK_xyzpp6A, ptr);

	// set initial step size (only for nongrad)
	std::vector<double> dx0(5);
	for (int i = 0; i < 5; i++) { dx0[i] = .01; }
	alg.set_initial_step(dx0);

	// set boundary constraints on x
	std::vector<double> limUp(5);
	std::vector<double> limDn(5);
	for (int i = 0; i < 5; i++) {
		limUp[i] = jntLims.up[i];
		limDn[i] = jntLims.dn[i];
	}
	alg.set_upper_bounds(limUp);
	alg.set_lower_bounds(limDn);

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec); //sec

	// solve
	nlopt::result res;
	double fmin = 1e3;
	std::vector<double> x(5);
	for (int i = 0; i < 5; i++) { x[i] = .1; }
	try {
		res = alg.optimize(x, fmin);
	}
	catch (nlopt::roundoff_limited e) {
		res = nlopt::FAILURE;
		std::cout << "Caught RoundoffLimited: " << e.what() << std::endl;
	}
	catch (nlopt::forced_stop e){
		res = nlopt::FAILURE;
		std::cout << "Caught ForcedStop: " << e.what() << std::endl;
	}
	catch (std::runtime_error e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	catch (std::invalid_argument e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	//printf("alg: res %d    fmin = %f    x = [%f %f %f %f %f]\n", res, fmin, x[0], x[1], x[2], x[3], x[4]);
	for (int i = 0; i < 5; i++) { qpVec(i) = x[i]; }
	pt.qps2point(&qpVec, &curXYZVec);
	pp.qps2phipsi(&qpVec, &curPPVec);
	
	// assign outputs
	Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	Eigen::Map<Eigen::VectorXd>(xyzGoal, curXYZVec.rows()) = curXYZVec;
	Eigen::Map<Eigen::VectorXd>(ppGoal, curPPVec.rows()) = curPPVec;

	ret = res;
	return ret;
}

double funIK_xyzpp11A(const std::vector<double> &x, std::vector<double> &grad, void *fData) {
	FDATA11A *ptr;
	ptr = (FDATA11A*)fData; // recover pt & tgt

	double qp[5] = { x[0],x[1],x[2],x[3],x[4] };
	Eigen::Vector3d xyzVec;
	ptr->pt.qps2point(qp, &xyzVec); // find current point

	Eigen::Vector2d ppVec;
	ptr->pp.qps2phipsi(qp, &ppVec);

	Eigen::VectorXd tgt(5);
	tgt << ptr->tgtXYZVec, ptr->tgtPPVec;
	Eigen::VectorXd cur(5);
	cur << xyzVec, ppVec;

	return (tgt - cur).norm();
}
InvKNLOpt_xyzpp11A::InvKNLOpt_xyzpp11A(KINEMATICPARAMS11A kParams, NLOPTPARAMS nParams, JOINTLIMITS jLims) {
	kinParams = kParams;
	nlParams = nParams;
	jntLims = jLims;
}
int InvKNLOpt_xyzpp11A::solve(double *qps, double *xyzGoal, double *ppGoal) {
	int ret = -1;
	Point11A pt(kinParams);
	PhiPsi11A pp(kinParams);

	Eigen::VectorXd qpVec(5);
	Eigen::Vector3d curXYZVec;
	Eigen::Vector2d curPPVec;
	Eigen::Vector3d tgtXYZVec;
	Eigen::Vector2d tgtPPVec;
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];
	tgtXYZVec << xyzGoal[0], xyzGoal[1], xyzGoal[2];
	tgtPPVec << ppGoal[0], ppGoal[1];

	FDATA11A ptpptgt;
	ptpptgt.pt = pt;
	ptpptgt.tgtXYZVec = tgtXYZVec;
	ptpptgt.pp = pp;
	ptpptgt.tgtPPVec = tgtPPVec;
	ptpptgt.nlParams = nlParams;
//printf("pp  %f %f\n", ppGoal[0], ppGoal[1]);
//printf("ppt %f %f\n", tgtPPVec(0), tgtPPVec(1));
//printf("ppt %f %f\n", ptpptgt.tgtPPVec(0), ptpptgt.tgtPPVec(1));

	nlopt::opt alg(ikTranslateNLOptAlg(nlParams.method), 5);

	// set objective
	void *ptr;
	ptr = &ptpptgt;
	alg.set_min_objective(funIK_xyzpp11A, ptr);

	// set initial step size (only for nongrad)
	std::vector<double> dx0(5);
	for (int i = 0; i < 5; i++) { dx0[i] = .01; }
	alg.set_initial_step(dx0);

	// set boundary constraints on x
	std::vector<double> limUp(5);
	std::vector<double> limDn(5);
	for (int i = 0; i < 5; i++) {
		limUp[i] = jntLims.up[i];
		limDn[i] = jntLims.dn[i];
	}
	alg.set_upper_bounds(limUp);
	alg.set_lower_bounds(limDn);

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec); //sec
//std::cout << "min: " << alg.get_stopval() << std::endl;
//std::cout << "ftl: " << alg.get_ftol_abs() << std::endl;
////std::cout << "xtl: " << alg.get_xtol_abs() << std::endl;
//std::cout << "mxe: " << alg.get_maxeval() << std::endl;
//std::cout << "mxt: " << alg.get_maxtime() << std::endl;
//std::cout << "mtd: " << alg.get_algorithm_name() << std::endl;
//printf("nlpms: %f %f %f %d %f %d\n", nlParams.minFunVal, nlParams.tolFunAbs, nlParams.tolXAbs, nlParams.maxIts, nlParams.maxTimeSec, nlParams.method);
//std::cout << nlParams.maxIts << "  " << nlParams.method << std::endl;

	//may need to set the local optimizer, eg MLSL http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
	if (ikTranslateNLOptAlg(nlParams.method) == GN_MLSL || ikTranslateNLOptAlg(nlParams.method) == GN_MLSL_LDS) {
		nlopt::opt lAlg(nlopt::LN_NELDERMEAD, alg.get_dimension());

		lAlg.set_min_objective(funIK_xyzpp11A, ptr);
		lAlg.set_upper_bounds(limUp);
		lAlg.set_lower_bounds(limDn);
		lAlg.set_stopval(nlParams.minFunVal);
		lAlg.set_ftol_abs(nlParams.tolFunAbs);
		lAlg.set_xtol_abs(nlParams.tolXAbs);
		lAlg.set_maxeval(nlParams.maxIts);
		lAlg.set_maxtime(nlParams.maxTimeSec);

		alg.set_local_optimizer(lAlg);
	}

	// solve
	nlopt::result res;
	double fmin = 1e3;
	std::vector<double> x(5);
	//set initial position within bounds
	for (int i = 0; i < 5; i++) {
		x[i] = qps[i];
		if (limUp[i] > 0 && x[i] > limUp[i])
			x[i] = limUp[i] * .9999;
		else if (limUp[i] < 0 && x[i] > limUp[i])
			x[i] = limUp[i] * 1.0001;
		if (limDn[i] < 0 && x[i] < limDn[i])
			x[i] = limDn[i] * .9999;
		else if (limDn[i] > 0 && x[i] < limDn[i])
			x[i] = limDn[i] * 1.0001;
	}
	
//printf("x : "); for (int i = 0; i < 5; i++) { printf("%f ", x[i]); } printf("\n");
//printf("up: "); for (int i = 0; i < 5; i++) { printf("%f ", limUp[i]); } printf("\n");
//printf("dn: "); for (int i = 0; i < 5; i++) { printf("%f ", limDn[i]); } printf("\n");
//printf("xy: "); for (int i = 0; i < 3; i++) { printf("%f ", ptpptgt.tgtXYZVec(i)); } printf("\n");
//printf("pp: "); for (int i = 0; i < 2; i++) { printf("%f ", ptpptgt.tgtPPVec(i)); } printf("\n");
	try {
		res = alg.optimize(x, fmin);
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
		std::cout << "Caught: " << e.what() << std::endl;
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
//printf("alg: res %d    fmin = %f    x = [%f %f %f %f %f]\n", res, fmin, x[0], x[1], x[2], x[3], x[4]);
	for (int i = 0; i < 5; i++) { qpVec(i) = x[i]; }
	pt.qps2point(&qpVec, &curXYZVec);
	pp.qps2phipsi(&qpVec, &curPPVec);

	// assign outputs
	Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	Eigen::Map<Eigen::VectorXd>(xyzGoal, curXYZVec.rows()) = curXYZVec;
	Eigen::Map<Eigen::VectorXd>(ppGoal, curPPVec.rows()) = curPPVec;

	ret = res;
	return ret;
}
int InvKNLOpt_xyzpp11A::getFval(double *qps, double *xyzGoal, double *ppGoal, double *fval) {
	int ret = -99;
	Point11A pt(kinParams);
	PhiPsi11A pp(kinParams);

	Eigen::VectorXd qpVec(5);
	Eigen::Vector3d curXYZVec;
	Eigen::Vector2d curPPVec;
	Eigen::Vector3d tgtXYZVec;
	Eigen::Vector2d tgtPPVec;
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];
	tgtXYZVec << xyzGoal[0], xyzGoal[1], xyzGoal[2];
	tgtPPVec << ppGoal[0], ppGoal[1];

	FDATA11A ptpptgt;
	ptpptgt.pt = pt;
	ptpptgt.tgtXYZVec = tgtXYZVec;
	ptpptgt.pp = pp;
	ptpptgt.tgtPPVec = tgtPPVec;
	ptpptgt.nlParams = nlParams;

	//set position
	std::vector<double> x(5);
	for (int i = 0; i < 5; i++) { x[i] = qps[i]; }

	//solve
	*fval = funIK_xyzpp11A(x, x, &ptpptgt);
	pt.qps2point(qps, xyzGoal);
	pp.qps2phipsi(qps, ppGoal);
	return ret;
}

// Point Center Sum
double funPointCenterSum(const std::vector<double> &x, std::vector<double> &grad, void *fData) {
	FDATA6A *ptr;
	ptr = (FDATA6A*)fData; // recover pt & tgt

	double qp[5] = { x[0],x[1],x[2],x[3],x[4] };
	Eigen::Vector3d xyzVec;
	Eigen::Vector2d csVec;
	ptr->pt.qps2point(qp, &xyzVec); // find current point
	ptr->cs.qps2ctrSum(qp, &csVec);
	
	double prim = (ptr->tgtXYZVec - xyzVec).norm();
	return prim + csVec.norm() * 1e-5 *ptr->nlParams.minFunVal / prim; // as in newton, increase cs influence while nearing target
}
InvKNLOpt_xyzCenterSum6A::InvKNLOpt_xyzCenterSum6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nParams, JOINTLIMITS jLims) {
	kinParams = kParams;
	nlParams = nParams;
	jntLims = jLims;
}
int InvKNLOpt_xyzCenterSum6A::solve(double *qps, double *xyz) {
	int ret = -1;
	Point6A pt(kinParams);
	CenterSum cs(jntLims);
	
	Eigen::VectorXd qpVec(5);
	Eigen::Vector3d curVec;
	Eigen::Vector3d tgtVec;
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];
	tgtVec << xyz[0], xyz[1], xyz[2];

	FDATA6A pttgt;
	pttgt.pt = pt;
	pttgt.tgtXYZVec = tgtVec;
	pttgt.cs = cs;

	nlopt::opt alg(ikTranslateNLOptAlg(nlParams.method), 5);

	// set objective
	void *ptr;
	ptr = &pttgt;
	alg.set_min_objective(funPointCenterSum, ptr);

	// set initial step size (only for nongrad)
	std::vector<double> dx0(5);
	for (int i = 0; i < 5; i++) { dx0[i] = .01; }
	alg.set_initial_step(dx0);

	// set boundary constraints on x
	std::vector<double> limUp(5);
	std::vector<double> limDn(5);
	for (int i = 0; i < 5; i++) {
		limUp[i] = jntLims.up[i];
		limDn[i] = jntLims.dn[i];
	}
	alg.set_upper_bounds(limUp);
	alg.set_lower_bounds(limDn);

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	// solve
	nlopt::result res;
	double fmin = 1e3;
	std::vector<double> x(5);
	for (int i = 0; i < 5; i++) { x[i] = .1; }
	try {
		res = alg.optimize(x, fmin);
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
		std::cout << "Caught: " << e.what() << std::endl;
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	//printf("alg: res %d    fmin = %f    x = [%f %f %f %f %f]\n", res, fmin, x[0], x[1], x[2], x[3], x[4]);
	for (int i = 0; i < 5; i++) { qpVec(i) = x[i]; }
	pt.qps2point(&qpVec, &curVec);

	// assign outputs
	Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	Eigen::Map<Eigen::VectorXd>(xyz, curVec.rows()) = curVec;

	ret = res;
	return ret;
}

// Point Phi Psi Center Sum
double funPointPhiPsiCenterSum(const std::vector<double> &x, std::vector<double> &grad, void *fData) {
	FDATA6A *ptr;
	ptr = (FDATA6A*)fData; // recover pt & tgt

	double qp[5] = { x[0],x[1],x[2],x[3],x[4] };
	Eigen::Vector3d xyzVec;
	Eigen::Vector2d ppVec;
	Eigen::Vector2d csVec;
	ptr->pt.qps2point(qp, &xyzVec); // find current point
	ptr->pp.qps2phipsi(qp, &ppVec);
	ptr->cs.qps2ctrSum(qp, &csVec);

	double prim = (ptr->tgtXYZVec - xyzVec).norm() + (ptr->tgtPPVec - ppVec).norm();
	return prim + csVec.norm() * 1e-5 *ptr->nlParams.minFunVal / prim; // as in newton, increase cs influence while nearing target
}
InvKNLOpt_xyzppCenterSum6A::InvKNLOpt_xyzppCenterSum6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nParams, JOINTLIMITS jLims) {
	kinParams = kParams;
	nlParams = nParams;
	jntLims = jLims;
}
int InvKNLOpt_xyzppCenterSum6A::solve(double *qps, double *xyzGoal, double *ppGoal) {
	int ret = -1;
	Point6A pt(kinParams);
	PhiPsi6A pp(kinParams);
	CenterSum cs(jntLims);

	Eigen::VectorXd qpVec(5);
	Eigen::Vector3d curXYZVec;
	Eigen::Vector3d tgtXYZVec;
	Eigen::Vector2d curPPVec;
	Eigen::Vector2d tgtPPVec;
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];
	tgtXYZVec << xyzGoal[0], xyzGoal[1], xyzGoal[2];
	tgtPPVec << ppGoal[0], ppGoal[1];

	FDATA6A pttgt;
	pttgt.pt = pt;
	pttgt.tgtXYZVec = tgtXYZVec;
	pttgt.pp = pp;
	pttgt.tgtPPVec = tgtPPVec;
	pttgt.cs = cs;

	nlopt::opt alg(ikTranslateNLOptAlg(nlParams.method), 5); // 5 = length of qpVec

	// set objective
	void *ptr;
	ptr = &pttgt;
	alg.set_min_objective(funPointPhiPsiCenterSum, ptr);

	// set initial step size (only for nongrad)
	std::vector<double> dx0(7);
	for (int i = 0; i < 5; i++) { dx0[i] = .01; }
	alg.set_initial_step(dx0);

	// set boundary constraints on x
	std::vector<double> limUp(5);
	std::vector<double> limDn(5);
	for (int i = 0; i < 5; i++) {
		limUp[i] = jntLims.up[i];
		limDn[i] = jntLims.dn[i];
	}
	alg.set_upper_bounds(limUp);
	alg.set_lower_bounds(limDn);

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	// solve
	nlopt::result res;
	double fmin = 1e3;
	std::vector<double> x(5);
	for (int i = 0; i < 5; i++) { x[i] = .1; }
	try {
		res = alg.optimize(x, fmin);
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
		std::cout << "Caught: " << e.what() << std::endl;
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	//printf("alg: res %d    fmin = %f    x = [%f %f %f %f %f]\n", res, fmin, x[0], x[1], x[2], x[3], x[4]);
	for (int i = 0; i < 5; i++) { qpVec(i) = x[i]; }
	pt.qps2point(&qpVec, &curXYZVec);

	// assign outputs
	Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	Eigen::Map<Eigen::VectorXd>(xyzGoal, curXYZVec.rows()) = curXYZVec;
	Eigen::Map<Eigen::VectorXd>(ppGoal, curPPVec.rows()) = curPPVec;
	
	ret = res;
	return ret;
}

// Single Point, minimize change from last
double funPointMinChange(const std::vector<double> &x, std::vector<double> &grad, void *fData) {
	FDATA6A *ptr;
	ptr = (FDATA6A*)fData; // recover pt & tgt

	Eigen::VectorXd qpsVec(5);
	qpsVec << x[0], x[1], x[2], x[3], x[4];
	Eigen::Vector3d xyzVec;

	ptr->pt.qps2point(&qpsVec, &xyzVec); // find current point
	
	//return (ptr->tgtXYZVec - xyzVec).norm() + (ptr->qpLastVec - qpsVec).norm(); // find error norm
	//return (ptr->tgtXYZVec - xyzVec).norm() + (ptr->qpLastVec - qpsVec).head(4).norm(); // don't care how the tip2tgt distance changes
	return (ptr->tgtXYZVec - xyzVec).norm() + (ptr->qpLastVec - qpsVec).head(4).norm() / (ptr->tgtXYZVec - xyzVec).norm();
}
InvKNLOpt_xyzMinChange6A::InvKNLOpt_xyzMinChange6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nParams, JOINTLIMITS jLims) {
	kinParams = kParams;
	nlParams = nParams;
	jntLims = jLims;
}
int InvKNLOpt_xyzMinChange6A::solve(double *qps, double *xyz, double *qpsLast) {
	int ret = -1;
	Point6A pt(kinParams);

	Eigen::VectorXd qpVec(5);
	Eigen::VectorXd qpLastVec(5);
	Eigen::Vector3d curVec;
	Eigen::Vector3d tgtVec;
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];
	qpLastVec << qpsLast[0], qpsLast[1], qpsLast[2], qpsLast[3], qpsLast[4];
	tgtVec << xyz[0], xyz[1], xyz[2];

	FDATA6A pttgt;
	pttgt.pt = pt;
	pttgt.tgtXYZVec = tgtVec;
	pttgt.qpLastVec = qpLastVec;

	nlopt::opt alg(ikTranslateNLOptAlg(nlParams.method), 5);

	// set objective
	void *ptr;
	ptr = &pttgt;
	alg.set_min_objective(funPointMinChange, ptr);

	// set initial step size (only for nongrad)
	std::vector<double> dx0(5);
	for (int i = 0; i < 5; i++) { dx0[i] = .01; }
	alg.set_initial_step(dx0);

	// set boundary constraints on x
	std::vector<double> limUp(5);
	std::vector<double> limDn(5);
	for (int i = 0; i < 5; i++) {
		limUp[i] = jntLims.up[i];
		limDn[i] = jntLims.dn[i];
	}
	alg.set_upper_bounds(limUp);
	alg.set_lower_bounds(limDn);

	// set start joint angles to be within bounds
	std::vector<double> x(5);
	for (int i = 0; i < 5; i++) {
		x[i] = qpVec(i);
		if (qpVec(i) <= limDn[i]) {
			x[i] = limDn[i];
		}
		if (limUp[i] <= qpVec(i)) {
			x[i] = limUp[i];
		}
	}

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	// solve
	nlopt::result res;
	double fmin = 1e3;
	try {
		res = alg.optimize(x, fmin);
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
		printf("Caught: %s for x=[ %4.4f %4.4f %4.4f %4.4f %4.4f ]\n", e.what(), x[0], x[1], x[2], x[3], x[4]);
		for (int i = 0; i < 5; i++) { printf("%d: %+08.3f < [%+08.3f] < %+08.3f\n", i, limDn[i], x[i], limUp[i]); }
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}
	//printf("alg: res %d    fmin = %f    x = [%f %f %f %f %f]\n", res, fmin, x[0], x[1], x[2], x[3], x[4]);
	for (int i = 0; i < 5; i++) { qpVec(i) = x[i]; }
	pt.qps2point(&qpVec, &curVec);

	// assign outputs
	Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	Eigen::Map<Eigen::VectorXd>(xyz, curVec.rows()) = curVec;

	ret = res;
	return ret;
}

// Single Point, minimize change from last, kalman filter in cost function
double funPoint_kalmanX(const std::vector<double> &x, std::vector<double> &grad, void *fData) {
	FDATA6A *ptr;
	ptr = (FDATA6A*)fData; // recover pt & tgt
	Eigen::VectorXd qpsVec(5);
	Eigen::Vector3d xyzVec;
	double xmea, xhat, vhat = 0;

	//filter
	for (int i = 0; i < 5; i++) {
		xmea = x[i];
		xhat = (ptr->qpLastVec(i));
		vhat = (ptr->qdLastVec(i));
		ptr->kl.filterNewMeasure( &xmea, &xhat, &vhat, (ptr->K11s)+i, (ptr->K21s)+i);
		qpsVec(i) = xhat; //only output for solver
	}

	ptr->pt.qps2point(&qpsVec, &xyzVec); // find location of x1 = filtered(x*)
	
	//return (ptr->tgtXYZVec - xyzVec).norm();// oscillates/magnifies sensor noise
	//return (ptr->tgtXYZVec - xyzVec).norm() + (ptr->qpLastVec - qpsVec).norm(); //similarity limits convergence
	return (ptr->tgtXYZVec - xyzVec).norm() + (ptr->qpLastVec - qpsVec).head(4).norm(); // don't care how the tip2tgt distance changes
}
InvKNLOpt_xyzKalX6A::InvKNLOpt_xyzKalX6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nParams, JOINTLIMITS jLims, double tsSeconds) {
	kinParams = kParams;
	nlParams = nParams;
	jntLims = jLims;
	tsSec = tsSeconds;
}
int InvKNLOpt_xyzKalX6A::solve(double *qps, double *xyz, double *qpsLast, double *qdsLast, double *K11s, double *K21s) {
	int ret = -1;
	Point6A pt(kinParams);
	Kalman2State kl(tsSec);

	Eigen::VectorXd qpVec(5);
	Eigen::VectorXd qpLastVec(5);
	Eigen::VectorXd qdLastVec(5);
	Eigen::Vector3d curVec;
	Eigen::Vector3d tgtVec;
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];
	qpLastVec << qpsLast[0], qpsLast[1], qpsLast[2], qpsLast[3], qpsLast[4];
	qdLastVec << qdsLast[0], qdsLast[1], qdsLast[2], qdsLast[3], qdsLast[4];
	tgtVec << xyz[0], xyz[1], xyz[2];

	FDATA6A pttgt;
	pttgt.pt = pt;
	pttgt.tgtXYZVec = tgtVec;
	pttgt.qpLastVec = qpLastVec;
	pttgt.qdLastVec = qdLastVec;
	pttgt.K11s = K11s;
	pttgt.K21s = K21s;
	pttgt.kl = kl;

	nlopt::opt alg(ikTranslateNLOptAlg(nlParams.method), 5);
	
	// set objective
	void *ptr;
	ptr = &pttgt;
	alg.set_min_objective(funPoint_kalmanX, ptr);
	
	// set initial step size (only used by nongradient methods)
	std::vector<double> dx0(5);
	for (int i = 0; i < 4; i++) { dx0[i] = .01; }
	dx0[4] = 1;
	alg.set_initial_step(dx0);

	//set boundary constraints
	std::vector<double> limUp(5);
	std::vector<double> limDn(5);
	for (int i = 0; i < 5; i++) {
		limUp[i] = jntLims.up[i];
		limDn[i] = jntLims.dn[i];
	}
	alg.set_upper_bounds(limUp);
	alg.set_lower_bounds(limDn);

	// set start joint angles to be within bounds
	std::vector<double> x(5);
	for (int i = 0; i < 5; i++) {
		x[i] = qpVec(i);
		if (qpVec(i) <= limDn[i]) {
			x[i] = limDn[i];
		}
		if (limUp[i] <= qpVec(i)) {
			x[i] = limUp[i];
		}
	}

	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	// solve
	nlopt::result res;
	double fmin = 1e3;
	try {
		res = alg.optimize(x, fmin);
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
		printf("Caught: %s for x=[ %4.4f %4.4f %4.4f %4.4f %4.4f ]\n", e.what(), x[0], x[1], x[2], x[3], x[4]);
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}

	//rerun the filter on the last x to find the correct qps and kalman state
	double xmea, xhat, vhat = 0;
	for (int i = 0; i < 5; i++) {
		xmea = x[i];
		xhat = qpLastVec(i);
		vhat = qdLastVec(i);
		kl.filterNewMeasure(&xmea, &xhat, &vhat, K11s+i, K21s+i);
		qpVec(i) = xhat; //output filtered qps
		qpsLast[i] = xhat; //output filter state
		qdsLast[i] = vhat;
	}
	//double res = *xmea - *xhat - *vhat * ts;
	//*xhat = *xhat + ts * *vhat + *K11 * res;
	//*vhat = *vhat + *K21*res;
	//double r = sqrt( pow(x[0] - qpLastVec(0) - qdLastVec(0)*tsSec, 2) //residual as |qerr| ?
	//	           + pow(x[1] - qpLastVec(1) - qdLastVec(1)*tsSec, 2)
	//	           + pow(x[2] - qpLastVec(2) - qdLastVec(2)*tsSec, 2)
	//               + pow(x[3] - qpLastVec(3) - qdLastVec(3)*tsSec, 2)
	//	           + pow(x[4] - qpLastVec(4) - qdLastVec(4)*tsSec, 2) );
	//for (int i = 0; i < 5; i++) {
	//	qpVec(i) = qpLastVec(i) + qdLastVec(i)*tsSec + K11s[i] * r;
	//	qpsLast[i] = qpVec(i);
	//	qdsLast[i] = qdLastVec(i) + K21s[i] * r;
	//}
	pt.qps2point(&qpVec, &curVec); //find converged xyz

	// assign outputs from vectors to arrays
	Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	Eigen::Map<Eigen::VectorXd>(xyz, curVec.rows()) = curVec;

	ret = res;
	return ret;
}

// Single Point, bound x by stdev*K11, filter solution
double funPoint_kalmanStdBounds(const std::vector<double> &x, std::vector<double> &grad, void *fData) {
	FDATA6A *ptr;
	ptr = (FDATA6A*)fData; // recover pt & tgt
	Eigen::VectorXd qpsVec(5);
	Eigen::Vector3d xyzVec;

	for (int i = 0; i < 5; i++) { qpsVec(i) = x[i]; }
	//qpsVec(4) = x[4];
	//double xmea, xhat, vhat, K11, K21 = 0;
	//for (int i = 0; i < 4; i++) {
	//	xmea = x[i];
	//	xhat = ptr->qpLastVec(i);
	//	vhat = ptr->qdLastVec(i);
	//	K11 = ptr->K11s[i];
	//	K21 = ptr->K21s[i];
	//	ptr->kl.filterNewMeasure(&xmea, &xhat, &vhat, &K11, &K21);
	//	qpsVec(i) = xhat;
	//}
	
	ptr->pt.qps2point(&qpsVec, &xyzVec); // find location of x1 = filtered(x*)

	//return (ptr->tgtXYZVec - xyzVec).norm();
	//return (ptr->tgtXYZVec - xyzVec).norm() + (ptr->qpLastVec - qpsVec).norm();
	return (ptr->tgtXYZVec - xyzVec).norm() + (ptr->qpLastVec - qpsVec).head(4).norm(); // don't care how the tip2tgt distance changes
}
InvKNLOpt_xyzKalBounds6A::InvKNLOpt_xyzKalBounds6A(KINEMATICPARAMS6A kParams, NLOPTPARAMS nParams, JOINTLIMITS jLims, double tsSeconds) {
	kinParams = kParams;
	nlParams = nParams;
	jntLims = jLims;
	tsSec = tsSeconds;
}
int InvKNLOpt_xyzKalBounds6A::solve(double *qps, double *xyz, double *qpsLast, double *qdsLast, double *K11s, double *K21s, double *stds) {
	int ret = -1;
	Point6A pt(kinParams);
	Kalman2State kl(tsSec);

	Eigen::VectorXd qpVec(5);
	Eigen::VectorXd qpLastVec(5);
	Eigen::VectorXd qdLastVec(5);
	Eigen::Vector3d curVec;
	Eigen::Vector3d tgtVec;
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];
	qpLastVec << qpsLast[0], qpsLast[1], qpsLast[2], qpsLast[3], qpsLast[4];
	qdLastVec << qdsLast[0], qdsLast[1], qdsLast[2], qdsLast[3], qdsLast[4];
	tgtVec << xyz[0], xyz[1], xyz[2];

	FDATA6A pttgt;
	pttgt.pt = pt;
	pttgt.tgtXYZVec = tgtVec;
	pttgt.qpLastVec = qpLastVec;
	pttgt.qdLastVec = qdLastVec;
	pttgt.K11s = K11s;
	pttgt.K21s = K21s;
	pttgt.kl = kl;

	nlopt::opt alg(ikTranslateNLOptAlg(nlParams.method), 5);

	// set objective
	void *ptr;
	ptr = &pttgt;
	alg.set_min_objective(funPoint_kalmanStdBounds, ptr);

	// set initial step size (only used by nongradient methods)
	std::vector<double> dx0(5);
	for (int i = 0; i < 5; i++) { dx0[i] = .001; }
	alg.set_initial_step(dx0);

	//set boundary constraints based on qpLast and stdev
	std::vector<double> limUp(5);
	std::vector<double> limDn(5);
	limUp[4] = jntLims.up[4]; // don't need to modify the tip2tgt bounds
	limDn[4] = jntLims.dn[4];
	double pos,neg, xhat, vhat = 0;
	for (int i = 0; i < 4; i++) {
		xhat = qpLastVec(i);
		vhat = qdLastVec(i);
		pos = xhat + vhat*tsSec + K11s[i] * stds[i];
		neg = xhat + vhat*tsSec - K11s[i] * stds[i];
		if (neg > pos) { //switch'em
			neg = xhat + vhat*tsSec + K11s[i] * stds[i];
			pos = xhat + vhat*tsSec - K11s[i] * stds[i];
		}
		if (pos < jntLims.up[i]  && jntLims.dn[i] < pos) {
			limUp[i] = pos;
		}else {
			limUp[i] = jntLims.up[i];
		}
		if (jntLims.dn[i] < neg && neg < jntLims.up[i]) {
			limDn[i] = neg;
		}else {
			limDn[i] = jntLims.dn[i];
		}
		if (limUp[i] - limDn[i] < 0.2) { // if the bounds are too narrow, revert, relying on the post-solver filtering to smooth
			limUp[i] = jntLims.up[i];
			limDn[i] = jntLims.dn[i];
		}
		//printf("%d: %+08.6f < %+08.6f < [%+08.6f +- %08.6f*%03.6f] < %+08.6f < %+08.6f = %+08.6f\n",i, jntLims.dn[i], neg, xhat+vhat*tsSec, K11s[i], stds[i], pos, jntLims.up[i], limUp[i]-limDn[i]);
	}
	alg.set_upper_bounds(limUp);
	alg.set_lower_bounds(limDn);

	// set start joint angles to be within bounds
	std::vector<double> x(5);
	for (int i = 0; i < 5; i++) {
		x[i] = qpVec(i);
		if (qpVec(i) <= limDn[i]) {
			x[i] = limDn[i];
		}
		if (limUp[i] <= qpVec(i)) {
			x[i] = limUp[i];
		}
	}
	
	// set stopping criteria
	alg.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	alg.set_ftol_abs(nlParams.tolFunAbs);
	alg.set_xtol_abs(nlParams.tolXAbs);
	alg.set_maxeval(nlParams.maxIts);
	alg.set_maxtime(nlParams.maxTimeSec);

	// solve
	nlopt::result res;
	double fmin = 1e3;
	try {
		res = alg.optimize(x, fmin);
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
		printf("Caught: %s for x=[ %4.4f %4.4f %4.4f %4.4f %4.4f ]\n", e.what(), x[0], x[1], x[2], x[3], x[4]);
		for (int i = 0; i < 5;  i++){ printf("%d: %+08.3f < [%+08.3f] < %+08.3f\n", i, limDn[i], x[i], limUp[i]); }
	}
	catch (std::bad_alloc e) {
		res = nlopt::FAILURE;
		std::cout << "Caught: " << e.what() << std::endl;
	}

	// unpack from solver & filter
	double xmea = 0; xhat = 0; vhat = 0;
	for (int i = 0; i < 5; i++) { 
		// plain
		//qpVec(i) = x[i];
		//qpsLast[i] = x[i];
		// filtered
		xmea = x[i];
		xhat = qpsLast[i];
		vhat = qdsLast[i];
		kl.filterNewMeasure(&xmea, &xhat, &vhat, K11s + i, K21s + i);
		qpVec(i) = xhat;
		qps[i] = xhat;
		qpsLast[i] = xhat;
		qdsLast[i] = vhat;
	}
	//qps[4] = x[4];
	//qpsLast[4] = x[4];
	//qpVec(4) = x[4];
	pt.qps2point(&qpVec, &curVec); //find converged xyz

	// assign outputs from vectors to arrays
	//Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	Eigen::Map<Eigen::VectorXd>(xyz, curVec.rows()) = curVec;

	ret = res;
	return ret;
}

*/


nlopt::algorithm ikTranslateNLOptAlg(nlMethod method) {
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
	default: return nlopt::GN_DIRECT;
	}
}

//useless definition to help linker, one for every possible template value
template class InvK_nlopt<TaskXYZ<FwdK6A>, FwdK6A>;
template class InvK_nlopt<TaskXYZ<FwdK11A>, FwdK11A>;
template class InvK_nlopt<TaskXYZUxUyUz<FwdK6A>, FwdK6A>;
template class InvK_nlopt<TaskXYZUxUyUz<FwdK11A>, FwdK11A>;