
#include "ik_newton.h"
#include <iostream>

InvKNewtonPoint::InvKNewtonPoint(KINEMATICPARAMS6A kParams, NEWTONPARAMS nParams, JOINTLIMITS jLims) {
	kinParams = kParams;
	nrParams = nParams;
	jntLims = jLims;
}
int InvKNewtonPoint::solve(double *qps, double *xyz) {
	int ret = -1;
	Point6A pt(kinParams);
	
	Eigen::VectorXd qpVec(5);
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];

	//printf("i qps[%f %f %f %f %f]\n", qpVec(0), qpVec(1), qpVec(2), qpVec(3), qpVec(4));

	Eigen::Vector3d err;
	Eigen::Vector3d cur;
	Eigen::Vector3d tgt;
	err.setZero();
	tgt << xyz[0], xyz[1], xyz[2];

	Eigen::MatrixXd J(3, 5);
	Eigen::MatrixXd JP(5, 3);
	J.setZero();

	//for (int it = 0; it < 5; it++){
	for (int it = 0; it < nrParams.maxIts; it++){

		pt.qps2point(&qpVec, &cur);
		err = tgt - cur;
		//printf("%d tgt [%f %f %f] - cur [%f %f %f] = err[%f %f %f] = %f\n", it, tgt(0), tgt(1), tgt(2), cur(0), cur(1), cur(2), err(0), err(1), err(2), err.norm());
		//printf("%d qps[%f %f %f %f %f] -- %f\n", it, qpVec(0), qpVec(1), qpVec(2), qpVec(3), qpVec(4), err.norm());

		pt.Jacobian(&J, &qpVec, &(nrParams.epslon) );
		JP = pseudoInverse(J);

		qpVec += JP * err * nrParams.fact01;

		// cath bends one way
		if (qpVec(3) < jntLims.dn[3]) {
			qpVec(3) = jntLims.dn[3] * 2;
		}

		// mod for readability
		qpVec(0) = fmod(qpVec(0), 2 * M_PI);
		qpVec(1) = fmod(qpVec(1), 2 * M_PI);
		qpVec(2) = fmod(qpVec(2), 2 * M_PI);
		qpVec(3) = fmod(qpVec(3), 2 * M_PI);
	}

	// assign outputs
	Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	Eigen::Map<Eigen::VectorXd>(xyz, cur.rows()) = cur;

	ret = (nrParams.errTol > (tgt - cur).norm())*2-1; // -1 = failed, 1 = success

	return ret;
}

InvKNewtonPriorityPointPhiPsi::InvKNewtonPriorityPointPhiPsi(KINEMATICPARAMS6A kParams, NEWTONPARAMS nParams, JOINTLIMITS jLims) {
	kinParams = kParams;
	nrParams = nParams;
	jntLims = jLims;
}
int InvKNewtonPriorityPointPhiPsi::solve(double *qps, double *xyzGoal, double *ppGoal) {
	int ret = -1;
	Point6A pt(kinParams);
	PhiPsi6A pp(kinParams);

	Eigen::VectorXd qpVec(5);
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];
	//printf(".....i qps[%f %f %f %f %f]\n", qpVec(0), qpVec(1), qpVec(2), qpVec(3), qpVec(4));
	

	// primary task: xyz position
	Eigen::Vector3d ptErr;
	Eigen::Vector3d ptCur;
	Eigen::Vector3d ptTgt;
	ptErr.setZero();
	ptTgt << xyzGoal[0], xyzGoal[1], xyzGoal[2];
	Eigen::MatrixXd ptJ(3, 5);
	Eigen::MatrixXd ptJP(5, 3);
	ptJ.setZero();
	ptJP.setZero();

	// secondary task: phipsi orientation
	Eigen::Vector2d ppErr;
	Eigen::Vector2d ppCur;
	Eigen::Vector2d ppTgt;
	ppErr.setZero();
	ppTgt << ppGoal[0], ppGoal[1];
	Eigen::MatrixXd ppJ(2, 5);
	Eigen::MatrixXd ppJP(5, 2);
	ptJ.setZero();
	ptJP.setZero();
	
	Eigen::MatrixXd eye(5,5); eye.setIdentity();

	//for (int it = 0; it < 5; it++){
	for (int it = 0; it < nrParams.maxIts; it++) {

		pt.qps2point(&qpVec, &ptCur);
		ptErr = ptTgt - ptCur;
		pp.qps2phipsi(&qpVec, &ppCur);
		ppErr = ppTgt - ppCur;
		//printf(".....%d tgt[%f %f %f] - cur[%f %f %f] = err[%f %f %f] = %f\n", it, ptTgt(0), ptTgt(1), ptTgt(2), ptCur(0), ptCur(1), ptCur(2), ptErr(0), ptErr(1), ptErr(2), ptErr.norm());
		//printf(".....%d ptCur[%f %f %f] = %f ppCur[%f %f] = %f\n", it, ptCur(0), ptCur(1), ptCur(2), ptErr.norm(), ppCur(0), ppCur(1), ppErr.norm());
		//printf(".....%d qps[%f %f %f %f %f] -- %f -- %f\n", it, qpVec(0), qpVec(1), qpVec(2), qpVec(3), qpVec(4), ptErr.norm(), ppErr.norm());
		pt.Jacobian(&ptJ, &qpVec, &nrParams.epslon);
		ptJP = pseudoInverse(ptJ);
		pp.Jacobian(&ppJ, &qpVec, &nrParams.epslon);
		ppJP = pseudoInverse(ppJ);

		//qpVec += ptJP * ptErr * factor;
		//qpVec += nrParams.factor*(ptJP*ptErr + (eye - ptJP*ptJ)*ppJP*ppErr); //eq 11.44 - 2008_ChiaveriniOrioloWalker - drops minimum velocity on 2nd, tracks components that don't conflict with primary
		qpVec += ptJP*ptErr*nrParams.fact01 + (eye - ptJP*ptJ)*ppJP*ppErr*nrParams.fact02;

		// cath bends one way
		if (qpVec(3) < jntLims.dn[3]) {
			qpVec(3) = jntLims.dn[3]*2;
		}

		// mod for readability
		qpVec(0) = fmod(qpVec(0), 2 * M_PI);
		qpVec(1) = fmod(qpVec(1), 2 * M_PI);
		qpVec(2) = fmod(qpVec(2), 2 * M_PI);
		qpVec(3) = fmod(qpVec(3), 2 * M_PI);
	}

	// assign outputs
	Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	Eigen::Map<Eigen::VectorXd>(xyzGoal, ptCur.rows()) = ptCur;
	Eigen::Map<Eigen::VectorXd>(ppGoal, ppCur.rows()) = ppCur;
	
	ret = (nrParams.errTol > (ptTgt - ptCur).norm()) * 2 - 1; // -1 = failed, 1 = success
	return ret;
}

InvKNewtonPriorityPointCenterSum::InvKNewtonPriorityPointCenterSum(KINEMATICPARAMS6A kParams, NEWTONPARAMS nParams, JOINTLIMITS jLims) {
	kinParams = kParams;
	nrParams = nParams;
	jntLims = jLims;
}
int InvKNewtonPriorityPointCenterSum::solve(double *qps, double *xyz) {
	int ret = -1;
	Point6A pt(kinParams);
	CenterSum cs(jntLims);

	Eigen::VectorXd qpVec(5);
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];

	//printf("i qps[%f %f %f %f %f]\n", qpVec(0), qpVec(1), qpVec(2), qpVec(3), qpVec(4));

	Eigen::Vector3d errVec;
	Eigen::Vector3d curVec;
	Eigen::Vector3d tgtVec;
	Eigen::Vector2d csVec;
	errVec.setZero();
	tgtVec << xyz[0], xyz[1], xyz[2];

	Eigen::MatrixXd ptJ(3, 5);
	Eigen::MatrixXd ptJP(5, 3);
	Eigen::MatrixXd csJ(2, 5);
	Eigen::MatrixXd csJP(5, 2);
	Eigen::MatrixXd eye(5, 5);

	ptJ.setZero();
	eye.setIdentity();

	for (int it = 0; it < nrParams.maxIts; it++) {

		pt.qps2point(&qpVec, &curVec);
		errVec = tgtVec - curVec;
		cs.qps2ctrSum(&qpVec, &csVec);

		//printf("%d tgt [%f %f %f] - cur [%f %f %f] = err[%f %f %f] = %f\n", it, tgt(0), tgt(1), tgt(2), cur(0), cur(1), cur(2), err(0), err(1), err(2), err.norm());
		//printf("%d qps[%f %f %f %f %f] -- %f %f %f = %f -- %f %f\n", it, qpVec(0), qpVec(1), qpVec(2), qpVec(3), qpVec(4), errVec(0),errVec(1),errVec(2), errVec.norm(), csVec(0),csVec(1));

		pt.Jacobian(&ptJ, &qpVec, &(nrParams.epslon));
		ptJP = pseudoInverse(ptJ);
		cs.Jacobian(&csJ, &qpVec, &(nrParams.epslon));
		csJP = pseudoInverse(csJ);
		
		csVec *= nrParams.errTol/errVec.norm(); // if task 1 converged, increase task 2
		//qpVec += ptJP * errVec * nrParams.fact01 + pseudoInverse(csJ*(eye - ptJP*ptJ)) * (csVec*nrParams.fact02 - csJ*ptJP*errVec*nrParams.fact01); //eq 11.43 - 2008_ChiaveriniOrioloWalker
		qpVec += ptJP * errVec * nrParams.fact01 + (eye - ptJP*ptJ)*csJP*csVec*nrParams.fact02; //like eq 11.44 - 2008_ChiaveriniOrioloWalker -- drops min velocity on 2nd task
		

		// cath bends one way
		if (qpVec(3) < jntLims.dn[3]) {
			qpVec(3) = jntLims.dn[3] * 2;
		}

		// mod for readability
		qpVec(0) = fmod(qpVec(0), 2 * M_PI);
		qpVec(1) = fmod(qpVec(1), 2 * M_PI);
		qpVec(2) = fmod(qpVec(2), 2 * M_PI);
		qpVec(3) = fmod(qpVec(3), 2 * M_PI);
	}

	// assign outputs
	Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	Eigen::Map<Eigen::VectorXd>(xyz, curVec.rows()) = curVec;

	ret = (nrParams.errTol > (tgtVec - curVec).norm()) * 2 - 1; // -1 = failed, 1 = success

	return ret;
}

// Like the above, but returns as soon as err1 < errTol
InvKNewtonPriorityPointCenterSumStops::InvKNewtonPriorityPointCenterSumStops(KINEMATICPARAMS6A kParams, NEWTONPARAMS nParams, JOINTLIMITS jLims) {
	kinParams = kParams;
	nrParams = nParams;
	jntLims = jLims;
}
int InvKNewtonPriorityPointCenterSumStops::solve(double *qps, double *xyz) {
	int ret = -1;
	Point6A pt(kinParams);
	CenterSum cs(jntLims);

	Eigen::VectorXd qpVec(5);
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];

	//printf("i qps[%f %f %f %f %f]\n", qpVec(0), qpVec(1), qpVec(2), qpVec(3), qpVec(4));

	Eigen::Vector3d errVec;
	Eigen::Vector3d curVec;
	Eigen::Vector3d tgtVec;
	Eigen::Vector2d csVec;
	errVec.setZero();
	tgtVec << xyz[0], xyz[1], xyz[2];

	Eigen::MatrixXd ptJ(3, 5);
	Eigen::MatrixXd ptJP(5, 3);
	Eigen::MatrixXd csJ(2, 5);
	Eigen::MatrixXd csJP(5, 2);
	Eigen::MatrixXd eye(5, 5);

	ptJ.setZero();
	eye.setIdentity();

	for (int it = 0; it < nrParams.maxIts; it++) {

		pt.qps2point(&qpVec, &curVec);
		errVec = tgtVec - curVec;
		cs.qps2ctrSum(&qpVec, &csVec);

		if (errVec.norm() < nrParams.errTol) { break; }

		//printf("%d tgt [%f %f %f] - cur [%f %f %f] = err[%f %f %f] = %f\n", it, tgt(0), tgt(1), tgt(2), cur(0), cur(1), cur(2), err(0), err(1), err(2), err.norm());
		//printf("%d qps[%f %f %f %f %f] -- %f %f %f = %f -- %f %f\n", it, qpVec(0), qpVec(1), qpVec(2), qpVec(3), qpVec(4), errVec(0),errVec(1),errVec(2), errVec.norm(), csVec(0),csVec(1));

		pt.Jacobian(&ptJ, &qpVec, &(nrParams.epslon));
		ptJP = pseudoInverse(ptJ);
		cs.Jacobian(&csJ, &qpVec, &(nrParams.epslon));
		csJP = pseudoInverse(csJ);

		//csVec *= nrParams.errTol / errVec.norm(); // if task 1 converged, increase task 2
		//qpVec += ptJP * errVec * nrParams.fact01 + pseudoInverse(csJ*(eye - ptJP*ptJ)) * (csVec*nrParams.fact02 - csJ*ptJP*errVec*nrParams.fact01); //eq 11.43 - 2008_ChiaveriniOrioloWalker
		qpVec += ptJP * errVec * nrParams.fact01 + (eye - ptJP*ptJ)*csJP*csVec*nrParams.fact02; //like eq 11.44 - 2008_ChiaveriniOrioloWalker -- drops min velocity on 2nd task
		
		// cath bends one way
		if (qpVec(3) < jntLims.dn[3]) {
			qpVec(3) = jntLims.dn[3] * 2;
		}

		// mod for readability
		qpVec(0) = fmod(qpVec(0), 2 * M_PI);
		qpVec(1) = fmod(qpVec(1), 2 * M_PI);
		qpVec(2) = fmod(qpVec(2), 2 * M_PI);
		qpVec(3) = fmod(qpVec(3), 2 * M_PI);
	}

	// assign outputs
	Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	Eigen::Map<Eigen::VectorXd>(xyz, curVec.rows()) = curVec;

	ret = (nrParams.errTol > (tgtVec - curVec).norm()) * 2 - 1; // -1 = failed, 1 = success

	return ret;
}


// Like the above, but returns as soon as err1 < errTol
InvKNewtonAugmentedPointCenterSumStops::InvKNewtonAugmentedPointCenterSumStops(KINEMATICPARAMS6A kParams, NEWTONPARAMS nParams, JOINTLIMITS jLims) {
	kinParams = kParams;
	nrParams = nParams;
	jntLims = jLims;
}
int InvKNewtonAugmentedPointCenterSumStops::solve(double *qps, double *xyz) {
	int ret = -1;
	PointCenterSum6A pcs(kinParams, jntLims);

	Eigen::VectorXd qpVec(5);
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];

	printf("i qps[%f %f %f %f %f]\n", qpVec(0), qpVec(1), qpVec(2), qpVec(3), qpVec(4));

	Eigen::VectorXd errVec(5);
	Eigen::VectorXd curVec(5);
	Eigen::VectorXd tgtVec(5);
	errVec.setZero();
	//pcs.qps2ptCtrSum(&qpVec, &tgtVec ); // qpVec is (typically) the previous solution, so conservative?
	//tgtVec *= .8;
	tgtVec[3] = 2; tgtVec[4] = 2;
	tgtVec[0] = xyz[0];
	tgtVec[1] = xyz[1];
	tgtVec[2] = xyz[2];
	std::cout << "tgtVec: " << tgtVec.transpose() << "\n";

	Eigen::MatrixXd J(5, 5);
	Eigen::MatrixXd JP(5, 5);
	
	J.setZero();

	for (int it = 0; it < nrParams.maxIts; it++) {
		pcs.qps2ptCtrSum(&qpVec, &curVec);
		errVec << tgtVec - curVec;
		
		if (errVec.norm() < nrParams.errTol) { break; }

		printf("%d qps[%8.3f %8.3f %8.3f %8.3f %8.3f] -- %8.3f %8.3f %8.3f %8.3f %8.3f = %8.3f\n", it, qpVec(0), qpVec(1), qpVec(2), qpVec(3), qpVec(4), errVec(0), errVec(1), errVec(2), errVec(3),errVec(4), errVec.norm());

		pcs.Jacobian(&J, &qpVec, &(nrParams.epslon));
		JP = pseudoInverse(J);
		//qpVec += J.inverse() * errVec * nrParams.fact01;
		//qpVec += JP * errVec * nrParams.fact01;
		qpVec += J.transpose() * errVec * nrParams.fact01;
		
		// cath bends one way
		if (qpVec(3) < jntLims.dn[3]) {
			qpVec(3) = jntLims.dn[3] * 2;
		}

		// mod for readability
		qpVec(0) = fmod(qpVec(0), 2 * M_PI);
		qpVec(1) = fmod(qpVec(1), 2 * M_PI);
		qpVec(2) = fmod(qpVec(2), 2 * M_PI);
		qpVec(3) = fmod(qpVec(3), 2 * M_PI);
	}

	// assign outputs
	Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	xyz[0] = curVec[0]; xyz[1] = curVec[1]; xyz[2] = curVec[2];

	ret = (nrParams.errTol > (tgtVec - curVec).norm()) * 2 - 1; // -1 = failed, 1 = success

	return ret;
}


InvKNewtonAugmentedPointPhiPsiStops::InvKNewtonAugmentedPointPhiPsiStops(KINEMATICPARAMS6A kParams, NEWTONPARAMS nParams, JOINTLIMITS jLims) {
	kinParams = kParams;
	nrParams = nParams;
	jntLims = jLims;
}
int InvKNewtonAugmentedPointPhiPsiStops::solve(double *qps, double *xyz, double *pp) {
	int ret = -1;
	PointPhiPsi6A ppp(kinParams, jntLims);

	Eigen::VectorXd qpVec(5);
	qpVec << qps[0], qps[1], qps[2], qps[3], qps[4];

	//printf("i qps[%f %f %f %f %f]\n", qpVec(0), qpVec(1), qpVec(2), qpVec(3), qpVec(4));

	Eigen::VectorXd errVec(5);
	Eigen::VectorXd curVec(5);
	Eigen::VectorXd tgtVec(5);
	errVec.setZero();
	tgtVec[0] = xyz[0];
	tgtVec[1] = xyz[1];
	tgtVec[2] = xyz[2];
	tgtVec[3] = pp[0];
	tgtVec[4] = pp[1];
	//std::cout << "tgtVec: " << tgtVec.transpose() << "\n";

	Eigen::MatrixXd J(5, 5);
	Eigen::MatrixXd JP(5, 5);

	J.setZero();

	for (int it = 0; it < nrParams.maxIts; it++) {
		ppp.qps2ptPhiPsi(&qpVec, &curVec);
		errVec << tgtVec - curVec;

		if (errVec.norm() < nrParams.errTol) { break; }

		//printf("%d qps[%f %f %f %f %f] -- %f %f %f %f %f = %f\n", it, qpVec(0), qpVec(1), qpVec(2), qpVec(3), qpVec(4), errVec(0), errVec(1), errVec(2), errVec(3), errVec(4), errVec.norm());

		ppp.Jacobian(&J, &qpVec, &(nrParams.epslon));
		JP = pseudoInverse(J);
		qpVec += J.inverse() * errVec * nrParams.fact01;
		//qpVec += JP * errVec * nrParams.fact01;
		//qpVec += J.transpose() * errVec * nrParams.fact01;

		// cath bends one way
		if (qpVec(3) < jntLims.dn[3]) {
			qpVec(3) = jntLims.dn[3] * 2;
		}

		// mod for readability
		qpVec(0) = fmod(qpVec(0), 2 * M_PI);
		qpVec(1) = fmod(qpVec(1), 2 * M_PI);
		qpVec(2) = fmod(qpVec(2), 2 * M_PI);
		qpVec(3) = fmod(qpVec(3), 2 * M_PI);
	}

	// assign outputs
	Eigen::Map<Eigen::VectorXd>(qps, qpVec.rows()) = qpVec;
	xyz[0] = curVec[0]; xyz[1] = curVec[1]; xyz[2] = curVec[2];
	pp[0] = curVec[3]; pp[1] = curVec[4];

	ret = (nrParams.errTol > (tgtVec - curVec).norm()) * 2 - 1; // -1 = failed, 1 = success

	return ret;
}