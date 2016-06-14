#include "taskDefinitions.h"
#include <iostream>

// Point
Point6A::Point6A(){
	KINEMATICPARAMS6A kinParams;
	FwdK6A fklocal(kinParams);
	fk = fklocal;
}
Point6A::Point6A(KINEMATICPARAMS6A inParams) {
	kinParams = inParams;
	FwdK6A fklocal(kinParams);
	fk = fklocal;
}
void Point6A::qps2point(double *qps, Eigen::Vector3d *xyzVec){
	Eigen::Matrix4d H;
	H.setZero();
	H = fk.qps2H05(qps);
	*xyzVec << H(0, 3), H(1, 3), H(2, 3);
}
void Point6A::qps2point(double *qps, double *xyz) {
	Eigen::Vector3d xyzVec;
	//xyzVec << xyz[0], xyz[1], xyz[2];
	Point6A::qps2point(qps, &xyzVec);
	Eigen::Map<Eigen::VectorXd>(xyz, xyzVec.rows()) = xyzVec;
}
void Point6A::qps2point(Eigen::VectorXd *qpVec, Eigen::Vector3d *xyzVec) {
	double qps[5];
	Eigen::Map<Eigen::VectorXd>(qps, qpVec->rows()) = *qpVec;
	Point6A::qps2point(qps, xyzVec);
}
void Point6A::Jacobian(Eigen::MatrixXd *J, double *qps, double *eps) {
	J->setZero(3,5); //rows:xyz, cols:qps

	Eigen::Vector3d pos; pos.setZero();
	Eigen::Vector3d neg; neg.setZero();
	double qpseps[5] = { qps[0],qps[1],qps[2],qps[3],qps[4] };
	
	for (int i = 0; i < 5; i++) {
		qpseps[i] += *eps;
		Point6A::qps2point(qpseps, &pos);
		qpseps[i] = qps[i];

		qpseps[i] -= *eps;
		Point6A::qps2point(qpseps, &neg);
		qpseps[i] = qps[i];

		J->col(i) << (pos - neg) / (2* (*eps));
	}
	//std::cout << "J=\n" << *J << std::endl;
}
void Point6A::Jacobian(Eigen::MatrixXd *J, Eigen::VectorXd *qpVec, double *eps){
	double qps[5];
	Eigen::Map<Eigen::VectorXd>(qps, qpVec->rows()) = *qpVec;

	Point6A::Jacobian(J, qps, eps);
}

Point11A::Point11A() {
	KINEMATICPARAMS11A kinParams;
	FwdK11A fklocal(kinParams);
	fk = fklocal;
}
Point11A::Point11A(KINEMATICPARAMS11A inParams) {
	kinParams = inParams;
	FwdK11A fklocal(kinParams);
	fk = fklocal;
}
void Point11A::qps2point(double *qps, Eigen::Vector3d *xyzVec) {
	Eigen::Matrix4d H;
	H.setZero();
	H = fk.qps2H05(qps);
	*xyzVec << H(0, 3), H(1, 3), H(2, 3);
}
void Point11A::qps2point(double *qps, double *xyz) {
	Eigen::Vector3d xyzVec(xyz[0], xyz[1], xyz[2]);
	Point11A::qps2point(qps, &xyzVec);
	Eigen::Map<Eigen::VectorXd>(xyz, xyzVec.rows()) = xyzVec;
}
void Point11A::qps2point(Eigen::VectorXd *qpVec, Eigen::Vector3d *xyzVec) {
	double qps[5];
	Eigen::Map<Eigen::VectorXd>(qps, qpVec->rows()) = *qpVec;
	Point11A::qps2point(qps, xyzVec);
}

// PhiPsi
PhiPsi6A::PhiPsi6A() {
	KINEMATICPARAMS6A kinParams;
	FwdK6A fklocal(kinParams);
	fk = fklocal;
}
PhiPsi6A::PhiPsi6A(KINEMATICPARAMS6A inParams) {
	kinParams = inParams;
	FwdK6A fklocal(kinParams);
	fk = fklocal;
}
void PhiPsi6A::qps2phipsi(double *qps, Eigen::Vector2d *pp) {
	Eigen::Matrix4d H = fk.qps2H05(qps);
	(*pp)(0) = atan2(H(2, 0), sqrt( pow(H(0, 0),2) + pow(H(1, 0),2) ));
	(*pp)(1) = atan2(H(1, 0), H(0, 0));
}
void PhiPsi6A::qps2phipsi(double *qps, double *pp) {
	Eigen::Vector2d ppVec;
	//ppVec << pp[0], pp[1];
	PhiPsi6A::qps2phipsi(qps, &ppVec);
	Eigen::Map<Eigen::VectorXd>(pp, ppVec.rows()) = ppVec;
}
void PhiPsi6A::qps2phipsi(Eigen::VectorXd *qpVec, Eigen::Vector2d *pp) {
	double qps[5];
	Eigen::Map<Eigen::VectorXd>(qps, qpVec->rows()) = *qpVec;

	PhiPsi6A::qps2phipsi(qps, pp);
}
void PhiPsi6A::Jacobian(Eigen::MatrixXd *J, double *qps, double *eps){
	J->setZero(2,5); //rows:phi,psi cols:qps

	Eigen::Vector2d pos; pos.setZero();
	Eigen::Vector2d neg; neg.setZero();

	double qpseps[5] = { qps[0], qps[1],qps[2],qps[3],qps[4] };

	for (int i = 0; i < 5; i++) {
		qpseps[i] += *eps;
		PhiPsi6A::qps2phipsi(qpseps, &pos);
		qpseps[i] = qps[i];

		qpseps[i] -= *eps;
		PhiPsi6A::qps2phipsi(qpseps, &neg);
		qpseps[i] = qps[i];

		J->col(i) << (pos - neg) / (2 * (*eps));
	}
	//std::cout << "J=\n" << *J << std::endl;
}
void PhiPsi6A::Jacobian(Eigen::MatrixXd *J, Eigen::VectorXd *qpVec, double *eps) {
	double qps[5];
	Eigen::Map<Eigen::VectorXd>(qps, qpVec->rows()) = *qpVec;

	PhiPsi6A::Jacobian(J, qps, eps);
}

PhiPsi11A::PhiPsi11A() {
	KINEMATICPARAMS11A kinParams;
	FwdK11A fklocal(kinParams);
	fk = fklocal;
}
PhiPsi11A::PhiPsi11A(KINEMATICPARAMS11A inParams) {
	kinParams = inParams;
	FwdK11A fklocal(kinParams);
	fk = fklocal;
}
void PhiPsi11A::qps2phipsi(double *qps, Eigen::Vector2d *pp) {
	Eigen::Matrix4d H = fk.qps2H05(qps);
	(*pp)(0) = atan2(H(2, 0), sqrt(pow(H(0, 0), 2) + pow(H(1, 0), 2)));
	(*pp)(1) = atan2(H(1, 0), H(0, 0));
}
void PhiPsi11A::qps2phipsi(double *qps, double *pp) {
	Eigen::Vector2d ppVec;
	ppVec << pp[0], pp[1];
	PhiPsi11A::qps2phipsi(qps, &ppVec);
	Eigen::Map<Eigen::VectorXd>(pp, ppVec.rows()) = ppVec;
}
void PhiPsi11A::qps2phipsi(Eigen::VectorXd *qpVec, Eigen::Vector2d *pp) {
	double qps[5];
	Eigen::Map<Eigen::VectorXd>(qps, qpVec->rows()) = *qpVec;

	PhiPsi11A::qps2phipsi(qps, pp);
}

// Center Sum
CenterSum::CenterSum() {
	JOINTLIMITS jl;
	jntLims = jl;
}
CenterSum::CenterSum(JOINTLIMITS jLims) {
	jntLims = jLims;
}
void CenterSum::qps2ctrSum(double *qps, double *cs) {
	// want large values away from center
	cs[0] = 0; cs[1] = 0;
	double mn;
	double rng;
	for (int i = 0; i < 5; i++) {
		mn = (jntLims.up[i] + jntLims.dn[i]) / 2;
		rng = jntLims.up[i] - jntLims.dn[i];
		cs[0] += pow((qps[i] - mn)/rng,2);

		cs[1] += qps[i]*qps[i]/rng;
	}
}
void CenterSum::qps2ctrSum(double *qps, Eigen::Vector2d *csVec) {
	double cs[5];
	CenterSum::qps2ctrSum(qps, cs);
	(*csVec)(0) = cs[0];
	(*csVec)(1) = cs[1];
}
void CenterSum::qps2ctrSum(Eigen::VectorXd *qpVec, Eigen::Vector2d *csVec) {
	double qps[5];
	Eigen::Map<Eigen::VectorXd>(qps, qpVec->rows()) = *qpVec;
	CenterSum::qps2ctrSum(qps, csVec);
}
void CenterSum::Jacobian(Eigen::MatrixXd *J, double *qps, double *eps) {
	J->setZero(2, 5); //rows:phi,psi cols:qps

	Eigen::Vector2d pos; pos.setZero();
	Eigen::Vector2d neg; neg.setZero();

	double qpseps[5] = { qps[0], qps[1],qps[2],qps[3],qps[4] };

	for (int i = 0; i < 5; i++) {
		qpseps[i] += *eps;
		CenterSum::qps2ctrSum(qpseps, &pos);
		qpseps[i] = qps[i];

		qpseps[i] -= *eps;
		CenterSum::qps2ctrSum(qpseps, &neg);
		qpseps[i] = qps[i];

		J->col(i) << (pos - neg) / (2 * (*eps));
	}
	//std::cout << "J=\n" << *J << std::endl;
}
void CenterSum::Jacobian(Eigen::MatrixXd *J, Eigen::VectorXd *qpVec, double *eps) {
	double qps[5];
	Eigen::Map<Eigen::VectorXd>(qps, qpVec->rows()) = *qpVec;
	CenterSum::Jacobian(J, qps, eps);
}

// Point Center Sum for Augmented
PointCenterSum6A::PointCenterSum6A() {
	KINEMATICPARAMS6A kn;
	kinParams = kn;
	JOINTLIMITS jl;
	jntLims = jl;
}
PointCenterSum6A::PointCenterSum6A(KINEMATICPARAMS6A kn, JOINTLIMITS jLims) {
	kinParams = kn;
	jntLims = jLims;
}
void PointCenterSum6A::qps2ptCtrSum(double *qps, double *pcs) {
	Point6A pt(kinParams);
	CenterSum cs(jntLims);
	
	pt.qps2point(qps, pcs);
	cs.qps2ctrSum(qps, pcs+3);
	//printf("%f %f %f %f %f\n", pcs[0], pcs[1], pcs[2], pcs[3], pcs[4]);
}
void PointCenterSum6A::qps2ptCtrSum(double *qps, Eigen::VectorXd *pcsVec) {
	double pcs[5];
	PointCenterSum6A::qps2ptCtrSum(qps, pcs);
	(*pcsVec)(0) = pcs[0];
	(*pcsVec)(1) = pcs[1];
	(*pcsVec)(2) = pcs[2];
	(*pcsVec)(3) = pcs[3];
	(*pcsVec)(4) = pcs[4];
}
void PointCenterSum6A::qps2ptCtrSum(Eigen::VectorXd *qpVec, Eigen::VectorXd *pcsVec) {
	double qps[5];
	Eigen::Map<Eigen::VectorXd>(qps, qpVec->rows()) = *qpVec;
	PointCenterSum6A::qps2ptCtrSum(qps, pcsVec);
}
void PointCenterSum6A::Jacobian(Eigen::MatrixXd *J, double *qps, double *eps) {
	J->setZero(5, 5); //rows:xyz,ctr,sum cols:qps

	Eigen::VectorXd pos(5); pos.setZero();
	Eigen::VectorXd neg(5); neg.setZero();

	double qpseps[5] = { qps[0], qps[1],qps[2],qps[3],qps[4] };

	for (int i = 0; i < 5; i++) {
		qpseps[i] += *eps;
		PointCenterSum6A::qps2ptCtrSum(qpseps, &pos);
		qpseps[i] = qps[i];

		qpseps[i] -= *eps;
		PointCenterSum6A::qps2ptCtrSum(qpseps, &neg);
		qpseps[i] = qps[i];

		J->col(i) << (pos - neg) / (2 * (*eps));
	}
	//std::cout << "J=\n" << *J << std::endl;
}
void PointCenterSum6A::Jacobian(Eigen::MatrixXd *J, Eigen::VectorXd *qpVec, double *eps) {
	double qps[5];
	Eigen::Map<Eigen::VectorXd>(qps, qpVec->rows()) = *qpVec;
	PointCenterSum6A::Jacobian(J, qps, eps);
}


// Point Phi Psi for Augmented
PointPhiPsi6A::PointPhiPsi6A() {
	KINEMATICPARAMS6A kn;
	kinParams = kn;
	JOINTLIMITS jl;
	jntLims = jl;
}
PointPhiPsi6A::PointPhiPsi6A(KINEMATICPARAMS6A kn, JOINTLIMITS jLims) {
	kinParams = kn;
	jntLims = jLims;
}
void PointPhiPsi6A::qps2ptPhiPsi(double *qps, double *ppp) {
	Point6A pt(kinParams);
	PhiPsi6A pp(kinParams);	

	pt.qps2point(qps, ppp);
	pp.qps2phipsi(qps, ppp + 3);
}
void PointPhiPsi6A::qps2ptPhiPsi(double *qps, Eigen::VectorXd *pppVec) {
	double ppp[5];
	PointPhiPsi6A::qps2ptPhiPsi(qps, ppp);
	(*pppVec)(0) = ppp[0];
	(*pppVec)(1) = ppp[1];
	(*pppVec)(2) = ppp[2];
	(*pppVec)(3) = ppp[3];
	(*pppVec)(4) = ppp[4];
}
void PointPhiPsi6A::qps2ptPhiPsi(Eigen::VectorXd *qpVec, Eigen::VectorXd *pppVec) {
	double qps[5];
	Eigen::Map<Eigen::VectorXd>(qps, qpVec->rows()) = *qpVec;
	PointPhiPsi6A::qps2ptPhiPsi(qps, pppVec);
}
void PointPhiPsi6A::Jacobian(Eigen::MatrixXd *J, double *qps, double *eps) {
	J->setZero(5, 5); //rows:xyz,ctr,sum cols:qps

	Eigen::VectorXd pos(5); pos.setZero();
	Eigen::VectorXd neg(5); neg.setZero();

	double qpseps[5] = { qps[0], qps[1],qps[2],qps[3],qps[4] };

	for (int i = 0; i < 5; i++) {
		qpseps[i] += *eps;
		PointPhiPsi6A::qps2ptPhiPsi(qpseps, &pos);
		qpseps[i] = qps[i];

		qpseps[i] -= *eps;
		PointPhiPsi6A::qps2ptPhiPsi(qpseps, &neg);
		qpseps[i] = qps[i];

		J->col(i) << (pos - neg) / (2 * (*eps));
	}
	//std::cout << "J=\n" << *J << std::endl;
}
void PointPhiPsi6A::Jacobian(Eigen::MatrixXd *J, Eigen::VectorXd *qpVec, double *eps) {
	double qps[5];
	Eigen::Map<Eigen::VectorXd>(qps, qpVec->rows()) = *qpVec;
	PointPhiPsi6A::Jacobian(J, qps, eps);
}


// Kalman 2 State filter
Kalman2State::Kalman2State() {
	ts = .01; // default constructor only here for FData struct, default objects replaced by IK constructors
}
Kalman2State::Kalman2State(double tsSec) {
	ts = tsSec;
}
// the new estimate replaces the old xhat and vhat
void Kalman2State::filterNewMeasure(double *xmea, double *xhat, double *vhat, double *K11, double *K21) {
	double res = *xmea - *xhat - *vhat * ts;
	*xhat = *xhat + ts * *vhat + *K11 * res;
	*vhat = *vhat + *K21*res;
}
