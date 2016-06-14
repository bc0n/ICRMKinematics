#pragma once

#include "kinematics_structs.h"
#include "forwardKinematics.h"

class Point6A{
public:
	KINEMATICPARAMS6A kinParams;
	FwdK6A fk;

	Point6A(void);
	Point6A(KINEMATICPARAMS6A inParams);

	void qps2point(double *qps, Eigen::Vector3d *xyzVec);
	void qps2point(double *qps, double *xyz);
	void qps2point(Eigen::VectorXd *qpVec, Eigen::Vector3d *xyzVec);

	void Jacobian(Eigen::MatrixXd *J, double *qps, double *eps);
	void Jacobian(Eigen::MatrixXd *J, Eigen::VectorXd *qpVec, double *eps);
};
class Point11A {
public:
	KINEMATICPARAMS11A kinParams;
	FwdK11A fk;

	Point11A(void);
	Point11A(KINEMATICPARAMS11A inParams);

	void qps2point(double *qps, Eigen::Vector3d *xyzVec);
	void qps2point(double *qps, double *xyz);
	void qps2point(Eigen::VectorXd *qpVec, Eigen::Vector3d *xyzVec);
};

class PhiPsi6A {
public:
	KINEMATICPARAMS6A kinParams;
	FwdK6A fk;

	PhiPsi6A();
	PhiPsi6A(KINEMATICPARAMS6A inParams);

	void qps2phipsi(double *qps, Eigen::Vector2d *ppVec);
	void qps2phipsi(double *qps, double *pp);
	void qps2phipsi(Eigen::VectorXd *qpVec, Eigen::Vector2d *ppVec);

	void Jacobian(Eigen::MatrixXd *J, double *qps, double *eps);
	void Jacobian(Eigen::MatrixXd *J, Eigen::VectorXd *qpVec, double *eps);
};
class PhiPsi11A {
public:
	KINEMATICPARAMS11A kinParams;
	FwdK11A fk;

	PhiPsi11A();
	PhiPsi11A(KINEMATICPARAMS11A inParams);

	void qps2phipsi(double *qps, Eigen::Vector2d *ppVec);
	void qps2phipsi(double *qps, double *pp);
	void qps2phipsi(Eigen::VectorXd *qpVec, Eigen::Vector2d *ppVec);
};

//center of joint ranges, minimize summed joint angles
class CenterSum {
private:
	JOINTLIMITS jntLims;

public:
	CenterSum();
	CenterSum(JOINTLIMITS jLims);

	void qps2ctrSum(double *qps, Eigen::Vector2d *csVec);
	void qps2ctrSum(double *qps, double *cs);
	void qps2ctrSum(Eigen::VectorXd *qpVec, Eigen::Vector2d *csVec);

	void Jacobian(Eigen::MatrixXd *J, double *qps, double *eps);
	void Jacobian(Eigen::MatrixXd *J, Eigen::VectorXd *qpVec, double *eps);
};

class PointCenterSum6A {
private:
	KINEMATICPARAMS6A kinParams;
	JOINTLIMITS jntLims;

public:
	PointCenterSum6A();
	PointCenterSum6A(KINEMATICPARAMS6A kn, JOINTLIMITS jLims);

	void qps2ptCtrSum(double *qps, Eigen::VectorXd *pcsVec);
	void qps2ptCtrSum(double *qps, double *pcs);
	void qps2ptCtrSum(Eigen::VectorXd *qpVec, Eigen::VectorXd *pcsVec);

	void Jacobian(Eigen::MatrixXd *J, double *qps, double *eps);
	void Jacobian(Eigen::MatrixXd *J, Eigen::VectorXd *qpVec, double *eps);
};

class PointPhiPsi6A {
private:
	KINEMATICPARAMS6A kinParams;
	JOINTLIMITS jntLims;

public:
	PointPhiPsi6A();
	PointPhiPsi6A(KINEMATICPARAMS6A kn, JOINTLIMITS jLims);

	void qps2ptPhiPsi(double *qps, Eigen::VectorXd *pppVec);
	void qps2ptPhiPsi(double *qps, double *ppp);
	void qps2ptPhiPsi(Eigen::VectorXd *qpVec, Eigen::VectorXd *pppVec);

	void Jacobian(Eigen::MatrixXd *J, double *qps, double *eps);
	void Jacobian(Eigen::MatrixXd *J, Eigen::VectorXd *qpVec, double *eps);
};

class Kalman2State {
private:
	double ts;
public:
	Kalman2State();
	Kalman2State(double tsSec);
	void filterNewMeasure(double *xmea, double *xhat, double *vhat, double *K11, double *K21); // new estimate placed in xhat, vhat
};





// http://eigen.tuxfamily.org/bz/show_bug.cgi?id=257#c14
template<typename _Matrix_Type_>_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon()) {
	Eigen::JacobiSVD< _Matrix_Type_ > svd(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
	double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
	return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}