#pragma once

#include "kinematics_structs.h"
#include "forwardKinematics.h"


//http://www.cplusplus.com/forum/general/186/
//Both the declaration and definition of templates must be included. This is because template functions cannot be compiled and linked independently, since they are generated on request for the specific types they are instantiated with.
template <class TFK> class TaskXYZ { // TFK = the type of FK model
private:
	TFK tfk;
	Eigen::Matrix4d H;
public:
	TaskXYZ(TFK fkArg) { tfk = fkArg; }
	void qps2task(double *qps, double *xyz) {
		H = tfk.qps2H05(qps);
		xyz[0] = H(0, 3);
		xyz[1] = H(1, 3);
		xyz[2] = H(2, 3);
	}
	double qps2taskError(double *qps, double *target) {
		//find task
		H = tfk.qps2H05(qps);

		//find the task error and return its norm
		double res = 0;
		for (int i = 0; i < 3; i++) {
			res += pow(target[i] - H(i,3), 2);
		}
		return sqrt(res);
	}
};

template <class TFK> class TaskPhiPsi {// TFK = the type of FK model
private:
	TFK tfk;
	Eigen::Matrix4d H;
public:
	TaskPhiPsi(TFK fkArg) { tfk = fkArg; }
	void qps2task(double *qps, double *pp) {
		H = tfk.qps2H05(qps);
		*(pp+0) = atan2(H(2, 0), sqrt(pow(H(0, 0), 2) + pow(H(1, 0), 2)));
		*(pp+1) = atan2(H(1, 0), H(0, 0));
	}
};

//center of joint ranges, minimize summed joint angles
template <class TFK> class TaskCenterSum {// TFK = the type of FK model
private:
	TFK tfk;
	JOINTLIMITS jntLims;
public:
	TaskCenterSum(JOINTLIMITS jLims, TFK fkArg) {
		jntLims = jLims;
		tfk = fkArg;
	}
	void qps2task(double *qps, double *cs) {
		// want large values away from center
		cs[0] = 0; cs[1] = 0;
		double mn;
		double rng;
		for (int i = 0; i < 5; i++) {
			mn = (jntLims.up[i] + jntLims.dn[i]) / 2;
			rng = jntLims.up[i] - jntLims.dn[i];
			cs[0] += pow((qps[i] - mn) / rng, 2);

			cs[1] += qps[i] * qps[i] / rng;
		}
	}
};

template <class TFK> class TaskXYZUxUyUz {// TFK = the type of FK model
private:
	TFK tfk;
	Eigen::Matrix4d H;
public:
	TaskXYZUxUyUz(TFK fkArg) { tfk = fkArg; }
	void qps2task(double *qps, double *xyzuxuyuz) {
		H = tfk.qps2H05(qps);
		xyzuxuyuz[0] = H(0, 3);
		xyzuxuyuz[1] = H(1, 3);
		xyzuxuyuz[2] = H(2, 3);
		xyzuxuyuz[3] = H(0, 0);
		xyzuxuyuz[4] = H(1, 0);
		xyzuxuyuz[5] = H(2, 0);
	}

	double qps2taskError(double *qps, double *target) {
		//find task
		H = tfk.qps2H05(qps);
		//find the task error and return its norm
		double res = 0;
		for (int i = 0; i < 3; i++) {
			res += pow(target[i] - H(i, 3), 2); //x,y,z
			res += pow(target[i+3] - H(i, 0), 2); //ux,uy,uz
		}
		return sqrt(res);
	}
};




//from http://eigen.tuxfamily.org/bz/show_bug.cgi?id=257#c14
template<typename _Matrix_Type_>_Matrix_Type_ pseudoInverse(const _Matrix_Type_ &a, double epsilon = std::numeric_limits<double>::epsilon()) {
	Eigen::JacobiSVD< _Matrix_Type_ > svd(a, Eigen::ComputeThinU | Eigen::ComputeThinV);
	double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
	return svd.matrixV() *  (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
}