#include "kn_estimator.h"
#include <iostream>


KNEstimate::KNEstimate(NLOPTPARAMS nParams) {
	nlParams = nParams;
}

double funKNE_AllSingleXYZ(const std::vector<double> &x, std::vector<double> &grad, void *fknData) {
	FKNDATA *ptr;
	ptr = (FKNDATA*)fknData;

	Eigen::Vector3d cmdVec;
	KINEMATICPARAMS6A params;
	params.cathL = x[0];
	params.rz01 = x[1];
	params.tx01 = x[2];
	params.ty01 = x[3];
	params.tz01 = x[4];
	params.tx23 = x[5];

	// FK with the new params
	Point6A pt(params);
	pt.qps2point(ptr->cmdQps, &cmdVec);

	return (cmdVec - ptr->meaVec).norm();
}
int KNEstimate::estimate_AllSingleXYZ(double *cmdQps, double *meaXYZ, double *kinArray) {
	int ret = -1;

	Eigen::Vector3d meaVec;
	meaVec << meaXYZ[0], meaXYZ[1], meaXYZ[2];
	FKNDATA fund;
	fund.cmdQps = cmdQps;
	fund.meaVec = meaVec;

	nlopt::opt lnbob(translateNLAlgForKN(nlParams.method), 6); //there are 6 kinematic params in knParams

	// set objective
	void *ptr;
	ptr = &fund;
	lnbob.set_min_objective(funKNE_AllSingleXYZ, ptr);

	// set initial step size (only used by nongradient methods)
	//std::vector<double> dx0(5);
	//for (int i = 0; i < 5; i++) { dx0[i] = .001; }
	//lnbob.set_initial_step(dx0);

	//set boundary constraints based on qpLast and stdev
	std::vector<double> limup(6);
	std::vector<double> limdn(6);
	limup[0] = 110;  //[mm] cathL
	limdn[0] = 90;
	limup[1] = 0;    //[rad] Rz01
	limdn[1] = -.5;
	limup[2] = 820;  //[mm] Tx01
	limdn[2] = 750;
	limup[3] = -100; //[mm] Ty01
	limdn[3] = -170;
	limup[4] = 0;    //[mm] Tz01
	limdn[4] = -10;
	limup[5] = 10;   //[mm] Tx23
	limdn[5] = 5;
	lnbob.set_upper_bounds(limup);
	lnbob.set_lower_bounds(limdn);

	// set start joint angles to be within bounds
	std::vector<double> x(6);
	for (int i = 0; i < 6; i++) {
		x[i] = kinArray[i];
		if (x[i] <= limdn[i]) {
			x[i] = limdn[i];
		}
		if (limup[i] <= x[i]) {
			x[i] = limup[i];
		}
	}

	// set stopping criteria
	lnbob.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	lnbob.set_ftol_abs(nlParams.tolFunAbs);
	lnbob.set_xtol_abs(nlParams.tolXAbs);
	lnbob.set_maxeval(nlParams.maxIts);
	lnbob.set_maxtime(nlParams.maxTimeSec);

	// solve
	nlopt::result res;
	double fmin = 1e3;
	try {
		res = lnbob.optimize(x, fmin);
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
	for (int i = 0; i < 6; i++) {
		kinArray[i] = x[i];
	}
	ret = res;
	return ret;
}

double funKNE_AllMultipleXYZ(const std::vector<double> &x, std::vector<double> &grad, void *fknData) {
	FKNDATA *ptr;
	ptr = (FKNDATA*)fknData;

	double tempQps[5];
	int off1 = (ptr->nSets) * 1;
	int off2 = (ptr->nSets) * 2;
	int off3 = (ptr->nSets) * 3;
	int off4 = (ptr->nSets) * 4;
	double tempXYZ[3];
	double accum = 0.0;
	double tascXYZ[3];
	double squared = 0;
	
	KINEMATICPARAMS6A params;
	params.cathL = x[0];
	params.rz01 =  x[1];
	params.tx01 =  x[2];
	params.ty01 =  x[3];
	params.tz01 =  x[4];
	params.tx23 =  x[5];

	// FK with the new params
	Point6A pt(params);
	//printf("accum = %7.3f\n", accum);
	for (int i = 0; i < ptr->nSets; i++) {
		tempQps[0] = *(ptr->cmdQpss + i);
		tempQps[1] = *(ptr->cmdQpss + i + off1);
		tempQps[2] = *(ptr->cmdQpss + i + off2);
		tempQps[3] = *(ptr->cmdQpss + i + off3);
		tempQps[4] = *(ptr->cmdQpss + i + off4);
		pt.qps2point(tempQps, tempXYZ);
		tascXYZ[0] = *(ptr->meaXsYsZs + i);
		tascXYZ[1] = *(ptr->meaXsYsZs + i + off1);
		tascXYZ[2] = *(ptr->meaXsYsZs + i + off2);
		//printf("%5d tqps[%+08.6f, %+08.6f, %+08.6f, %+08.6f, %+08.6f]\n", i, tempQps[0], tempQps[1], tempQps[2], tempQps[3], tempQps[4]);
		//printf("%5d txyz[%+08.6f, %+08.6f, %+08.6f]\n", i, tempXYZ[0], tempXYZ[1], tempXYZ[2]);
		//printf("%5d mxyz[%+08.6f, %+08.6f, %+08.6f]\n", i, *(ptr->meaXsYsZs + i), *(ptr->meaXsYsZs + i + off1), *(ptr->meaXsYsZs + i + off2));
		//printf("%5d mxyz[%+08.6f, %+08.6f, %+08.6f]\n", i, tascXYZ[0],tascXYZ[1],tascXYZ[2]);


		//printf("accum = %7.3f\n", accum);
		//accum += pow(tempXYZ[0] - *(ptr->meaXsYsZs + i), 2);
		//accum += abs(tempXYZ[0] - *(ptr->meaXsYsZs + i));
		//printf("accum = %7.3f\n", accum);
		squared = pow(tempXYZ[0] - tascXYZ[0], 2);
		accum += squared;
		//printf("%5d     [%+08.6f] = %+08.6f\n",i, squared, accum);

		//accum += pow(tempXYZ[1] - *(ptr->meaXsYsZs + i + off1), 2);
		//accum += abs(tempXYZ[1] - *(ptr->meaXsYsZs + i + off1));
		//printf("accum = %7.3f\n", accum);
		squared = pow(tempXYZ[1] - tascXYZ[1], 2);
		accum += squared;
		//printf("%5d     [%+08.6f] = %+08.6f\n",i, squared, accum);


		//accum += pow(tempXYZ[2] - *(ptr->meaXsYsZs + i + off2), 2);
		//accum += abs(tempXYZ[2] - *(ptr->meaXsYsZs + i + off2));
		//printf("accum = %7.3f\n", accum);
		squared = pow(tempXYZ[2] - tascXYZ[2], 2);
		accum += squared;
		//printf("%5d     [%+08.6f] = %+08.6f\n",i, squared, accum);

		//printf("%5d x[%+08.6f %+08.6f %+08.6f %+08.6f %+08.6f %+08.6f] q[%+05.4f %+05.4f %+05.4f %+05.4f %+05.4f] = [%+05.4f %+05.4f %+05.4f] = %8.3f\n", i, x[0],x[1],x[2],x[3],x[4],x[5], tempQps[0], tempQps[1], tempQps[2], tempQps[3], tempQps[4], tempXYZ[0], tempXYZ[1], tempXYZ[2], accum);
	}
	//printf("Gives %8.7e\n\n", sqrt(accum));
	return sqrt(accum);
}
int KNEstimate::estimate_AllMultipleXYZ(int nSets, double *cmdQpss, double *meaXsYsZs, double *kinArray, double *kinLim, double *fmin) {
	int ret = -99;

	FKNDATA fund;
	fund.cmdQpss = cmdQpss;
	fund.meaXsYsZs = meaXsYsZs;
	fund.nSets = (const int)nSets;

	nlopt::opt lnbob(translateNLAlgForKN(nlParams.method), 6); //there are 6 kinematic params in knParams

	// set objective
	void *ptr;
	ptr = &fund;
	lnbob.set_min_objective(funKNE_AllMultipleXYZ, ptr);

	// set initial step size (only used by nongradient methods)
	//std::vector<double> dx0(6);
	////for (int i = 0; i < 6; i++) { dx0[i] = .1; }
	//dx0[0] = 1;
	//dx0[1] = .1;
	//dx0[2] = 5;
	//dx0[3] = 5;
	//dx0[4] = 1;
	//dx0[5] = .1;
	//lnbob.set_initial_step(dx0);

	//set boundary constraints based on qpLast and stdev
	std::vector<double> limup(6);
	std::vector<double> limdn(6);
	for (int i = 0; i < 6; i++) {
		limup[i] = kinLim[i];
		limdn[i] = kinLim[i+6];
	}
	//limup[0] = 110;  //[mm] cathL
	//limdn[0] = 90;
	//limup[1] = 0;    //[rad] Rz01
	//limdn[1] = -.5;
	//limup[2] = 820;  //[mm] Tx01
	//limdn[2] = 750;
	//limup[3] = -100; //[mm] Ty01
	//limdn[3] = -170;
	//limup[4] = 0;    //[mm] Tz01
	//limdn[4] = -10;
	//limup[5] = 10;   //[mm] Tx23
	//limdn[5] = 5;
	lnbob.set_upper_bounds(limup);
	lnbob.set_lower_bounds(limdn);
	
	// set start joint angles to be within bounds
	std::vector<double> x(6);
	for (int i = 0; i < 6; i++) {
		x[i] = kinArray[i];
		if (x[i] <= limdn[i]) {
			x[i] = limdn[i];
		}
		if (limup[i] <= x[i]) {
			x[i] = limup[i];
		}
	}

	// set stopping criteria
	lnbob.set_stopval(nlParams.minFunVal); // stop when fmin is less than this
	lnbob.set_ftol_abs(nlParams.tolFunAbs);
	lnbob.set_xtol_abs(nlParams.tolXAbs);
	lnbob.set_maxeval(nlParams.maxIts);
	lnbob.set_maxtime(nlParams.maxTimeSec);
	
	//double temp = 0;
	//printf("fun on points\n");
	//temp = funKNE_AllMultipleXYZ(x, x, ptr);
	//printf("\n\n");
	
	// solve
	nlopt::result res;
	double fmn = 1e9;
	try {
		res = lnbob.optimize(x, fmn);
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
	
	//printf("fun on points\n");
	//temp = funKNE_AllMultipleXYZ(x, x, ptr);
	//printf("\n\n");

	
	// unpack from solver
	for (int i = 0; i < 6; i++) {
		kinArray[i] = x[i];
	}


	//printf("knf [%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f]\n", kinArray[0], kinArray[1], kinArray[2], kinArray[3], kinArray[4], kinArray[5]);
	*fmin = fmn;
	//*fmin = funKNE_AllMultipleXYZ(x, x, ptr);
	ret = res;

	//printf("knx [%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f]\n", x[0], x[1], x[2], x[3], x[4], x[5]);
	//printf("res %d fmin %f\n", res, *fmin);
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//const char fname[] = "D:\\20151008_inter2d_forPaper\\ICRMKinematics\\Release\\funOnPoints.txt";
	//FILE *fid;
	////fid = fopen(fname, "w");
	//fopen_s(&fid, fname, "w");
	//fprintf(fid, "knx [%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f]\n", x[0], x[1], x[2], x[3], x[4], x[5]);
	//fprintf(fid, "res %d fmin %f\n", res, *fmin);
	//

	//double tempQps[5];
	//int off1 = fund.nSets * 1;
	//int off2 = fund.nSets * 2;
	//int off3 = fund.nSets * 3;
	//int off4 = fund.nSets * 4;
	//double tempXYZ[3];
	//double accum = 0.0;
	//double tascXYZ[3];
	//double squared = 0;

	//KINEMATICPARAMS params;
	//params.cathL = x[0];
	//params.rz01 = x[1];
	//params.tx01 = x[2];
	//params.ty01 = x[3];
	//params.tz01 = x[4];
	//params.tx23 = x[5];
	//
	//// FK with the new params
	//Point pt(params);
	////printf("accum = %7.3f\n", accum);
	//for (int i = 0; i <fund.nSets; i++) {
	//	tempQps[0] = *(fund.cmdQpss + i);
	//	tempQps[1] = *(fund.cmdQpss + i + off1);
	//	tempQps[2] = *(fund.cmdQpss + i + off2);
	//	tempQps[3] = *(fund.cmdQpss + i + off3);
	//	tempQps[4] = *(fund.cmdQpss + i + off4);
	//	pt.qps2point(tempQps, tempXYZ); 
	//	tascXYZ[0] = *(fund.meaXsYsZs + i);
	//	tascXYZ[1] = *(fund.meaXsYsZs + i + off1);
	//	tascXYZ[2] = *(fund.meaXsYsZs + i + off2);
	//	fprintf(fid, "%5d tqps[%+08.6f, %+08.6f, %+08.6f, %+08.6f, %+08.6f]\n", i, tempQps[0], tempQps[1], tempQps[2], tempQps[3], tempQps[4]);
	//	fprintf(fid, "%5d txyz[%+08.6f, %+08.6f, %+08.6f]\n", i, tempXYZ[0], tempXYZ[1], tempXYZ[2]);
	//	fprintf(fid, "%5d mxyz[%+08.6f, %+08.6f, %+08.6f]\n", i, tascXYZ[0], tascXYZ[1], tascXYZ[2]);


	//	//printf("accum = %7.3f\n", accum);
	//	//accum += pow(tempXYZ[0] - *(ptr->meaXsYsZs + i), 2);
	//	//accum += abs(tempXYZ[0] - *(ptr->meaXsYsZs + i));
	//	//printf("accum = %7.3f\n", accum);
	//	squared = pow(tempXYZ[0] - tascXYZ[0], 2);
	//	accum += squared;
	//	fprintf(fid, "%5d     [%+08.6f] = %+08.6f\n", i, squared, accum);

	//	//accum += pow(tempXYZ[1] - *(ptr->meaXsYsZs + i + off1), 2);
	//	//accum += abs(tempXYZ[1] - *(ptr->meaXsYsZs + i + off1));
	//	//printf("accum = %7.3f\n", accum);
	//	squared = pow(tempXYZ[1] - tascXYZ[1], 2);
	//	accum += squared;
	//	fprintf(fid, "%5d     [%+08.6f] = %+08.6f\n", i, squared, accum);


	//	//accum += pow(tempXYZ[2] - *(ptr->meaXsYsZs + i + off2), 2);
	//	//accum += abs(tempXYZ[2] - *(ptr->meaXsYsZs + i + off2));
	//	//printf("accum = %7.3f\n", accum);
	//	squared = pow(tempXYZ[2] - tascXYZ[2], 2);
	//	accum += squared;
	//	fprintf(fid, "%5d     [%+08.6f] = %+08.6f\n", i, squared, accum);

	//	fprintf(fid, "%5d x[%+08.6f %+08.6f %+08.6f %+08.6f %+08.6f %+08.6f] q[%+05.4f %+05.4f %+05.4f %+05.4f %+05.4f] = [%+05.4f %+05.4f %+05.4f] = %8.3f\n", i, x[0], x[1], x[2], x[3], x[4], x[5], tempQps[0], tempQps[1], tempQps[2], tempQps[3], tempQps[4], tempXYZ[0], tempXYZ[1], tempXYZ[2], accum);
	//}
	//printf("Gives %8.7e\n\n", sqrt(accum));
	//fclose(fid);
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	return ret;
}



nlopt::algorithm translateNLAlgForKN(nlMethod method) {
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