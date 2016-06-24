%% Ben Conrad -- Test the ICRMKinematics DLL in Matlab -- 20160623
setComputer;

%% Load Library
% dpath = '..\x64\Release\ICRMKinematics.dll';
% hpath = '..\ICRMKinematics\kinematicsDLL.h';
dpath = 'D:\interleavedCatheter\ICRMKinematics\x64\Release\ICRMKinematics.dll';
hpath = 'D:\interleavedCatheter\ICRMKinematics\ICRMKinematics\kinematicsDLL.h';

%if : "the specified module could not be found" then it can't find libnlopt.dll; copy into release
%[notfound,warnings] = loadlibrary(dpath, hpath, 'mfilename','i2dx64', 'alias', 'i2dll') to generate a prototype file...bfm
if ~libisloaded('i2dll');
    [notfound,warnings] = loadlibrary(dpath, hpath, 'alias', 'i2dll')
end

%% Test FKs
i06 = inter2D_kn6A;
i11 = inter2D_kn11A;

qps = [1,-.2,.1,3,10]';

kn6a = [i06.pms.cathL, i06.pms.rz01,i06.pms.tx01,i06.pms.ty01,i06.pms.tz01,i06.pms.tx23]';
kn11a = i11.paramStruct2Array( i11.pms );
H06 = reshape(eye(4),16,1); H11 = H06;
%int get6AH01(double *qps, double *kinArray, double *arrayH01)
[ret,qps,kinArray,H06] = calllib( 'i2dll', 'get6AH05', qps, kn6a, H06 );
[ret,qps,kinArray,H11] = calllib( 'i2dll', 'get11AH05', qps, kn11a, H11 );

%unpack from 12 element array
H06(13:15) = H06(10:12);
H06(9:11) = H06(7:9);
H06(5:7) = H06(4:6);
H06(1:3) = H06(1:3);
H06([4,8,12]) = [0,0,0];
H11(13:15) = H11(10:12);
H11(9:11) = H11(7:9);
H11(5:7) = H11(4:6);
H11(1:3) = H11(1:3);
H11([4,8,12]) = [0,0,0];

% [H, h, reshape(i06.forwardK(qps),16,1 )]
H06 = reshape(H06,4,4);
H11 = reshape(H11,4,4);
disp('Compare FKs');
disp([H06, i06.forwardK(qps)])
disp([H11, i11.forwardK(qps)])

%% Test Tasks
%int getTask6A_xyz(double *qps, double *kinArray, double *xyz)
[ret, ~, ~, xyz6a] = calllib( 'i2dll', 'getTask6A_xyz', qps, kn6a, zeros(3,1) );
[ret, ~, ~, xyz11a] = calllib( 'i2dll', 'getTask11A_xyz', qps, kn11a, zeros(3,1) );
disp('task6a_xyz');
disp([xyz6a, xyz11a]);
%int getTask6A_phiPsi(double *qps, double *kinArray, double *pp)
[ret, ~, ~, pp6a] = calllib( 'i2dll', 'getTask6A_phiPsi', qps, kn6a, zeros(2,1) );
[ret, ~, ~, pp11a] = calllib( 'i2dll', 'getTask11A_phiPsi', qps, kn11a, zeros(2,1) );
disp('task6a_pp');
disp([pp6a, pp11a]);
%int getTask6A_xyzuxuyuz(double *qps, double *kinArray, double *xyz, double *uxyz)
[ret, ~, ~, xyz6a, uxyz6a] = calllib( 'i2dll', 'getTask6A_xyzuxuyuz', qps, kn6a, zeros(3,1), zeros(3,1) );
[ret, ~, ~, xyz11a, uxyz11a] = calllib( 'i2dll', 'getTask11A_xyzuxuyuz', qps, kn11a, zeros(3,1), zeros(3,1) );
disp('task6a_xyz');
disp([xyz6a, xyz11a]);

%% Test IKs
nlArray = [1e9,60,0,1e-9,1e-11,1e-9]; %[maxIts,maxTime,method, minFunVal, tolFun, tolX]
jntArray = [-2,-.8,-1,1e-3,0, 2,.8,1,5,500];
qpsG = qps; xyzG = xyz6a; qps0 = zeros(5,1);
% 	// returns -1 = failure
% 	//         -2 = invalid args
% 	//         -3 = out of memory
% 	//         -4 = roundoff limited
% 	//         -5 = forced stop
% 	//          1 = success
% 	//          2 = fun stopval reached
% 	//          3 = fun tol reached
% 	//          4 = x tol reached
% 	//          5 = max evals
% 	//          6 = max time reached
%int getQps_IKnlopt_xyz6A(double *qps, double *kinArray, double *nlArray, double *jntArray, double *xyz)
[ret6, qps6a, ~, ~, ~, xyz6a] = calllib( 'i2dll', 'getQps_IKnlopt_xyz6A', qps0, kn6a, nlArray, jntArray, xyzG);
[ret11, qps11a, ~, ~, ~, xyz11a] = calllib( 'i2dll', 'getQps_IKnlopt_xyz11A', qps0, kn11a, nlArray, jntArray, xyzG);
disp('getQps_IKnlopt_xyz')
disp([0,ret6,ret11; qpsG, qps6a, qps11a; xyzG, xyz6a, xyz11a]);

%% Test IPs
n = 100;
qps = [linspace(0,2,n); linspace(-.8,.8,n); linspace(-1,1,n); linspace(1,3,n); 100*ones(1,n)]';
kn11G = kn11a * 1.05-.1;
HG = zeros(n,4,4);
for i = 1:n;
    HG(i,:,:) = i11.forwardK( qps(i,:), kn11G );
end
stackedQ = reshape(qps', n*5, 1);
stackedU = reshape(HG(:,1:3,1)',n*3,1);
stackedX = reshape(HG(:,1:3,4)',n*3,1);
kn11up = kn11a + .2 + abs(kn11a)*.1;
kn11dn = kn11a - .2 - abs(kn11a)*.1;

%int estimatePmsQ_IPNLOpt_xyzdotu11A(int nSamps, double *stackedQ, double *stackedU, double *stackedX,
% double *k11up, double *k11dn, double *q0Lims, double *nlArray, double *qps0, double *kps0, double *fmin)
[ret, ~,~,~, ~,~, ~, ~, qpsE, knE, fmin] = calllib('i2dll', 'estimatePmsQ_IPNLOpt_xyzdotu11A', ...
    n, stackedQ, stackedU, stackedX, kn11up, kn11dn, jntArray, nlArray, qps0, kn11a, 0 )


%% Unload
unloadlibrary i2dll;