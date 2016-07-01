%% Ben Conrad -- Inter2D Tester -- 20160614
setComputer;

%% Create objects
% i2d{1} = inter2D_kn5A_dll('..\x64\Release\ICRMKinematics.dll','..\ICRMKinematics\kinematicsDLL.h');
i2d{1} = inter2D_kn11A_dll('..\x64\Release\ICRMKinematics.dll','..\ICRMKinematics\kinematicsDLL.h');
% i2d{end+1} = inter2D_kn11A;
% i2d{end+1} = inter2D_kn6A;
n = length(i2d);
    i05a.opts.maxIts = 1e9;
    i05a.opts.maxTime = 60;
    i05a.opts.minFunVal = 1e-12;
    i05a.opts.tolFun = 1e-12;
    i05a.opts.tolX = 1e-12;
    i05a.opts.method = 00; % GN_DIRECT
    %i05a.opts.method = 01; % GN_DIRECT_L --locally biased
    %i05a.opts.method = 02; % GN_DIRECT_L_NOSCAL
    %i05a.opts.method = 03; % GN_DIRECT_L_RAND
    %i05a.opts.method = 04; % GN_DIRECT_L_RAND_NOSCAL
    %i05a.opts.method = 05; % GN_DIRECT_NOSCAL
    %i05a.opts.method = 06; % GN_ESCH
    %i05a.opts.method = 07; % GN_ISRES
    %i05a.opts.method = 08; % GN_MLSL -- slow due to local searches
    %i05a.opts.method = 09; % GN_MLSL_LDS
    %i05a.opts.method = 12; % LN_BOBYQA
    i05a.opts.method = 13; % LN_COBYLA
    %i05a.opts.method = 14; % LN_NelderMead
    %i05a.opts.method = 15; % LN_NEWUOA
    %i05a.opts.method = 16; % LN_NEWUOA_BOUND
    %i05a.opts.method = 17; % LN_PRAXIS
    %i05a.opts.method = 18; % LN_SUBPLX
    i11a.opts = i05a.opts;


%% Check point results
qp0 = [0,0,0,0,10];
for i = 1:n;
    disp(i2d{i}.name)
    disp(i2d{i}.forwardK(qp0));
end

%% Check Inverse Parameter Search
% create 'measured' via altered params
i2d{1}.kns.lCath = 90;
i2d{1}.kns.rz01 = .3;
nms = 100;
qs = [linspace(-.5,.5,nms)',linspace(.2,-.5,nms)',linspace(-.3,.3,nms)', linspace(1,2,nms)', 13*ones(nms,1)];
for i = nms:-1:1;
    Hs(i,:,:) = i2d{1}.forwardK( qs(i,:) + [.1,.2,-.2,.3,5] );
end
i2d{1}.kns.lCath = 95;
%estimate
qpup = [ .5; .3; .3;   1; 10];
qp0  = [  0;  0;  0;.2501;  0];
% qp0  = [  0;  0;  0;.251;  0];%something odd in 11A; fmin = inf unless qp0(4) => .251; not seen in the cpp test or when running from actual measurements
qpdn = [-.5;-.3;-.3;1e-3;-10];
knup = [ 826.000  -46.000   -8.000    0.200   -0.040    0.200    0.200    1.200    1.200  110.000    0.200 ]';
kn0  = [ 806.000  -66.000  -28.000    0.000   -0.240    0.000    0.000    1.000    1.000   95.000    0.000 ]';
kndn = [ 786.000  -86.000  -48.000   -0.200   -0.440   -0.200   -0.200    0.800    0.800   90.000   -0.200 ]';

for i = n:-1:1;
% %     kn0 = i2d{i}.knStruct2Array( i2d{i}.kns );
% %     knup = kn0 + .1*abs(kn0) + .1;
% %     kndn = kn0 - .1*abs(kn0) - .1;
    res{i} = i2d{i}.estimate_qp0kn0_xyzuxuyuz(qs, Hs,  kn0,knup,kndn, qp0,qpup,qpdn );
    disp(res{1})
end

% fprintf('ret = '); for i = 1:n; fprintf('%d ',res{i}.ret); end; fprintf('\n');
% fprintf('fmin = '); for i = 1:n; fprintf('%8.3f ',res{i}.fmin); end; fprintf('\n');
fprintf('Estimating q0:\n');
disp(qp0')
for i = 1:n; disp( res{i}.qp0'); end
fprintf('Estimating kn0:\n');
for i = 1:n; disp( [i2d{i}.knStruct2Array(i2d{i}.kns)'; i2d{i}.knStruct2Array(res{i}.kn0)']); end

%% Check Inverse Kinematic Search
% i2d{1}.kns.cathL = 90;
% goalQps = [.1;.2;.1;2;33];
% qp0 = zeros(5,1);
% for i = 1:n;
%     goal = i2d{i}.taskXYZUxUyUz( i2d{i}.forwardK(goalQps) );
%     res{i} = i2d{i}.findQps_xyzuxuyuz( qp0, goal );
% end
% fprintf('Found qps:\n');
% disp([[0,goalQps']; [res{1}.fval, res{1}.qps']; [res{2}.fval, res{2}.qps']]);


