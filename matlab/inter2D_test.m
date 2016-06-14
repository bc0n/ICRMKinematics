%% Ben Conrad -- Inter2D Tester -- 20160614
setComputer;

%% Create objects
i2d{2} = inter2D_kn11A;
i2d{1} = inter2D_kn6A;
n = length(i2d);

%% Check point results
qp0 = [0,0,0,0,10];

for i = 1:n;
    disp(i2d{i}.name)
    disp(i2d{i}.forwardK(qp0));
end

%% Check Inverse Parameter Search
q0 = [0;.1;0;.1;0];
qs = [linspace(-.5,.5,10)',linspace(.2,-.5,10)',linspace(-.3,.3,10)', linspace(1,2,10)', 13*ones(10,1)];
for i = 10:-1:1;
    Hs(i,:,:) = i2d{1}.forwardK( qs(i,:)'+q0 );
end

qp0 = zeros(5,1);
qpPM = [.2;.2;.2;.5;10];
for i = n:-1:1;
    pm0 = i2d{i}.paramStruct2Array( i2d{i}.pms );
    pm0 = pm0*1.001;
    pmPM = abs(pm0*.05);
    res{i} = i2d{i}.findQp0Pm0_xyzuxuyuz(Hs, qs, qp0,qpPM, pm0,pmPM);
    %res{i} = i2d{i}.findQp0Pm0_xyzdotu(Hs, qs, qp0,qpPM, pm0,pmPM);
end

fprintf('Estimating q0:\n');
disp([q0, res{1}.qp0, res{2}.qp0]');
fprintf('Estimating pms0:\n');
for i = 1:n;
    disp( [i2d{i}.paramStruct2Array( i2d{i}.pms)'; res{i}.pm0']);
end

%% Check Inverse Kinematic Search
i2d{1}.pms.cathL = 90;
goalQps = [.1;.2;.1;2;33];
qp0 = zeros(5,1);
for i = 1:n;
    goal = i2d{i}.taskXYZUxUyUz( i2d{i}.forwardK(goalQps) );
    res{i} = i2d{i}.findQps_xyzuxuyuz( qp0, goal );
end
fprintf('Found qps:\n');
disp([[0,goalQps']; [res{1}.fval, res{1}.qps']; [res{2}.fval, res{2}.qps']]);


