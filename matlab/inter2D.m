% This class creates an inter2D model and provides for generating trajectories, solving the kinematics, and plotting the results
classdef inter2D
  
    properties(GetAccess = 'public', SetAccess = 'public')
        pms;    %kinematic params
        eps;    %Numeric Jacobian delta
        maxIts; %maximum Newton Raphson iterations
        errTol; %error tolerance
        factor; %NR update factor
        
        qp0;    %initial joint positions
        jlims;  %joint limits
        
        opts;
        
        name;
        nums; %numbers of things
        
    end
    
    methods
        function obj = inter2D
            obj.eps = 1e-3;
            obj.maxIts = 1e2;
            obj.errTol = 1e-2;
            obj.factor = 1;
            obj.qp0 = .1*ones(5,1);
            obj.jlims.up = [pi,1,pi,5,300];
            obj.jlims.dn = [-pi,-1,-pi,1e-3,1];

            obj.opts = optimset('fmincon');
            obj.opts.Display = 'off';
        end
        
%---task funs---
        function task = taskXYZ(obj, H) %#ok<INUSL>
            task = H(1:3,4);
        end
        function task = taskXYZUxUyUz(obj, H) %#ok<INUSL>
            task = [H(1:3,4);H(1:3,1)];
        end
        function task = taskXYZPP(obj, H)
            task = obj.H2xyzpp( H );
        end
        
        %Compute the numeric Jacobian for taskFun about operating point qp
        function J = numQJ(obj,qp,taskFun)
            J = [ taskFun([qp(1)+obj.dt, qp(2), qp(3), qp(4), qp(5)]) - taskFun([qp(1)-obj.dt, qp(2), qp(3), qp(4), qp(5)]), ...
                  taskFun([qp(1), qp(2)+obj.dt, qp(3), qp(4), qp(5)]) - taskFun([qp(1), qp(2)-obj.dt, qp(3), qp(4), qp(5)]), ...
                  taskFun([qp(1), qp(2), qp(3)+obj.dt, qp(4), qp(5)]) - taskFun([qp(1), qp(2), qp(3)-obj.dt, qp(4), qp(5)]), ...
                  taskFun([qp(1), qp(2), qp(3), qp(4)+obj.dt, qp(5)]) - taskFun([qp(1), qp(2), qp(3), qp(4)-obj.dt, qp(5)]), ...
                  taskFun([qp(1), qp(2), qp(3), qp(4), qp(5)+obj.dt]) - taskFun([qp(1), qp(2), qp(3), qp(4), qp(5)-obj.dt]) ] /2/obj.dt;
        end %numJ
        
        
%---kinematic parameter error functions---
        function errSum = errQp0Pm0_xyzuxuyuz( obj, qp0, pm0, meaHs, meaQs )
            qp0 = reshape(qp0, 1,5);
            n = size(meaHs,1);
            errSum = 0;
            for it = 1:n;
                H = obj.forwardK( qp0 + meaQs(it,:), pm0 );
                errSum = errSum + norm(H(1:3,4)-meaHs(it,1:3,4)') + norm(H(1:3,1)-meaHs(it,1:3,1)');
            end
            
        end
        
        function errSum = errQp0Pm0_xyzdotu( obj, qp0, pm0, meaHs, meaQs )
            qp0 = reshape(qp0, 1,5);
            n = size(meaHs,1);
            His = zeros(n,4,4);
            for it = 1:n;
                His(it,:,:) = obj.forwardK( qp0 + meaQs(it,:), pm0 );
            end
            errSum = sum(abs(1-dot( meaHs(:,1:3,1), His(:,1:3,1), 2)) ) / n;
            errSum = errSum + sum( sqrt(sum((meaHs(:,1:3,4)-His(:,1:3,4)).^2,2)) ) / n;
        end
        
%---kinematic parameter search---
        function res = findQp0Pm0_xyzuxuyuz(obj, Hs, qs, qp0,qpPM, pm0,pmPM)
            x0 = [qp0;pm0];
            xp = [qp0 + qpPM; pm0 + pmPM];
            xn = [qp0 - qpPM; pm0 - pmPM];
            
            fun = @(x)obj.errQp0Pm0_xyzuxuyuz( x(1:obj.nums.qps), x(obj.nums.qps+(1:obj.nums.pms)), Hs, qs);
            
            tic
            [out.x,out.f,out.flag, out.message] = fmincon( fun, x0, [],[],[],[], xn,xp, [], obj.opts );
            res.telapsed = toc;
            
            res.qp0 = out.x(1:obj.nums.qps);
            res.pm0 = out.x(obj.nums.qps+(1:obj.nums.pms));
            res.fval = out.f;
            res.flag = out.flag;
            res.msg = out.message;
        end

        function res = findQp0Pm0_xyzdotu(obj, Hs, qs, qp0,qpPM, pm0,pmPM)
            x0 = [qp0;pm0];
            xp = [qp0 + qpPM; pm0 + pmPM];
            xn = [qp0 - qpPM; pm0 - pmPM];
            
            fun = @(x)obj.errQp0Pm0_xyzdotu( x(1:obj.nums.qps), x(obj.nums.qps+(1:obj.nums.pms)), Hs, qs);
            
            tic
            [out.x,out.f,out.flag, out.message] = fmincon( fun, x0, [],[],[],[], xn,xp, [], obj.opts );
            res.telapsed = toc;
            
            res.qp0 = out.x(1:obj.nums.qps);
            res.pm0 = out.x(obj.nums.qps+(1:obj.nums.pms));
            res.fval = out.f;
            res.flag = out.flag;
            res.msg = out.message;
        end

%---inverse kinematic error functions
        function err = errQps_xyzuxuyuz( obj, qps, goalXYZUxUyUz )
            err = norm(goalXYZUxUyUz - obj.taskXYZUxUyUz( obj.forwardK(qps) ));            
        end
        
%---inverse kinematic search---
        function res = findQps_xyzuxuyuz(obj, qp0, goalXYZUxUyUz)
            x0 = qp0;
            xn = obj.jlims.dn;
            xp = obj.jlims.up;
            
            fun = @(x)obj.errQps_xyzuxuyuz( x, goalXYZUxUyUz);
            
            tic
            [out.x,out.f,out.flag, out.message] = fmincon( fun, x0, [],[],[],[], xn,xp, [], obj.opts );
            res.telapsed = toc;
            
            res.qps = out.x(1:obj.nums.qps);
            res.fval = out.f;
            res.flag = out.flag;
            res.msg = out.message;
            res.goal = goalXYZUxUyUz;
        end

        
%---plotting stuff---
        % cathPts = obj.xyzCathPts( qp )
        % cathPts = [ [base],...[x;y;z;1],...[tip] ]
        function cathPts = xyzCathPts(obj,qp)
            [H01,H02,H03,H04,H05] = obj.forwardK( qp );
            r = obj.cathL/qp(4); %radius of catheter arc
            cathPts = [ H03* [0;0;0;1], ...
                  H03* obj.Ty(r*(1-cos(qp(4)*.1)))*obj.Tx(r*sin(qp(4)*.1))*obj.Rz(qp(4)*.1)*[0;0;0;1], ...
                  H03* obj.Ty(r*(1-cos(qp(4)*.2)))*obj.Tx(r*sin(qp(4)*.2))*obj.Rz(qp(4)*.2)*[0;0;0;1], ...
                  H03* obj.Ty(r*(1-cos(qp(4)*.3)))*obj.Tx(r*sin(qp(4)*.3))*obj.Rz(qp(4)*.3)*[0;0;0;1], ...
                  H03* obj.Ty(r*(1-cos(qp(4)*.4)))*obj.Tx(r*sin(qp(4)*.4))*obj.Rz(qp(4)*.4)*[0;0;0;1], ...
                  H03* obj.Ty(r*(1-cos(qp(4)*.5)))*obj.Tx(r*sin(qp(4)*.5))*obj.Rz(qp(4)*.5)*[0;0;0;1], ...
                  H03* obj.Ty(r*(1-cos(qp(4)*.6)))*obj.Tx(r*sin(qp(4)*.6))*obj.Rz(qp(4)*.6)*[0;0;0;1], ...
                  H03* obj.Ty(r*(1-cos(qp(4)*.7)))*obj.Tx(r*sin(qp(4)*.7))*obj.Rz(qp(4)*.7)*[0;0;0;1], ...
                  H03* obj.Ty(r*(1-cos(qp(4)*.8)))*obj.Tx(r*sin(qp(4)*.8))*obj.Rz(qp(4)*.8)*[0;0;0;1], ...
                  H03* obj.Ty(r*(1-cos(qp(4)*.9)))*obj.Tx(r*sin(qp(4)*.9))*obj.Rz(qp(4)*.9)*[0;0;0;1], ...
                  H03* obj.Ty(r*(1-cos(qp(4)*1.)))*obj.Tx(r*sin(qp(4)*1.))*obj.Rz(qp(4)*1.)*[0;0;0;1] ];
        end
        
        function drawConfig(obj, qps )
            [h01,h02,h03,h04,h05] = obj.forwardK(qps);
            
            %proximal: rectangles in local xy and xz, starting at h01 (pitch axis)
            hxyTR = h01*obj.Ty(6); hxyTL = h01*obj.Ty(6)*obj.Tx(-30); hxyBL = h01*obj.Ty(-6)*obj.Tx(-30); hxyBR = h01*obj.Ty(-6);
            xy = [h01(1:3,4), hxyTR(1:3,4), hxyTL(1:3,4), hxyBL(1:3,4), hxyBR(1:3,4), h01(1:3,4) ];
            plot3(xy(1,:), xy(2,:), xy(3,:), '-b');
            hxzTR = h01*obj.Tz(6); hxzTL = h01*obj.Tz(6)*obj.Tx(-30); hxzBL = h01*obj.Tz(-6)*obj.Tx(-30); hxzBR = h01*obj.Tz(-6);
            xz = [h01(1:3,4), hxzTR(1:3,4), hxzTL(1:3,4), hxzBL(1:3,4), hxzBR(1:3,4), h01(1:3,4) ];
            plot3(xz(1,:), xz(2,:), xz(3,:), '-b');
            %pitch: rectangles in local xy and xz between pitch axis and roll axis
            lPitch = 4;
            hxyTL = h02*obj.Ty(6); hxyTR = h02*obj.Ty(6)*obj.Tx(lPitch); hxyBR = h02*obj.Ty(-6)*obj.Tx(lPitch); hxyBL = h02*obj.Ty(-6);
            xy = [h02(1:3,4), hxyTL(1:3,4), hxyTR(1:3,4), hxyBR(1:3,4), hxyBL(1:3,4), h02(1:3,4)];
            plot3(xy(1,:), xy(2,:), xy(3,:), '-c');
            hxzTL = h02*obj.Tz(6); hxzTR = h02*obj.Tz(6)*obj.Tx(lPitch); hxzBR = h02*obj.Tz(-6)*obj.Tx(lPitch); hxzBL = h02*obj.Tz(-6);
            xz = [h02(1:3,4), hxzTL(1:3,4), hxzTR(1:3,4), hxzBR(1:3,4), hxzBL(1:3,4), h02(1:3,4)];
            plot3(xz(1,:), xz(2,:), xz(3,:), '-c');
            %roll: rectangles in local xy and xz from the roll axis to the catheter start
            lRoll = 8.2-lPitch;
            hxyTL = h03*obj.Ty(6); hxyTR = h03*obj.Ty(6)*obj.Tx(lRoll); hxyBR = h03*obj.Ty(-6)*obj.Tx(lRoll); hxyBL = h03*obj.Ty(-6);
            xy = [h03(1:3,4), hxyTL(1:3,4), hxyTR(1:3,4), hxyBR(1:3,4), hxyBL(1:3,4), h03(1:3,4)];
            plot3(xy(1,:), xy(2,:), xy(3,:), '-g');
            hxzTL = h03*obj.Tz(6); hxzTR = h03*obj.Tz(6)*obj.Tx(lRoll); hxzBR = h03*obj.Tz(-6)*obj.Tx(lRoll); hxzBL = h03*obj.Tz(-6);
            xz = [h03(1:3,4), hxzTL(1:3,4), hxzTR(1:3,4), hxzBR(1:3,4), hxzBL(1:3,4), h03(1:3,4)];
            plot3(xz(1,:), xz(2,:), xz(3,:), '-g');
        end

    end %dynamic methods
    
    methods (Static)
        % Rth = @(th)[cos(th),0,-sin(th); 0,1,0; sin(th),0,cos(th)]; %about x
        % Rps = @(ps)[cos(ps),sin(ps),0; -sin(ps),cos(ps),0; 0,1,0]; %about y
        % Rph = @(ph)[cos(ph),sin(ph),0; -sin(ph),cos(ph),0; 0,0,1]; %about z
        function H = Rx(a)
            H = [1,0,0,0; 0,cos(a),-sin(a),0; 0,sin(a),cos(a),0; 0,0,0,1]; %rotation about x
        end
        function H = Ry(b)
            H = [cos(b),0,sin(b),0; 0,1,0,0; -sin(b),0,cos(b),0; 0,0,0,1]; %about y
        end
        function H = Rz(c)
            H = [cos(c),-sin(c),0,0; sin(c),cos(c),0,0; 0,0,1,0; 0,0,0,1]; %about z
        end
        function H = Tx(a)
            H = [1,0,0,a; 0,1,0,0; 0,0,1,0; 0,0,0,1]; %translation along x
        end
        function H = Ty(b)
            H = [1,0,0,0; 0,1,0,b; 0,0,1,0; 0,0,0,1]; %along y
        end
        function H = Tz(c)
            H = [1,0,0,0; 0,1,0,0; 0,0,1,c; 0,0,0,1]; %along z
        end
        function H = Rtk(th, k)
            if round(norm(k)*1e3)/1e3 ~= 1; error('MATLAB:Rtk','k must have unit norm'); end;
            H = [ k(1)^2*(1-cos(th))+cos(th),         k(1)*k(2)*(1-cos(th))-k(3)*sin(th), k(1)*k(3)*(1-cos(th))+k(2)*sin(th), 0;... %the axis(k)&angle(th) tranformation matrix
                  k(1)*k(2)*(1-cos(th))+k(3)*sin(th), k(2)^2*(1-cos(th))+cos(th),         k(2)*k(3)*(1-cos(th))-k(1)*sin(th), 0;...
                  k(1)*k(3)*(1-cos(th))-k(2)*sin(th), k(2)*k(3)*(1-cos(th))+k(1)*sin(th), k(3)^2*(1-cos(th))+cos(th),         0; 0,0,0,1];
        end
        function H = Raer( a, e, r) %ZYZ Euler Rotation
            %azimuth = phi = about Z0
            %elevation = theta = about  Y1
            %roll = psi = about Z2
            H = inter2D.Rz(a)*inter2D.Ry(e)*inter2D.Rz(r);
        end
        function H = Rasc( a, e, r) %The unknown Ascension rotation matrix
            %this is similar to an RPY rotation Rx(r)Ry(e)Rz(a) except for some signs...NOT ZYZ Euler
            H = eye(4);
            ca = cos(a); sa = sin(a);
            ce = cos(e); se = sin(e);
            cr = cos(r); sr = sin(r);
            H(1,1) =  ce*ca;
            H(2,1) = -cr*sa+sr*se*ca;
            H(3,1) =  sr*sa+cr*se*ca;
            H(1,2) =  ce*sa;
            H(2,2) =  cr*ca+sr*se*sa;
            H(2,3) = -sr*ca+cr*se*sa;
            H(3,1) = -se;
            H(3,2) =  sr*ce;
            H(3,3) =  cr*ce;
        end
        function H = Rquat( q )
            %from Robot Modelling Planning and Control...differs in signs from Kammer's
%             n = q(1); ex = q(2); ey = q(3); ez = q(4);
%             H = eye(4);
%             H(1,1) = 2*(n^2+ex^2)-1; %q(1)^2+q(2)^2-q(3)^2-q(4)^2
%             H(2,1) = 2*(ex*ey+n*ez); %2*(q(2)*q(3)-q(1)*q(4))
%             H(3,1) = 2*(ex*ez-n*ey); %2*(q(2)*q(4)+q(1)*q(3))
%             H(1,2) = 2*(ex*ey-n*ez); %2*(q(2)*q(3)+q(1)*q(4))
%             H(2,2) = 2*(n^2+ey^2)-1; %q(1)^2-q(2)^2+q(3)^2-q(4)^2
%             H(3,2) = 2*(ey*ez+n*ex); %2*(q(3)*q(4)-q(1)*q(2))
%             H(1,3) = 2*(ex*ez+n*ey); %2*(q(2)*q(4)-q(1)*q(3))
%             H(2,3) = 2*(ey*ez-n*ex); %2*(q(3)*q(4)+q(1)*q(2))
%             H(3,3) = 2*(n^2+ez^2)-1; %q(1)^2-q(2)^2-q(3)^2+q(4)^2
            
            %from Kammer, 2008_Markley_UnitQuaternionFromRotationMatrix
            H=  [ q(1)^2+q(2)^2-q(3)^2-q(4)^2, 2*(q(2)*q(3)+q(1)*q(4)),     2*(q(2)*q(4)-q(1)*q(3)),     0;...
                  2*(q(2)*q(3)-q(1)*q(4)),     q(1)^2-q(2)^2+q(3)^2-q(4)^2, 2*(q(3)*q(4)+q(1)*q(2)),     0;...
                  2*(q(2)*q(4)+q(1)*q(3)),     2*(q(3)*q(4)-q(1)*q(2)),     q(1)^2-q(2)^2-q(3)^2+q(4)^2, 0; 0,0,0,1];
        end
        
        function q = R2Quat( R )
            q(1) = .5 * sqrt( R(1,1)+R(2,2)+R(3,3) +1 );
            q(2) = .5 * sign(R(3,2)-R(2,3)) * sqrt(R(1,1)-R(2,2)-R(3,3)+1);
            q(3) = .5 * sign(R(1,3)-R(3,1)) * sqrt(R(2,2)-R(3,3)-R(1,1)+1);
            q(4) = .5 * sign(R(2,1)-R(1,2)) * sqrt(R(3,3)-R(1,1)-R(2,2)+1);
        end
        
        function pp = phiPsi(varg)
            ux = varg(1); uy = varg(2); uz = varg(3);
            pp(1) = atan2( uz, norm([ux,uy]) ); %phi = angle above the xy plane
            pp(2) = atan2( uy, ux ); %psi = angle in the xy plane
        end
        function xyzpp = H2xyzpp(H)
            xyzpp = zeros(5,1);
            xyzpp(1:3) = H(1:3,4);
            xyzpp(4:5) = inter2D.phiPsi(H(1:3,1));
        end
        
        
        function txt = addDisplayNameDatatip(~,event_obj,hdl)
            pos = get(event_obj,'Position');
            for ind = 1:length(hdl);
                %if hdl(ind).XData(2) == pos(1) && hdl(ind).YData(2) == pos(2)
                if hdl(ind).XData(1) == pos(1) && hdl(ind).YData(1) == pos(2)
                    break;
                end;
            end;
            txt = sprintf('%s @[%5.4f,%5.4f,%5.4f]',hdl(ind).DisplayName,pos(1),pos(2),pos(3));
        end %datatip
        
        %plot single pose
        function h = quiver3Pose( xyz, vec, style )
            vec = vec*style.scl; %scale here and 0 below to avoid 'prevent overlapping'
            h = quiver3( xyz(1),xyz(2),xyz(3), vec(1),vec(2),vec(3), 0,...
                    'LineWidth',style.lwd, 'Color',style.lcol, 'LineStyle',style.lsty );
        end
        function h = quiver3Hx(varargin)
            if nargin == 2;
                H = varargin{1};
                in = varargin{2};
                h = inter2D.quiver3Pose( H(1:3,4),H(1:3,1), in );
            else
                h.scl = 1;
                h.lwd = 1;
                h.lcol = [1,0,0];
                h.lsty = '-';
                h.msz = 5;
                h.mcol = 'none';
                h.msty = 'o';
            end
        end %quiver3Hx
        function h = quiver3Hy(varargin)
            if nargin == 2;
                H = varargin{1};
                in = varargin{2};
                h = inter2D.quiver3Pose( H(1:3,4),H(1:3,2), in );
            else
                h.scl = 1;
                h.lwd = 1;
                h.lcol = [0,1,0];
                h.lsty = '-';
                h.msz = 5;
                h.mcol = 'none';
                h.msty = 'o';
            end
        end %quiver3Hy
        function h = quiver3Hz(varargin)
            if nargin == 2;
                H = varargin{1};
                in = varargin{2};
                h = inter2D.quiver3Pose( H(1:3,4),H(1:3,3), in );
            else
                h.scl = 1;
                h.lwd = 1;
                h.lcol = [0,0,1];
                h.lsty = '-';
                h.msz = 5;
                h.mcol = 'none';
                h.msty = 'o';
            end
        end %quiver3Hz
        function h = plot3H(varargin)
            if nargin == 2;
                H = varargin{1};
                style = varargin{2};
                h = line( H(1,4),H(2,4),H(3,4),...
                    'LineWidth',style.lwd, 'Color',style.lcol, 'LineStyle',style.lsty, ...
                    'MarkerSize',style.msz, 'MarkerEdgeColor',style.lcol, 'MarkerFaceColor',style.mcol, 'Marker',style.msty );
            else
                h.scl = 1;
                h.lwd = 1;
                h.lcol = [0,0,1];
                h.lsty = '-';
                h.msz = 5;
                h.mcol = 'none';
                h.msty = 'o';
            end
        end %plot3H
        function h = drawTriad( varargin )
            if nargin == 2;
                H = varargin{1};
                style = varargin{2};
                
                h = quiver3( H(1,4),H(2,4),H(3,4), H(1,1)*style.scl,H(2,1)*style.scl,H(3,1)*style.scl, 0, ...
                             'LineWidth',style.lwd, 'Color','r', 'LineStyle',style.lsty );
                set(h,'DisplayName',style.name);
                h = quiver3( H(1,4),H(2,4),H(3,4), H(1,2)*style.scl,H(2,2)*style.scl,H(3,2)*style.scl, 0, ...
                             'LineWidth',style.lwd, 'Color','g', 'LineStyle',style.lsty );
                set(h,'DisplayName',style.name);
                h = quiver3( H(1,4),H(2,4),H(3,4), H(1,3)*style.scl,H(2,3)*style.scl,H(3,3)*style.scl, 0, ...
                             'LineWidth',style.lwd, 'Color','b', 'LineStyle',style.lsty );
                set(h,'DisplayName',style.name);
                
                h = plot3( H(1,4),H(2,4),H(3,4), ...
                    'LineWidth',style.lwd, 'Color',style.lcol, 'LineStyle',style.lsty, ...
                    'MarkerSize',style.msz, 'MarkerEdgeColor',style.lcol, 'MarkerFaceColor',style.mcol, 'Marker',style.msty );
                set(h,'DisplayName',style.name);
            else
                h.scl = 1;
                h.lwd = 1;
                h.lcol = [1,0,0];
                h.lsty = '-';
                h.msz = 5;
                h.mcol = 'none';
                h.msty = 'o';
                h.name = 'unnamed';
            end
        end %drawTriad
        
        %plot array of poses
        function h = quiver3poseArray( xyzs, vecs, style )
            vecs = vecs*style.scl;
            xyzs = reshape(xyzs', 3, length(xyzs)); %transpose b/c matlab reshapes over rows then columns
            vecs = reshape(vecs', 3, length(vecs));
            h = quiver3( xyzs(1,:),xyzs(2,:),xyzs(3,:), vecs(1,:),vecs(2,:),vecs(3,:), 0,...
                    'LineWidth',style.lwd, 'Color',style.lcol, 'LineStyle',style.lsty );
            
        end %quiver3poseArray
        
        function p = perDiff( a, b )
            if (a ~= 0)
                p = (a-b)./a;
            elseif (b~=0)
                p = (b-a)./b;
            else
                p = 0;
            end
        end
        
    end %static methods
end