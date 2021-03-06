% This class creates an inter2D that calls the 11A methods from ICRMKinematics.dll
classdef inter2D_kn11A_dll < inter2D
    
    
    methods
        
        function obj = inter2D_kn11A_dll(varargin)
            obj@inter2D();
            obj.name = 'kn11Adll';
            
            obj.kns.tx01 = 806;
            obj.kns.ty01 = -66;
            obj.kns.tz01 = -28;
            obj.kns.ry01 = 0;
            obj.kns.rz01 = -.24;
            obj.kns.ry34 = 0; %rotation at the catheter base
            obj.kns.rz34 = 0; %redundant with alpha near 0
            obj.kns.kAlpha = 1;%alpha gain
            obj.kns.eAlpha = 1;%alpha raised to
            obj.kns.lCath = obj.drw.lCath; %as fabricated
            obj.kns.ry45 = 0; %rotation at the catheter tip
            
            obj.nums.kns = 11;
            obj.nums.qps = 5;
            
            % Algs
            %          00;  GN_DIRECT
            %          01;  GN_DIRECT_L --locally biased
            %          02;  GN_DIRECT_L_NOSCAL
            %          03;  GN_DIRECT_L_RAND
            %          04;  GN_DIRECT_L_RAND_NOSCAL
            %          05;  GN_DIRECT_NOSCAL
            %          06;  GN_ESCH
            %          07;  GN_ISRES
            %          08;  GN_MLSL -- slow due to local searches
            %          09;  GN_MLSL_LDS
            %          10;  GN_ORIG_DIRECT
            %          11;  GN_ORIG_DIRECT_L
            %          12;  LN_BOBYQA
            %          13;  LN_COBYLA
            %          14;  LN_NelderMead
            %          17;  LN_PRAXIS
            %          18;  LN_SUBPLX
            %   returns -1 = failure
            %           -2 = invalid args
            %           -3 = out of memory
            %           -4 = roundoff limited
            %           -5 = forced stop
            %            1 = success
            %            2 = fun stopval reached
            %            3 = fun tol reached
            %            4 = x tol reached
            %            5 = max evals
            %            6 = max time reached
            obj.opts = obj.getDefaultNLOptParams; %replace fminunc options with NLOpt options
            
            %Load Library
            if nargin == 2;
                dpath = varargin{1};
                hpath = varargin{2};
            else
                dpath = 'ICRMKinematics.dll';
                hpath = 'kinematicsDLL.h';
            end
            %if : "the specified module could not be found" then it can't find libnlopt.dll; copy into release
            if ~libisloaded(obj.name);
                loadlibrary(dpath,hpath, 'alias', obj.name);
            end
            %get processed signatures with % libfunctions('kn11Adll','-full')
        end
        
        function delete(obj)
            unloadlibrary(obj.name);
        end
        
        function varargout = forwardK(obj, qps, kn)
            if nargin < 3
                kn = obj.kns;
            end
            %if kns struct, translate into array
            if isstruct(kn)
                kns = obj.knStruct2Array(kn);
            else
                kns = kn;
            end

            [~,~,~,H] = calllib( obj.name, 'get11AH05', qps, kns, zeros(12,1));
            varargout = {[reshape(H,3,4);0,0,0,1]};
            if nargout > 1;
                [~,~,~,H] = calllib( obj.name, 'get11AH01', qps, kns, zeros(12,1) );
                H01 = [reshape(H,3,4);0,0,0,1];
                [~,~,~,H] = calllib( obj.name, 'get11AH02', qps, kns, zeros(12,1) );
                H02 = [reshape(H,3,4);0,0,0,1];
                [~,~,~,H] = calllib( obj.name, 'get11AH03', qps, kns, zeros(12,1) );
                H03 = [reshape(H,3,4);0,0,0,1];
                [~,~,~,H] = calllib( obj.name, 'get11AH04', qps, kns, zeros(12,1) );
                H04 = [reshape(H,3,4);0,0,0,1];
                [~,~,~,H] = calllib( obj.name, 'get11AH05', qps, kns, zeros(12,1) );
                H05 = [reshape(H,3,4);0,0,0,1];
                varargout = [{H01}, {H02}, {H03}, {H04}, {H05}];
            end
        end %FK
        
        
        function res = fun_kn0_xyzuxuyuz(obj, qps, Hs, kn0)
            n = numel(qps)/obj.nums.qps;
            
            stackedQ = reshape(qps', n*obj.nums.qps,1);
            stackedX = reshape(Hs(:,1:3,4)', n*3,1);
            stackedU = reshape(Hs(:,1:3,1)', n*3,1);
        
            if isstruct(kn0); kn0 = obj.knStruct2Array(kn0); end;
            
            fmin = 1111111.1;
            
            tic
            [ret,~,~,~,~,fmin] = calllib(obj.name, 'fun_kn0_xyzuxuyuz11A', n, stackedQ, stackedX, stackedU, kn0, fmin);
            res.telapsed = toc;
            res.ret = ret;
            res.fmin = fmin;
        end
        
        % res[ret,qp0,fmin] = obj.estimate_qp0_xyzuxuyuz( qps, Hs, kn0, qp0,qpup,qpdn )
        function res = estimate_qp0_xyzuxuyuz(obj, qps, Hs, kn0, qp0,qpup,qpdn )
            n = numel(qps)/obj.nums.qps;
            
            stackedQ = reshape(qps', n*obj.nums.qps,1);
            stackedX = reshape(Hs(:,1:3,4)', n*3,1);
            stackedU = reshape(Hs(:,1:3,1)', n*3,1);
            qpup = reshape(qpup,5,1); qpdn = reshape(qpdn,5,1);
            if isstruct(kn0); kn0 = obj.knStruct2Array(kn0); end;
            nlArray = obj.nlStruct2Array(obj.opts);
            
            tic
            [ret,~,~,~,knEst,qpEst,~,~,~,fmin] = calllib(obj.name, 'estimate_qp0_xyzuxuyuz11A', n, stackedQ, stackedX, stackedU, kn0, qp0,qpup,qpdn, nlArray, 1e3);
            res.telapsed = toc;
            res.ret = ret;
            res.kn0 = obj.knArray2Struct(knEst);
            res.qp0 = qpEst;
            res.fmin = fmin;
        end
        
        % res[ret,kn0,fmin] = obj.estimate_kn0_xyzuxuyuz( qps, Hs, kn0,knup,kndn )
        function res = estimate_kn0_xyzuxuyuz(obj, qps, Hs, kn0,knup,kndn )
            n = numel(qps)/obj.nums.qps;
            
            stackedQ = reshape(qps', n*obj.nums.qps,1);
            stackedX = reshape(Hs(:,1:3,4)', n*3,1);
            stackedU = reshape(Hs(:,1:3,1)', n*3,1);
            if isstruct(kn0); kn0 = obj.knStruct2Array(kn0); end;
            if isstruct(knup); knup = obj.knStruct2Array(knup); end;
            if isstruct(kndn); kndn = obj.knStruct2Array(kndn); end;
            
            nlArray = obj.nlStruct2Array(obj.opts);
            fmin = 1111111.1;
            
            tic
            [ret,~,~,~,knEst,~,~,~,fmin] = calllib(obj.name, 'estimate_kn0_xyzuxuyuz11A', n, stackedQ, stackedX, stackedU, kn0, knup, kndn, nlArray, fmin);
            res.telapsed = toc;
            res.ret = ret;
            res.kn0 = obj.knArray2Struct(knEst);
            res.fmin = fmin;
        end
        % res[ret,kn0,fmin] = obj.estimate_kn0_xyzuxuyuz( qps, Hs, kn0,knup,kndn )
        function res = estimate_kn0_xyzuxuyuz_sub(obj, qps, Hs, kn0,knup,kndn, subset )
            n = numel(qps)/obj.nums.qps;
            
            stackedQ = reshape(qps', n*obj.nums.qps,1);
            stackedX = reshape(Hs(:,1:3,4)', n*3,1);
            stackedU = reshape(Hs(:,1:3,1)', n*3,1);
            if isstruct(kn0); kn0 = obj.knStruct2Array(kn0); end;
            if isstruct(knup); knup = obj.knStruct2Array(knup); end;
            if isstruct(kndn); kndn = obj.knStruct2Array(kndn); end;
            
            nlArray = obj.nlStruct2Array(obj.opts);
            fmin = 1111111.1;
            
            tic
            [ret,~,~,~,knEst,~,~,~,~,fmin] = calllib(obj.name, 'estimate_kn0_xyzuxuyuz11A_subset', n, stackedQ, stackedX, stackedU, kn0, knup, kndn, subset, nlArray, fmin);
            res.telapsed = toc;
            res.ret = ret;
            res.kn0 = obj.knArray2Struct(knEst);
            res.fmin = fmin;
        end
        
        % res[ret,kn0,qp0,fmin] = obj.estimate_qp0kn0_xyzuxuyuz( qps, Hs, kn0, knup,kndn, qp0,qpup,qpdn )
        function res = estimate_qp0kn0_xyzuxuyuz(obj, qps, Hs, qp0,qpup,qpdn, kn0, knup,kndn )
            n = numel(qps)/obj.nums.qps;
            
            stackedQ = reshape(qps', n*obj.nums.qps,1);
            stackedX = reshape(Hs(:,1:3,4)', n*3,1);
            stackedU = reshape(Hs(:,1:3,1)', n*3,1);
            qpup = reshape(qpup,5,1); qpdn = reshape(qpdn,5,1);
            if isstruct(kn0); kn0 = obj.knStruct2Array(kn0); end;
            if isstruct(knup); knup = obj.knStruct2Array(knup); end;
            if isstruct(kndn); kndn = obj.knStruct2Array(kndn); end;
            
            nlArray = obj.nlStruct2Array(obj.opts);
            
            tic
            % [int32, doublePtr, doublePtr, doublePtr, doublePtr, doublePtr, doublePtr, doublePtr, doublePtr, doublePtr, doublePtr] estimate_qp0kn0_xyzuxuyuz11A(                   int32,  doublePtr, doublePtr, doublePtr, doublePtr,   doublePtr, doublePtr, doublePtr, doublePtr, doublePtr, doublePtr)
            [    ret,         ~,         ~,         ~,         d,         ~,         f,         ~,         ~,         ~,         j] = calllib(obj.name, 'estimate_qp0kn0_xyzuxuyuz11A', n,   stackedQ,  stackedX,  stackedU,       qp0, [qpdn;qpup],       kn0,      knup,      kndn,   nlArray, 1e3);
            %    ret    stackedQ  stackedX   stackedU        qp0     qpupdn        kn0       knup       kndn        nla       fmin
            res.telapsed = toc;
            res.ret = ret;
            res.kn0 = obj.knArray2Struct(f);
            res.qp0 = d;
            res.fmin = j;
        end
        function res = estimate_qp0kn0_xyzuxuyuz_subset(obj, qps, Hs, qp0,qpup,qpdn, kn0,knup,kndn,subset)
            n = numel(qps)/obj.nums.qps;
            
            stackedQ = reshape(qps', n*obj.nums.qps,1);
            stackedX = reshape(Hs(:,1:3,4)', n*3,1);
            stackedU = reshape(Hs(:,1:3,1)', n*3,1);
            qpup = reshape(qpup,5,1); qpdn = reshape(qpdn,5,1);
            if isstruct(kn0); kn0 = obj.knStruct2Array(kn0); end;
            if isstruct(knup); knup = obj.knStruct2Array(knup); end;
            if isstruct(kndn); kndn = obj.knStruct2Array(kndn); end;
            nlArray = obj.nlStruct2Array(obj.opts);
            
            tic
            % [int, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl] estimate_qp0kn0_xyzuxuyuz11A_subset(                     int,      dbl,      dbl,      dbl, dbl,         dbl, dbl,  dbl,  dbl,    dbl,    dbl, dbl)
            [  ret,   ~,   ~,   ~,   d,   ~,   f,   ~,   ~,   ~,   ~,   k] = calllib(obj.name, 'estimate_qp0kn0_xyzuxuyuz11A_subset', n, stackedQ, stackedX, stackedU, qp0, [qpdn;qpup], kn0, knup, kndn, subset, nlArray, 1e3);
            %  ret    q    x    u  qp0 updn  kn0 knup kndn  sub  nla fmin
            
            res.telapsed = toc;
            res.ret = ret;
            res.kn0 = obj.knArray2Struct(f);
            res.qp0 = d;
            res.fmin = k;
        end
        function res = estimate_qp0kn0_xyzuxuyuz_subset_filter(obj, qps, Hs, qp0,qpup,qpdn, kn0,knup,kndn,subset, fcutHz )
            n = numel(qps)/obj.nums.qps;
sty = bplot3; sty.lco = [0,0,0]; sty.lst = 'none'; sty.msz = 10;
bplot3(Hs(:,1:3,4), sty); hold on; title('filtering');
            % Filter ASC
            fs = 1/10e-3; %10ms loop
            wn = fcutHz/fs/2;
            [z,p,k] = butter(3, wn, 'low');
            [sos,g] = zp2sos(z,p,k);
            for i = 3:-1:1;
                for j = 4:-1:1;
                    Hs(:,i,j) = filtfilt(sos,g, Hs(:,i,j));
                end
            end
sty.lco = [1,0,0]; bplot3(Hs(:,1:3,4), sty);

            stackedQ = reshape(qps', n*obj.nums.qps,1);
            stackedX = reshape(Hs(:,1:3,4)', n*3,1);
            stackedU = reshape(Hs(:,1:3,1)', n*3,1);
            qpup = reshape(qpup,5,1); qpdn = reshape(qpdn,5,1);
            if isstruct(kn0); kn0 = obj.knStruct2Array(kn0); end;
            if isstruct(knup); knup = obj.knStruct2Array(knup); end;
            if isstruct(kndn); kndn = obj.knStruct2Array(kndn); end;
            nlArray = obj.nlStruct2Array(obj.opts);

            tic
            % [int, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl] estimate_qp0kn0_xyzuxuyuz11A_subset(                     int,      dbl,      dbl,      dbl, dbl,         dbl, dbl,  dbl,  dbl,    dbl,    dbl, dbl)
            [  ret,   ~,   ~,   ~,   d,   ~,   f,   ~,   ~,   ~,   ~,   k] = calllib(obj.name, 'estimate_qp0kn0_xyzuxuyuz11A_subset', n, stackedQ, stackedX, stackedU, qp0, [qpdn;qpup], kn0, knup, kndn, subset, nlArray, 1e3);
            %  ret    q    x    u  qp0 updn  kn0 knup kndn  sub  nla fmin
            
            res.telapsed = toc;
            res.ret = ret;
            res.kn0 = obj.knArray2Struct(f);
            res.qp0 = d;
            res.fmin = k;
        end
        function res = estimate_qp0kn0_xyzuxuyuz_subset_mlsl(obj, qps, Hs, qp0,qpup,qpdn, kn0,knup,kndn,subset, nlMLSL,nlLDS )
            n = numel(qps)/obj.nums.qps;
            
            stackedQ = reshape(qps', n*obj.nums.qps,1);
            stackedX = reshape(Hs(:,1:3,4)', n*3,1);
            stackedU = reshape(Hs(:,1:3,1)', n*3,1);
            qpup = reshape(qpup,5,1); qpdn = reshape(qpdn,5,1);
            if isstruct(kn0); kn0 = obj.knStruct2Array(kn0); end;
            if isstruct(knup); knup = obj.knStruct2Array(knup); end;
            if isstruct(kndn); kndn = obj.knStruct2Array(kndn); end;
            if isstruct(nlMLSL); nlMLSL = obj.nlStruct2Array(nlMLSL); end;
            if isstruct(nlLDS); nlLDS = obj.nlStruct2Array(nlLDS); end;
            tic
            % [int, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl, dbl] estimate_qp0kn0_xyzuxuyuz11A_subset_mlsllds(                     int,      dbl,      dbl,      dbl, dbl,         dbl, dbl,  dbl,  dbl,    dbl,    dbl,   dbl, dbl)
            [  ret,   ~,   ~,   ~,   d,   ~,   f,   ~,   ~,   ~,   ~,   ~,   l] = calllib(obj.name, 'estimate_qp0kn0_xyzuxuyuz11A_subset_mlsllds', n, stackedQ, stackedX, stackedU, qp0, [qpdn;qpup], kn0, knup, kndn, subset, nlMLSL, nlLDS, 1e3);
            %  ret    q    x    u  qp0 updn  kn0 knup kndn  sub  nla lnla fmin
            res.telapsed = toc;
            res.ret = ret;
            res.kn0 = obj.knArray2Struct(f);
            res.qp0 = d;
            res.fmin = l;
        end
        
        function h = drawConfig(obj, qps, lcol )
            if nargin < 3;
                lcol = [1,1,1];
            end
            
            [h01,h02,h03,h04,h05] = obj.forwardK(qps);
            col = [0;0;0;1]; %H*col = XYZ1
            
            %01 too large
            %12
            xyz = [ h01*col, h01*obj.Ty( obj.drw.rProx)*col, h01*obj.Ty( obj.drw.rProx)*obj.Tx(obj.drw.lProx)*col, ...
                    h01*obj.Ty(-obj.drw.rProx)*obj.Tx(obj.drw.lProx)*col, h01*obj.Ty(-obj.drw.rProx)*col, h01*col,h02*col];
            hs(1) = plot3(xyz(1,:), xyz(2,:), xyz(3,:), '-', 'color', lcol*.7);
            xyz = [ h01*col, h01*obj.Tz( obj.drw.rProx)*col, h01*obj.Tz( obj.drw.rProx)*obj.Tx(obj.drw.lProx)*col, ...
                    h01*obj.Tz(-obj.drw.rProx)*obj.Tx(obj.drw.lProx)*col, h01*obj.Tz(-obj.drw.rProx)*col, h01*col];
            hs(end+1) = plot3(xyz(1,:), xyz(2,:), xyz(3,:), '-', 'color', lcol*.7);
            %23
            xyz = [ h02*col, h02*obj.Ty( obj.drw.rProx)*col, h02*obj.Ty( obj.drw.rProx)*obj.Tx(obj.drw.lPtch)*col, ...
                    h02*obj.Ty(-obj.drw.rProx)*obj.Tx(obj.drw.lPtch)*col, h02*obj.Ty(-obj.drw.rProx)*col, h02*col,h03*col];
            hs(end+1) = plot3(xyz(1,:), xyz(2,:), xyz(3,:), '-', 'color', lcol*.5);
            xyz = [ h02*col, h02*obj.Tz( obj.drw.rProx)*col, h02*obj.Tz( obj.drw.rProx)*obj.Tx(obj.drw.lPtch)*col, ...
                    h02*obj.Tz(-obj.drw.rProx)*obj.Tx(obj.drw.lPtch)*col, h02*obj.Tz(-obj.drw.rProx)*col, h02*col];
            hs(end+1) = plot3(xyz(1,:), xyz(2,:), xyz(3,:), '-', 'color', lcol*.5);
            
            %34
            qps(4) = qps(4) + obj.qp0(4);
%             if qps(4) < 1e-3; qps(4) = 1e-3; end;
%             r = obj.kns.lCath/qps(4); %radius of catheter arc
            al = (qps(4)*obj.kns.kAlpha)^obj.kns.eAlpha; %allow slope,exponent in the articulation angle
            if al < 1e-3; al = 1e-3; end
            r = obj.kns.lCath/al; %radius of catheter arc
            
            h3 = h03*obj.Tx(obj.drw.lRoll)*obj.Ry(obj.kns.ry34)*obj.Rz(obj.kns.rz34);
            cathPts = [ h3*col, h03*obj.Tz(obj.drw.rCath)*obj.Tx(obj.drw.lRoll)*col, h03*obj.Tz(obj.drw.rCath)*col,...
                        h03*col,h03*obj.Ty(obj.drw.rCath)*col, h03*obj.Ty(obj.drw.rCath)*obj.Tx(obj.drw.lRoll)*col, h3*col, ...
                        h3* obj.Ty(r*(1-cos(al*.1)))*obj.Tx(r*sin(al*.1))*obj.Rz(al*.1)*col, ...
                        h3* obj.Ty(r*(1-cos(al*.2)))*obj.Tx(r*sin(al*.2))*obj.Rz(al*.2)*col, ...
                        h3* obj.Ty(r*(1-cos(al*.3)))*obj.Tx(r*sin(al*.3))*obj.Rz(al*.3)*col, ...
                        h3* obj.Ty(r*(1-cos(al*.4)))*obj.Tx(r*sin(al*.4))*obj.Rz(al*.4)*col, ...
                        h3* obj.Ty(r*(1-cos(al*.5)))*obj.Tx(r*sin(al*.5))*obj.Rz(al*.5)*col, ...
                        h3* obj.Ty(r*(1-cos(al*.6)))*obj.Tx(r*sin(al*.6))*obj.Rz(al*.6)*col, ...
                        h3* obj.Ty(r*(1-cos(al*.7)))*obj.Tx(r*sin(al*.7))*obj.Rz(al*.7)*col, ...
                        h3* obj.Ty(r*(1-cos(al*.8)))*obj.Tx(r*sin(al*.8))*obj.Rz(al*.8)*col, ...
                        h3* obj.Ty(r*(1-cos(al*.9)))*obj.Tx(r*sin(al*.9))*obj.Rz(al*.9)*col, ...
                        h3* obj.Ty(r*(1-cos(al*1.)))*obj.Tx(r*sin(al*1.))*obj.Rz(al*1.)*col ];
            hs(end+1) = plot3( cathPts(1,:),cathPts(2,:),cathPts(3,:), '.-', 'color', lcol*1);
            
            %45
            xyz = [ h04*col, h05*col];
            hs(end+1) = plot3(xyz(1,:), xyz(2,:), xyz(3,:), '-', 'color', lcol*.7);
            
            hs(1).DisplayName = 'drawCath';
            for i = 2:length(hs); hs(i).HandleVisibility = 'off'; end
            he = linkprop(hs, 'Visible');
            h = hs(1);
        end

        
    end %public methods
    
    methods (Static)
        function nlpms = getDefaultNLOptParams
            nlpms.maxIts = 1e9;
            nlpms.maxTime = 60;
            nlpms.method = 13; %ln_cobyla
            nlpms.minFunVal = 1e-6; %desired minimum function value
            nlpms.tolFun = 1e-6; %tolerance to ensure continued decrease in the function value
            nlpms.tolX = 1e-6; %tolerance on the smallest step
        end
        function nls = nlArray2Struct(nla)
            nls.maxIts =    nla(1);
            nls.maxTime =   nla(2);
            nls.method =    nla(3);
            nls.minFunVal = nla(4);
            nls.tolFun =    nla(5);
            nls.tolX =      nla(6);
        end
        function nla = nlStruct2Array(nls)
            nla(1) = nls.maxIts;
            nla(2) = nls.maxTime;
            nla(3) = nls.method;
            nla(4) = nls.minFunVal;
            nla(5) = nls.tolFun;
            nla(6) = nls.tolX;
        end
        
        function kns = knArray2Struct(kna)
            if length(kna) == 11;
                kns.tx01 = kna(1);
                kns.ty01 = kna(2);
                kns.tz01 = kna(3);
                kns.ry01 = kna(4);
                kns.rz01 = kna(5);
                kns.ry34 = kna(6);
                kns.rz34 = kna(7);
                kns.kAlpha = kna(8);
                kns.eAlpha = kna(9);
                kns.lCath = kna(10);
                kns.ry45 = kna(11);
            else
                error('Matlab:inter2D_kn11A_dll','kna does not have 11 elements');
            end
        end
        function kna = knStruct2Array(kns)
            if isstruct(kns)
                kna(1,1) = kns.tx01;
                kna(2,1) = kns.ty01;
                kna(3,1) = kns.tz01;
                kna(4,1) = kns.ry01;
                kna(5,1) = kns.rz01;
                kna(6,1) = kns.ry34;
                kna(7,1) = kns.rz34;
                kna(8,1) = kns.kAlpha;
                kna(9,1) = kns.eAlpha;
                kna(10,1) = kns.lCath;
                kna(11,1) = kns.ry45;
            else
                error('Matlab:inter2D_kn11A_dll','kns is not a struct');
            end
        end %paramStruct2Array
        
    end %statics
    
    
end