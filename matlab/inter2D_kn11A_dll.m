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
            %          03;  GN_DIRECT_L_RAND
            %          04;  GN_ESCH
            %          05;  GN_ISRES
            %          06;  GN_MLSL -- slow due to local searches
            %          07;  GN_MLSL_LDS
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
            %get processed signatures with % libfunctions('i2dll','-full')
        end
        
        function delete(obj)
            unloadlibrary(obj.name);
        end
        
        function varargout = forwardK(obj, qps, kn)
            disp(nargout)
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
        
        % res[ret,qp0,fmin] = obj.estimate_qp0_xyzuxuyuz( qps, Hs, kn0, qp0,qpup,qpdn )
        function res = estimate_qp0_xyzuxuyuz(obj, qps, Hs, kn0, qp0,qpup,qpdn )
            n = numel(qps)/obj.nums.qps;
            
            stackedQ = reshape(qps, n*obj.nums.qps,1);
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
            
            stackedQ = reshape(qps, n*obj.nums.qps,1);
            stackedX = reshape(Hs(:,1:3,4)', n*3,1);
            stackedU = reshape(Hs(:,1:3,1)', n*3,1);
            if isstruct(kn0); kn0 = obj.knStruct2Array(kn0); end;
            if isstruct(knup); knup = obj.knStruct2Array(knup); end;
            if isstruct(kndn); kndn = obj.knStruct2Array(kndn); end;
            
            nlArray = obj.nlStruct2Array(obj.opts);
            
            tic
            [ret,~,~,~,knEst,a,b,c,fmin] = calllib(obj.name, 'estimate_kn0_xyzuxuyuz11A', n, stackedQ, stackedX, stackedU, kn0, knup, kndn, nlArray, 1e3)
            res.telapsed = toc;
            res.ret = ret;
            res.kn0 = obj.knArray2Struct(knEst);
            res.fmin = fmin;
        end
        
        % res[ret,kn0,qp0,fmin] = obj.estimate_qp0kn0_xyzuxuyuz( qps, Hs, kn0, knup,kndn, qp0,qpup,qpdn )
        function res = estimate_qp0kn0_xyzuxuyuz(obj, qps, Hs, kn0, knup,kndn, qp0,qpup,qpdn )
            n = numel(qps)/obj.nums.qps;
            
            stackedQ = reshape(qps, n*obj.nums.qps,1);
            stackedX = reshape(Hs(:,1:3,4)', n*3,1);
            stackedU = reshape(Hs(:,1:3,1)', n*3,1);
            qpup = reshape(qpup,5,1); qpdn = reshape(qpdn,5,1);
            if isstruct(kn0); kn0 = obj.knStruct2Array(kn0); end;
            if isstruct(knup); knup = obj.knStruct2Array(knup); end;
            if isstruct(kndn); kndn = obj.knStruct2Array(kndn); end;
            
            nlArray = obj.nlStruct2Array(obj.opts);
            
            tic
            [ret,~,~,~,knEst,~,~,qpEst,~,~,fmin] = calllib(obj.name, 'estimate_qp0kn0_xyzuxuyuz11A', n,   stackedQ, stackedX, stackedU, kn0, knup, kndn, qp0, [qpdn;qpup], nlArray, 1e3);
            res.telapsed = toc;
            res.ret = ret;
            res.kn0 = obj.knArray2Struct(knEst);
            res.qp0 = qpEst;
            res.fmin = fmin;
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