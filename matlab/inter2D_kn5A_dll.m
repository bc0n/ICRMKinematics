% This class creates an inter2D that calls the 5A methods from ICRMKinematics.dll
classdef inter2D_kn5A_dll < inter2D
    
    
    methods
        
        function obj = inter2D_kn5A_dll(varargin)
            obj@inter2D();
            obj.name = 'kn5Adll';
            
            obj.kns.tx01 = 806;
            obj.kns.ty01 = -66;
            obj.kns.tz01 = -28;
            obj.kns.rz01 = -.24;
            obj.kns.lCath = obj.drw.lCath; %as fabricated
            
            obj.nums.kns = 5;
            obj.nums.qps = 5;
            
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
            if nargin < 3
                kn = obj.kns;
            end
            %if kns struct, translate into array
            if isstruct(kn)
                kns = obj.knStruct2Array(kn);
            else
                kns = kn;
            end
            if nargout == 1; %only H05
                [~,~,~,Hc] = calllib( obj.name, 'get5AH05', qps, kns, zeros(12,1));
                varargout{1} = [reshape(Hc,3,4);0,0,0,1];
            else %all
                [~,~,~,H] = calllib( obj.name, 'get5AH01', qps, kns, zeros(12,1) );
                varargout{1} = [reshape(H,3,4);0,0,0,1];
                [~,~,~,H] = calllib( obj.name, 'get5AH02', qps, kns, zeros(12,1) );
                varargout{2} = [reshape(H,3,4);0,0,0,1];
                [~,~,~,H] = calllib( obj.name, 'get5AH03', qps, kns, zeros(12,1) );
                varargout{3} = [reshape(H,3,4);0,0,0,1];
                [~,~,~,H] = calllib( obj.name, 'get5AH04', qps, kns, zeros(12,1) );
                varargout{4} = [reshape(H,3,4);0,0,0,1];
                [~,~,~,H] = calllib( obj.name, 'get5AH05', qps, kns, zeros(12,1) );
                varargout{5} = [reshape(H,3,4);0,0,0,1];
            end
        end %FK
        
        % res[telapsed,ret,kn0,qp0,fmin] = obj.estimate_qp0kn0_xyzuxuyuz( qps, Hs, kn0, knup,kndn, qp0,qpup,qpdn )
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
            [ret,~,~,~,knEst,~,~,qpEst,~,~,fmin] = calllib(obj.name, 'estimate_qp0kn0_xyzuxuyuz5A', n,   stackedQ, stackedX, stackedQ, kn0, knup, kndn, qp0, [qpdn;qpup], nlArray, 1e3);
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
            if length(kna) == 5;
                kns.tx01 = kna(1);
                kns.ty01 = kna(2);
                kns.tz01 = kna(3);
                kns.rz01 = kna(4);
                kns.lCath = kna(5);
            else
                error('Matlab:inter2D_kn5A_dll','kna does not have 5 elements');
            end
        end
        function kna = knStruct2Array(kns)
            if isstruct(kns)
                kna(1,1) = kns.tx01;
                kna(2,1) = kns.ty01;
                kna(3,1) = kns.tz01;
                kna(4,1) = kns.rz01;
                kna(5,1) = kns.lCath;
            else
                error('Matlab:inter2D_kn5A_dll','kns is not a struct');
            end
        end %paramStruct2Array
        
    end %statics
    
    
end