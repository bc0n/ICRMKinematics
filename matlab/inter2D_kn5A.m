% This class creates an inter2D with a 5 parameter kinematic model
% Frames:
% 0 - Global reference, X directed along the proximal segment (before bend), Y right, Z up
% 1 - the distal face of the coaxial input
% 2 - the pitch axis, rotating about Z2
% 3 - the proximal face of the roll, rotating about X2/3
% 4 - the catheter tip, remote center of rotation about Z4, X4 points out the catheter tip
% 5 - a projection of distance q5 along X4
classdef inter2D_kn5A < inter2D
    methods
        
        function obj = inter2D_kn5A()
            obj@inter2D();
            obj.name = 'kn5A';
            
            obj.kns.tx01 = 806;
            obj.kns.ty01 = -66;
            obj.kns.tz01 = -28;
            obj.kns.rz01 = -.24;
            obj.kns.lCath = obj.drw.lCath; %as fabricated
            
            obj.nums.kns = 5;
            obj.nums.qps = 5;
        end
        
         function varargout = forwardK(obj, qp, kn)
            if nargin < 3
                kn = obj.kns;
            end
            %if kns array, translate into struct
            if isstruct(pm)
                kns = kn;
            else
                kns = obj.knArray2Struct(pm);
            end

            % add initial joint positions...
            qp = reshape(qp, 5,1);
            qp = qp + obj.qp0;
            
            %check q bounds
            if qp(4) < 1e-3; qp(4) = 1e-3; end %the catheter actuation angle cannot be less than 1e-3 (~0)
            
            %compute transforms
            H01 = obj.Tx(kns.tx01)*obj.Ty(kns.ty01)*obj.Tz(kns.tz01)*obj.Rz(kns.rz01) * obj.Rx(qp(1)); %#ok<*PROP> %prox roll
            H12 = obj.Tx(obj.drw.lProx) * obj.Rz(qp(2)); %pitch
            H23 = obj.Tx(obj.drw.lPtch) * obj.Rx(qp(3)); %roll
            r = kns.lCath/qp(4); %radius of catheter arc
            H34 = obj.Tx(obj.drw.lRoll) * obj.Ty(r*(1-cos(qp(4))))*obj.Tx(r*sin(qp(4))) * obj.Rz(qp(4));
            H45 = obj.Tx( qp(5) ); %translation along x4 to the target
            
            varargout = {H01*H12*H23*H34*H45}; %H05
            if nargout == 5;
                varargout = [{H01}, {H01*H12}, {H01*H12*H23}, {H01*H12*H23*H34}, {H01*H12*H23*H34*H45}]; %return H01,H02,H03,H04,H05
            end
         end %fk
    end % public methods
    
    methods (Static)
        function kns = knArray2Struct(pma)
            if length(pma) == 5;
                kns.tx01 = pma(1);
                kns.ty01 = pma(2);
                kns.tz01 = pma(3);
                kns.rz01 = pma(4);
                kns.lCath = pma(5);
            else
                error('Matlab:inter2D_kn5A','pma does not have 5 elements');
            end
        end %knArray2Struct
        function pma = knStruct2Array(kns)
            if isstruct(kns)
                pma(1,1) = kns.tx01;
                pma(2,1) = kns.ty01;
                pma(3,1) = kns.tz01;
                pma(4,1) = kns.rz01;
                pma(5,1) = kns.lCath;
            else
                error('Matlab:inter2D_kn5A','kns is not a struct');
            end
        end %knStruct2Array
        
    end %statics
    
end%class