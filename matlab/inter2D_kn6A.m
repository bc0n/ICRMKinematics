% This class creates an inter2D model and provides for generating trajectories, solving the kinematics, and plotting the results
% Frames:
% 0 - Global reference, X directed along proximal segment, Z upward
% 1 - 
% 2 - 
% 3 -   
% 4 - 
% 5 - Translation along the X4
classdef inter2D_kn6A < inter2D
%     properties (Access = public)
%         pms;        
%     end %public props
    
    methods
        
        function obj = inter2D_kn6A()
            obj@inter2D();
            obj.name = 'kn6A';
            
            obj.pms.tx01 = 806;
            obj.pms.ty01 = -66;
            obj.pms.tz01 = -28;
            obj.pms.rz01 = -.24;
            obj.pms.tx23 = 8.2;
            obj.pms.cathL = 98;%95+3.75; %sldrw
            
            obj.nums.pms = 6;
            obj.nums.qps = 5;
        end
        
         function varargout = forwardK(obj, qp, pm)
            if nargin < 3
                pm = obj.pms;
            end
            %if pms array, translate into struct
            if isstruct(pm)
                pms = pm;
            else
                pms = obj.paramArray2Struct(pm);
            end

            
            %check q bounds
            if qp(4) < 1e-3; qp(4) = 1e-3; end %the catheter actuation angle cannot be less than 1e-3 (~0)
            
            %compute transforms
            H01 = obj.Tx(pms.tx01)*obj.Ty(pms.ty01)*obj.Tz(pms.tz01)*obj.Rz(pms.rz01) * obj.Rx(qp(1)); %#ok<*PROP> %prox roll
            H12 = obj.Rz(qp(2)); %pitch
            H23 = obj.Tx(pms.tx23) * obj.Rx(qp(3)); %roll
            r = pms.cathL/qp(4); %radius of catheter arc
            H34 = obj.Ty(r*(1-cos(qp(4)))) * obj.Tx(r*sin(qp(4))) * obj.Rz(qp(4));
            H45 = obj.Tx( qp(5) ); %translation along x4 to the target
            
            varargout = {H01*H12*H23*H34*H45}; %H05
            if nargout == 5;
                varargout = [{H01}, {H01*H12}, {H01*H12*H23}, {H01*H12*H23*H34}, {H01*H12*H23*H34*H45}]; %return H01,H02,H03,H04,H05
            end
         end %fk
    end % public methods
    
    methods (Static)
        function pms = paramArray2Struct(pma)
            if length(pma) == 6;
                pms.tx01 = pma(1);
                pms.ty01 = pma(2);
                pms.tz01 = pma(3);
                pms.rz01 = pma(4);
                pms.tx23 = pma(5);
                pms.cathL = pma(6);
            else
                error('Matlab:inter2D_kn6A','pma does not have 6 elements');
            end
        end %paramArray2Struct
        function pma = paramStruct2Array(pms)
            if isstruct(pms)
                pma(1,1) = pms.tx01;
                pma(2,1) = pms.ty01;
                pma(3,1) = pms.tz01;
                pma(4,1) = pms.rz01;
                pma(5,1) = pms.tx23;
                pma(6,1) = pms.cathL;
            else
                error('Matlab:inter2D_kn6A','pms is not a struct');
            end
        end %paramStruct2Array
        
    end %statics
    
end%class