% This class creates an inter2D with an 6 parameter kinematic model
% Frames:
% 0 - Global reference, X directed along proximal segment, Z upward
% 1 - 
% 2 - 
% 3 -   
% 4 - 
% 5 - Translation along the X4
classdef inter2D_kn6A < inter2D

    methods
        
        function obj = inter2D_kn6A()
            obj@inter2D();
            obj.name = 'kn6A';
            
            obj.kns.tx01 = 806;
            obj.kns.ty01 = -66;
            obj.kns.tz01 = -28;
            obj.kns.rz01 = -.24;
            obj.kns.ry34 = 0; %angle of the articulation plane, to account for drooping
            obj.kns.lCath = obj.drw.lCath; %as fabricated
            
            obj.nums.kns = 6;
            obj.nums.qps = 5;
        end
        
         function varargout = forwardK(obj, qp, kn)
            if nargin < 3
                kn = obj.kns;
            end
            %if kns array, translate into struct
            if isstruct(kn)
                kns = kn;
            else
                kns = obj.knArray2Struct(kn);
            end
            
            %check q bounds
            if qp(4) < 1e-3; qp(4) = 1e-3; end %the catheter actuation angle cannot be less than 1e-3 (~0)
            
            %compute transforms
            H01 = obj.Tx(kns.tx01)*obj.Ty(kns.ty01)*obj.Tz(kns.tz01)*obj.Rz(kns.rz01) * obj.Rx(qp(1)); %#ok<*PROP> %prox roll
            H12 = obj.Tx(obj.drw.lProx) * obj.Rz(qp(2)); %pitch
            H23 = obj.Tx(obj.drw.lPtch) * obj.Rx(qp(3)); %roll
            r = kns.lCath/qp(4); %radius of catheter arc
            H34 = obj.Tx(obj.drw.lRoll)*obj.Ry(kns.ry34) * obj.Ty(r*(1-cos(qp(4))))*obj.Tx(r*sin(qp(4))) * obj.Rz(qp(4)); %ry34 for out-of-articulation plane
            H45 = obj.Tx( qp(5) ); %translation along x4 to the target
            
            varargout = {H01*H12*H23*H34*H45}; %H05
            if nargout == 5;
                varargout = [{H01}, {H01*H12}, {H01*H12*H23}, {H01*H12*H23*H34}, {H01*H12*H23*H34*H45}]; %return H01,H02,H03,H04,H05
            end
         end %fk
    end % public methods
    
    methods (Static)
        function kns = knArray2Struct(kna)
            if length(kna) == 6;
                kns.tx01 = kna(1);
                kns.ty01 = kna(2);
                kns.tz01 = kna(3);
                kns.rz01 = kna(4);
                kns.ry34 = kna(5);
                kns.lCath = kna(6);
            else
                error('Matlab:inter2D_kn6A','kna does not have 6 elements');
            end
        end %knArray2Struct
        function kna = knStruct2Array(kns)
            if isstruct(kns)
                kna(1,1) = kns.tx01;
                kna(2,1) = kns.ty01;
                kna(3,1) = kns.tz01;
                kna(4,1) = kns.rz01;
                kna(5,1) = kns.ry34;
                kna(6,1) = kns.lCath;
            else
                error('Matlab:inter2D_kn6A','kns is not a struct');
            end
        end %paramStruct2Array
        
    end %statics
    
end%class
