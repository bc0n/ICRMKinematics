% This class creates an inter2D with an 11 parameter kinematic model
% Frames:
% 0 - Global reference, X directed along proximal segment, Z upward
% 1 - 
% 2 - 
% 3 -   
% 4 - 
% 5 - Translation along the X4
classdef inter2D_kn11A < inter2D
%     properties (Access = public)
%         kns; 
%     end %public props
    
    methods (Access = public)
        
        function obj = inter2D_kn11A()
            obj@inter2D();
            obj.name = 'kn11A';
            
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
            
            %don't divide by 0
            al = (qp(4)*obj.kns.kAlpha)^obj.kns.eAlpha; %allow slope,exponent in the articulation angle
            if al < 1e-3; al = 1e-3; end
            r = kns.lCath/al; %radius of catheter arc
            
            %compute transforms
            H01 = obj.Tx(kns.tx01)*obj.Ty(kns.ty01)*obj.Tz(kns.tz01)*obj.Ry(kns.ry01)*obj.Rz(kns.rz01) * obj.Rx(qp(1)); %prox roll
            H12 = obj.Tx(obj.drw.lProx) * obj.Rz(qp(2)); %pitch
            H23 = obj.Tx(obj.drw.lPtch) * obj.Rx(qp(3)); %roll
            H34 = obj.Tx(obj.drw.lRoll)*obj.Ry(kns.ry34)*obj.Rz(kns.rz34) * obj.Ty(r*(1-cos(al))) * obj.Tx(r*sin(al)) * obj.Rz(al); %ry34 for out-of-articulation plane
            H45 = obj.Ry(kns.ry45) * obj.Tx( qp(5) ); %translation along x4 to the target

            varargout = {H01*H12*H23*H34*H45}; %H05
            if nargout == 5;
                varargout = [{H01}, {H01*H12}, {H01*H12*H23}, {H01*H12*H23*H34}, {H01*H12*H23*H34*H45}]; %return H01,H02,H03,H04,H05
            end
        end %FK
    end %public methods
     
    methods (Static)
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
                error('Matlab:inter2D_kn11A','kna does not have 11 elements');
            end
        end %paramArray2Struct
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
                error('Matlab:inter2D_kn11A','kns is not a struct');
            end
        end %paramStruct2Array
        
    end %statics
    
    
end