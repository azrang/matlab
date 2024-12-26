%% function  y= uhlmeas(x,t,Qmeasch)
%
% Measurement Equation
% Computes measured y
%
function y= uhlmeasure(x,t,Qmeasch)

% Fixed parameters
% Radar Center
ctr= [6375;0];

% Default noise Cholesky factors
if nargin<3
    % Default Qmeasch
        Qmeasch= diag(sqrt([1,17e-3]));
end

y= [norm(x(1:2)-ctr);atan2(x(2)-ctr(2),x(1)-ctr(1))]+Qmeasch*randn(2,1);

end
