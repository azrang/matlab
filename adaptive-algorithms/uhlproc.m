%% function  xdot= uhlproc(x,t,Qprocch)
%
% Process Equation
% Computes xdot value
%
%
% dt= time increment
% x= [rx,ry,vx,vy,logbeta]:  (rx,ry)= 2D position, (vx,vy)= 2D velocity
% 
% True drag coefficient beta unknown.
% Modeled as beta= beta0*exp(logbeta)  with nominal (known) beta0
%     Parametrized by logbeta because beta0>0 must be sustained
%     for stability.
%
% Qprocch= Cholesky factor of process noise
%
function xdot= uhlproc(x,t,Qprocch)

N= size(x,1); % # of states
UseEuler= 1;   % if 0, use midpoint

% Fixed parameters
% ballistic coefficient
beta0= 0.597983;
% Other stuff
H0= 13.406;
Gm0= 3.986e5;
R0= 6374;

% Default noise Cholesky factors
if nargin<3
    % Default Qprocch
        %Qprocch= diag(sqrt([0,0,2.4064e-5,2.4064e-5,0.005]));
        Qprocch= diag(sqrt([0,0,2.4064e-5,2.4064e-5,0]));
        % No direct disturbance on position, only velocity
        % No disturbance on logbeta, so beta is in theory
        %   (unknown) constant
end

% Drag Coefficient
beta= beta0*exp(x(5));

% Computing quantities
R= norm(x(1:2)); % distance from center of earth
V= norm(x(3:4)); % speed
D= -beta*V*exp((R0-R)/H0);  % drag
G= -Gm0/R^3;  % gravity force related term
if 0
    [D,G]
end

v= Qprocch*randn(N,1);
xdot= [x(3);x(4);D*x(3)+G*x(1);D*x(4)+G*x(2);0]+v;

end
