%% function xupdate= uhlprocsim(x,t,dt,method,Qprocch)
%
% Simulates xdot= f(x) using Euler or midpoint methods
% x= state, t= time, dt= time increment
% method= 'Euler' or 'midpoint'
% Qprocch= Cholesky factor of proc noise; if not present,
%    uses default value defined internally to uhlproc
%
% Noise free proc equation:  xdot= uhlproc(x,t,0)
%
% v= disturbance 
%
% Euler method:
%   xupdate= x+dt*ulhproc(x,t,0)+dt*v
% Midpoint method:
%   xupdate= x+dt*f(t+dt/2,x+dt/2*f(x,t,0))+dt*v;
%
function xupdate= uhlprocsim(x,t,dt,method,Qprocch)

Qdefault= 0;
if nargin<4
    Qdefault= 1;
    UseMidpoint= 0;
    UseEuler= 1;   % default
    % Will use default Qprocch inside uhlproc function
elseif ischar(method(1))   % this is method
    if strcmp(lower(method(1)),'e') % use Euler
        UseEuler= 1;
        UseMidpoint= 0;
    elseif strcmp(lower(method(1)),'m') % use midpoint
        UseMidpoint= 1;
        UseEuler= 0;
    else
        disp('Error: specify Euler or midpoint');
        return
    end
    if nargin<5
        Qdefault= 1;
    end
elseif nargin<5   % 4th value is Q, no method specified
    Qprocch= method;
    UseEuler=1; % default
    UseMidpoint= 0;
else % 4th value is Qprocch, 5th value is method
    tmp= method;
    method= Qprocch;
    Qprocch= tmp;   % swap method and Qprocch
    if strcmp(lower(method(1)),'e') % use Euler
        UseEuler= 1;
        UseMidpoint= 0;
    elseif strcmp(lower(method(1)),'m') % use midpoint
        UseMidpoint= 1;
        UseEuler= 0;
    else
        disp('Error: specify Euler or midpoint');
        return
    end
end

if 0
   UseEuler
   UseMidpoint
   Qdefault
end

if UseEuler
    if Qdefault
        xupdate= x+dt*uhlproc(x,t);
    else
        xupdate= x+dt*uhlproc(x,t,Qprocch);
    end
elseif UseMidpoint
    if Qdefault
        xupdate= x+dt*uhlproc(x+dt/2*uhlproc(x,t,0),t+dt/2);
    else
        xupdate= x+dt*uhlproc(x+dt/2*uhlproc(x,t,0),t+dt/2,Qprocch);
    end
else
    disp('Error: only choices are Euler or midpoint');
    return
end

