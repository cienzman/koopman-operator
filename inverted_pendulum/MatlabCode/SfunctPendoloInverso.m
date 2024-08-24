function [sys,x0,str,ts,simStateCompliance] = ...
    SfunctPendoloInverso(t, x, u, flag, x0,mb,mp,Lb,Lp,b1,b2,g,Rm,kt,km,Kgi,etag,etam)

switch flag

  %%%%%%%%%%%%%%%%%%
  % Initialization %
  %%%%%%%%%%%%%%%%%%
  case 0
    [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(x0);

  %%%%%%%%%%%%%%%
  % Derivatives %
  %%%%%%%%%%%%%%%
  case 1
    sys=mdlDerivatives(t,x,u,mb,mp,Lb,Lp,b1,b2,g,Rm,kt,km,Kgi,etag,etam);

  %%%%%%%%%%
  % Update %
  %%%%%%%%%%
  case 2
    sys=mdlUpdate(t,x,u);

  %%%%%%%%%%%
  % Outputs %
  %%%%%%%%%%%
  case 3
    sys=mdlOutputs(t,x,u); %in questo caso, non servono parametri aggiuntivi

  %%%%%%%%%%%%%%%%%%%%%%%
  % GetTimeOfNextVarHit %
  %%%%%%%%%%%%%%%%%%%%%%%
  case 4
    sys=mdlGetTimeOfNextVarHit(t,x,u);

  %%%%%%%%%%%%%
  % Terminate %
  %%%%%%%%%%%%%
  case 9
    sys=mdlTerminate(t,x,u);

  %%%%%%%%%%%%%%%%%%%%
  % Unexpected flags %
  %%%%%%%%%%%%%%%%%%%%
  otherwise
    DAStudio.error('Simulink:blocks:unhandledFlag', num2str(flag));

end

% end sfuntmpl

%
%=============================================================================
% mdlInitializeSizes
% Return the sizes, initial conditions, and sample times for the S-function.
%=============================================================================
%
function [sys,x0,str,ts,simStateCompliance]=mdlInitializeSizes(x0)

%
% call simsizes for a sizes structure, fill it in and convert it to a
% sizes array.
%
% Note that in this example, the values are hard coded.  This is not a
% recommended practice as the characteristics of the block are typically
% defined by the S-function parameters.
%
sizes = simsizes;

sizes.NumContStates  = 4; % DA SETTARE
sizes.NumDiscStates  = 0;
sizes.NumOutputs     = 4; % DA SETTARE
sizes.NumInputs      = 1; % DA SETTARE
sizes.DirFeedthrough = 0;
sizes.NumSampleTimes = 1;   % at least one sample time is needed %controlla questo

sys = simsizes(sizes);

%
% initialize the initial conditions
%
%x0  = [0 0]';


%
% str is always an empty matrix
%
str = [];

%
% initialize the array of sample times
%
ts  = [0 0];

% Specify the block simStateCompliance. The allowed values are:
%    'UnknownSimState', < The default setting; warn and assume DefaultSimState
%    'DefaultSimState', < Same sim state as a built-in block
%    'HasNoSimState',   < No sim state
%    'DisallowSimState' < Error out when saving or restoring the model sim state
simStateCompliance = 'UnknownSimState';

% end mdlInitializeSizes

%
%=============================================================================
% mdlDerivatives
% Return the derivatives for the continuous states.
%=============================================================================
%
function sys=mdlDerivatives(t,x,u,mb,mp,Lb,Lp,b1,b2,g,Rm,kt,km,Kgi,etag,etam)

% {
 alfa= 6*Lp^3*Rm*mp+8*Lb^2*Lp*Rm*mb+24*Lb^2*Lp*Rm*mp-6*Lp^3*Rm*mp*cos(x(2))^2-18*Lb^2*Lp*Rm*mp*cos(x(2))^2;
 beta= 4*Lb^2*mb+12*Lb^2*mp-3*Lp^2*mp*cos(x(2))^2+3*Lp^2*mp;
gamma= Lp^2*mp-( (9*Lb^2*Lp^2*mp^2*cos(x(2))^2) / (4*Lb^2*mb+12*Lb^2*mp-3*Lp^2*mp*cos(x(2))^2+3*Lp^2*mp) );

f3_1= ( 9*Lb*Lp^2*mp*Rm*( sin(x(2))-( sin(x(2))^3 ) ) ) * x(3)^3;
f3_2= ( 6*Lp^3*mp*Rm*sin(2*x(2)) )  *x(3)*x(4);
f3_3= ( -24*b1*Lp*Rm-24*Kgi^2*km*kt*Lp*etag*etam )  *x(3);
f3_4= ( 12*Lb*Lp^2*mp*Rm*sin(x(2)) ) * x(4)^2;
f3_5= ( 36*b2*Lb*Rm*cos(x(2)) ) * x(4);
f3_6= ( 24*etag*etam*Kgi*kt*Lp*u + 9*g*Lb*Lp*mp*Rm*sin(2*x(2)) )  ;

f3= (f3_1 + f3_2 + f3_3+ f3_4+ f3_5+ f3_6)*(1/alfa);


f4_1= 3*Lp^2*mp*cos(x(2))*sin(x(2))*beta *x(3)^3;
f4_2= -36*Lb*Lp^3*mp^2*cos(x(2))^2*sin(x(2)) *x(3)*x(4);
f4_3= -82*Lb*Lp*mp*cos(x(2)) *(b1 + (etag*etam*Kgi^2*km*kt)/Rm  ) *x(3);
f4_4= -36*Lb^2*Lp^2*mp^2*sin(x(2))*cos(x(2)) * x(4)^2;
f4_5= -12*b2*beta*x(4);
f4_6= ( ( ( ( 6*etag*etam*Kgi*kt*Lb*Lp*mp*cos(x(2))*u ) / ( Rm*(4*Lb^2*mb+12*Lb^2*mp-3*Lp^2*mp*cos(x(2))^2+3*Lp^2*mp) ) ) )+0.5*g*Lp*mp*sin(x(2)) )*12*beta; 

f4= (f4_1 + f4_2 + f4_3+ f4_4 + f4_5+ f4_6)* (1/(4*beta*gamma));
sys(1)=x(3);
sys(2)=x(4);
sys(3)=f3;
sys(4)=f4;
% }


%{
alfa=  6*Lp^3*Rm*mp+8*Lb^2*Lp*Rm*mb+24*Lb^2*Lp*Rm*mp-6*Lp^3*Rm*mp*cos(x(2)).^2-18*Lb^2*Lp*Rm*mp*cos(x(2)).^2;
beta=  4*Lb^2*mb+12*Lb^2*mp-3*Lp^2*mp*cos(x(2)).^2+3*Lp^2*mp;
gamma= Lp^2*mp-( (9*Lb^2*Lp^2*mp^2*cos(x(2)).^2) / (4*Lb^2*mb+12*Lb^2*mp-3*Lp^2*mp*cos(x(2)).^2+3*Lp^2*mp) );

f3_1= ( 9*Lb*Lp^2*mp*Rm*( sin(x(2))-( sin(x(2)).^3 ) ) ) .* x(3).^3;
f3_2= ( 6*Lp^3*mp*Rm*sin(2*x(2)) )  .*x(3).*x(4);
f3_3=( -24*b1*Lp*Rm-24*Kgi^2*km*kt*Lp*etag*etam )  *x(3);
f3_4= ( 12*Lb*Lp^2*mp*Rm*sin(x(2)) ) .* x(4).^2;
f3_5= ( 36*b2*Lb*Rm*cos(x(2)) ) .* x(4);
f3_6= ( 24*etag*etam*Kgi*kt*Lp*u + 9*g*Lb*Lp*mp*Rm*sin(2*x(2)) )  ;

f3=@(t,x,u) (f3_1 + f3_2 + f3_3+ f3_4+ f3_5+ f3_6).*(1./alfa);


f4_1= 3*Lp^2*mp*cos(x(2)).*sin(x(2)).*beta(t,x,u) .*x(3).^3;
f4_2= -36*Lb*Lp^3*mp^2*cos(x(2)).^2.*sin(x(2)) .*x(3).*x(4);
f4_3= -82*Lb*Lp*mp*cos(x(2)) .*(b1 + (etag*etam*Kgi^2*km*kt)/Rm  ) .*x(3);
f4_4= -36*Lb^2*Lp^2*mp^2*sin(x(2)).*cos(x(2)) .* x(4).^2;
f4_5= -12*b2.*beta(t,x,u).*x(4);
f4_6= ( ( ( ( 6*etag*etam*Kgi*kt*Lb*Lp*mp*cos(x(2)).*u ) ./ ( Rm*(4*Lb^2*mb+12*Lb^2*mp-3*Lp^2*mp*cos(x(2)).^2+3*Lp^2*mp) ) ) )+0.5*g*Lp*mp*sin(x(2)) )*12.*beta(t,x,u);

f4=(f4_1 + f4_2 + f4_3+ f4_4+ f4_5+ f4_6).* (1./(4*beta.*gamma));

sys(1)=x(3);
sys(2)=x(4);
sys(3)=f3;
sys(4)=f4;
%}
 %end mdlDerivatives

%
%=============================================================================
% mdlUpdate
% Handle discrete state updates, sample time hits, and major time step
% requirements.
%=============================================================================
%
function sys=mdlUpdate(t,x,u)

sys = [];

% end mdlUpdate

%
%=============================================================================
% mdlOutputs
% Return the block outputs.
%=============================================================================
%
function sys=mdlOutputs(t,x,u)


    sys(1) = x(1); % DEFINIRE y(t) %%y(t) = x1 = x(1)
    sys(2) = x(2);
    sys(3) = x(3);
    sys(4) = x(4);
% end mdlOutputs

%
%=============================================================================
% mdlGetTimeOfNextVarHit
% Return the time of the next hit for this block.  Note that the result is
% absolute time.  Note that this function is only used when you specify a
% variable discrete-time sample time [-2 0] in the sample time array in
% mdlInitializeSizes.
%=============================================================================
%
function sys=mdlGetTimeOfNextVarHit(t,x,u)
 
sampleTime = 1;    %  Example, set the next hit to be one second later.
sys = t + sampleTime;

% end mdlGetTimeOfNextVarHit

%
%=============================================================================
% mdlTerminate
% Perform any end of simulation tasks.
%=============================================================================
%
function sys=mdlTerminate(t,x,u)

sys = [];

% end mdlTerminate
