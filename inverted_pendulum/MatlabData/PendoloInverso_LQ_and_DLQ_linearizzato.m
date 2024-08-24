%% PENDOLO INVERSO LQ CONTROL
clear all
close all
addpath("Tesi\")
rng(2141444)
%% *************************** Dynamics ***********************************
mb=0.257;
mp=0.127;
Lb=0.216;
Lp=0.337;
b1=0.009;
b2=5.2*10^(-4);
g=9.81;
Rm=2.6;
kt=0.00767;
km=0.00767;
Kgi=14;
etag=0.90;
etam=0.69;
%%
xOrigine=[0;0;0;0];
%%

alfa=@(t,x,u)  6*Lp^3*Rm*mp+8*Lb^2*Lp*Rm*mb+24*Lb^2*Lp*Rm*mp-6*Lp^3*Rm*mp*cos(x(2,:)).^2-18*Lb^2*Lp*Rm*mp*cos(x(2,:)).^2;
beta=@(t,x,u)  4*Lb^2*mb+12*Lb^2*mp-3*Lp^2*mp*cos(x(2,:)).^2+3*Lp^2*mp;
gamma=@(t,x,u) Lp^2*mp-( (9*Lb^2*Lp^2*mp^2*cos(x(2,:)).^2) / (4*Lb^2*mb+12*Lb^2*mp-3*Lp^2*mp*cos(x(2,:)).^2+3*Lp^2*mp) );

f3_1=@(t,x,u) ( 9*Lb*Lp^2*mp*Rm*( sin(x(2,:))-( sin(x(2,:)).^3 ) ) ) .* x(3,:).^2;
f3_2=@(t,x,u) ( 6*Lp^3*mp*Rm*sin(2*x(2,:)) )  .*x(3,:).*x(4,:);
f3_3=@(t,x,u) ( -24*b1*Lp*Rm-24*Kgi^2*km*kt*Lp*etag*etam )  *x(3,:);
f3_4=@(t,x,u) ( 12*Lb*Lp^2*mp*Rm*sin(x(2,:)) ) .* x(4,:).^2;
f3_5=@(t,x,u) ( 36*b2*Lb*Rm*cos(x(2,:)) ) .* x(4,:);
f3_6=@(t,x,u) ( 24*etag*etam*Kgi*kt*Lp*u + 9*g*Lb*Lp*mp*Rm*sin(2*x(2,:)) )  ;

f3=@(t,x,u) (f3_1(t,x,u) + f3_2(t,x,u) + f3_3(t,x,u)+ f3_4(t,x,u)+ f3_5(t,x,u)+ f3_6(t,x,u)).*(1./alfa(t,x,u));


f4_1=@(t,x,u) 3*Lp^2*mp*cos(x(2,:)).*sin(x(2,:)).*beta(t,x,u) .*x(3,:).^2;
f4_2=@(t,x,u) -36*Lb*Lp^3*mp^2*cos(x(2,:)).^2.*sin(x(2,:)) .*x(3,:).*x(4,:);
f4_3=@(t,x,u) -82*Lb*Lp*mp*cos(x(2,:)) .*(b1 + (etag*etam*Kgi^2*km*kt)/Rm  ) .*x(3,:);
f4_4=@(t,x,u) -36*Lb^2*Lp^2*mp^2*sin(x(2)).*cos(x(2,:)) .* x(4,:).^2;
f4_5=@(t,x,u) -12*b2.*beta(t,x,u).*x(4,:);
f4_6=@(t,x,u) ( ( ( ( 6*etag*etam*Kgi*kt*Lb*Lp*mp*cos(x(2,:)).*u ) ./ ( Rm*(4*Lb^2*mb+12*Lb^2*mp-3*Lp^2*mp*cos(x(2,:)).^2+3*Lp^2*mp) ) ) )+0.5*g*Lp*mp*sin(x(2,:)) )*12.*beta(t,x,u);

f4=@(t,x,u) (f4_1(t,x,u) + f4_2(t,x,u) + f4_3(t,x,u)+ f4_4(t,x,u)+ f4_5(t,x,u)+ f4_6(t,x,u)).* (1./(4*beta(t,x,u).*gamma(t,x,u)));
f_u = @(t,x,u)([ x(3,:) ; x(4,:) ; f3(t,x,u) ; f4(t,x,u)] );
max_x1=@(t,x,u) max(x(3,:));
max_x2=@(t,x,u) max(x(4,:));
max_x3=@(t,x,u) max(f3(t,x,u));
max_x4=@(t,x,u) max(f4(t,x,u));

n = 4;
m = 1; % number of control inputs.
%% ************************** Discretization ******************************

deltaT = 5/1000;
%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );

f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   ); %f_ud rappresenta la dinamica del sistema discretizzata.

%% Linearizzazione standard e LQ/DLQ CONTROL

% Definizione dei punti di equilibrio e dei relativi controlli
xbar= [0;0;0;0]; ubar=0;
xbar1=[0;pi/20;0;0]; ubar1=0;
xbar2=[0;pi/10;0;0]; ubar2=0;
xbar3=[pi/10;pi/20;0;0]; ubar3=0;
xbar4=[pi/10;pi/10;0;0]; ubar4=0;
xbar5=[pi;pi/10;0;0]; ubar5=0;
xbar6=[pi;pi/8;0;0]; ubar6=0;
xbar7=[pi;pi/5;0;0]; ubar7=0;
%xbar8=[pi;pi/4;0;0]; ubar8=0;
%xbar9=[pi;pi/3;0;0]; ubar9=0;
points = {xbar1,xbar2,xbar3,xbar4,xbar5,xbar6,xbar7};
controls = {ubar,ubar1,ubar2,ubar3,ubar4,ubar5,ubar6,ubar7};
num_points = numel(points);
% MATRIX Q AND R
Q = eye(4); %2 stati--> Q=2x2
R = eye(1); %1 ingresso--> R scalare

% Inizializzazione delle matrici per i controlli


% Calcolo dei controlli per ogni punto di equilibrio
simOut = cell(1, num_points);
    % Linearizzazione attorno L'origine
    [Ac, Bc, c, Ad, Bd, cd] = linearizza(xOrigine, controls{1}, f_u, deltaT);
    [Ad,Bd]=c2d(Ac,Bc,deltaT)

    % Calcolo dei controlli LQ
    %B_completa = [Bc, c];
    K_LQ = LQControl(Ac, Bc, Q, R, points{1});

    % Calcolo dei controlli DLQ
    K_DLQ = DLQControl(Ad, Bd, Q, R, points{1}, deltaT);
    InitialConditions=[];
for i = 1:num_points
    t_sim=6;
    simOut{i} = sim('PendoloInverso_LQControl_simulation', 'SimulationMode', 'normal', 'SaveOutput', 'on');
    InitialConditions=[InitialConditions simOut{i}.X'];
end
%%
save Data_for_KoopmanDLQ_InvertedPendolum InitialConditions K_DLQ deltaT

















%% ****************************  Function linearizza  ***********************************
function [Ac_xbar, Bc_xbar, c_xbar, Ad_xbar, Bd_xbar,cd_xbar]=linearizza(xbar,ubar,f_u,deltaT)

% Local linearization predictor at xbar
x = sym('x',[4;1]);
u = sym('u',[1;1]);
Ac_xbar = double(subs(jacobian(f_u(0,x,u),x),[x;u],[xbar;ubar]));
Bc_xbar = double(subs(jacobian(f_u(0,x,u),u),[x;u],[xbar;ubar]));

%c_bar = double(subs(f_u(0,x,u),[x;u],[xbar;ubar])) - Ac_x0*xbar - %Bc_x0*ubar; %f(xbar,ubar)-Ac_xbar-Bc_ubar
c_xbar= double(subs(f_u(0,x,u),[x;u],[xbar;ubar]));
ABc = expm([Ac_xbar [Bc_xbar c_xbar] ; zeros(2,6)]*deltaT); % discretize
Ad_xbar = ABc(1:4,1:4); Bd_xbar = ABc(1:4,3); cd_xbar = ABc(1:4,4);
%Ad_x0= expm((Ac_x0)*deltaT); ok coerente con la teoria
%Bd_x0 ? Non chiarissimo
end


%% FUNCTION LQControl
function [K]=LQControl(A,B,Q,R,xbar)
% %{
[K, ~, Eigenvalues] = lqr(A, B, Q, R);


Acl= A-(B*K);
Bcl=zeros(size(B));
Ccl=[0,1,0,0];
Dcl=0;
syscl=ss(Acl,Bcl,Ccl,Dcl);

syscl_NOTCTRL=ss(A,B,Ccl,0);
tfinal=100;
figure();
initial(syscl, xbar, tfinal)
hold on;
initial(syscl_NOTCTRL, xbar, tfinal)
grid on
legend('Con controllo', 'Senza controllo');
% %}

end

%% FUNCTION DLQControl
function [K]=DLQControl(A,B,Q,R,xbar,deltaT);
% %{
[K, ~, Eigenvalues] = dlqr(A, B, Q, R);
%{
Acl= A-(B*K);
Bcl=zeros(size(B));
Ccl=[1,zeros(1,size(Acl,1)-1)];
Dcl=0;
syscl=ss(Acl,Bcl,Ccl,Dcl,deltaT);

syscl_NOTCTRL=ss(A,B,[1,zeros(1,size(Acl,1)-1)],0,deltaT);
tfinal=40;
figure();
initial(syscl, xbar, tfinal)
hold on;
initial(syscl_NOTCTRL, xbar, tfinal)
grid on
legend('Con controllo', 'Senza controllo');
disp(max(abs(K)));
% %}
%}
end

