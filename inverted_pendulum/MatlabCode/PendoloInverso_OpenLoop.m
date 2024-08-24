clear all;
close all;
load Data_for_KoopmanDLQ_InvertedPendolum;
rng(102837);
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

n = 4;
m = 1; % number of control inputs.
%% ************************** Discretization ******************************

%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );

f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   ); %f_ud rappresenta la dinamica del sistema discretizzata.

%% ************************** Collect data ********************************

tic
disp('Starting data collection')
% Numero di secondi per traiettoria = 20 come per VDP
%Nsim/(Ntraj*deltaT)
Nsim = 10;
Ntraj = 500;
Nsim=200;
Ntraj=1000;

% Random initial conditions
Xcurrent = [];
for i = 1:Ntraj
    colonna = randi([1, size(InitialConditions,2)]);
    Xrandom=InitialConditions(:,colonna);
    Xcurrent = [Xcurrent Xrandom];
end

X = []; Y = []; U = [];

for i = 1:Nsim
    Xnext = f_ud(0,Xcurrent,-K_DLQ*Xcurrent);
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U -K_DLQ*Xcurrent];
    Xcurrent = Xnext; %aggiorna lo stato corrente del sistema per la prossima iterazione del ciclo.
end
fprintf('Data collection DONE, time = %1.2f s \n', toc);
%% ************************** Basis functions *****************************
basisFunction = 'rbf';
rbf_type = 'thinplate';
for Nrbf = 100:100
    % RBF centers
    cent = rand(n,Nrbf)*2 - 1;
    % Lifting mapping - RBFs + the state itself
    liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)] ); % funzione che prende in input lo stato e lo solleva tramite le rbf
    Nlift = Nrbf + n;

    %% ******************************* Lift ***********************************

    disp('Starting LIFTING')
    tic
    Xlift = liftFun(X);
    Ylift = liftFun(Y);

    fprintf('Lifting DONE, time = %1.2f s \n', toc);
    %% ********************** Build predictor *********************************

    disp('Starting REGRESSION')
    tic
    W = [Ylift ; X];
    V = [Xlift; U];
    VVt = V*V';
    WVt = W*V';
    M = WVt * pinv(VVt); % Matrix [A B; C 0] rappresenta il predittore lineare del sistema.
    Alift = M(1:Nlift,1:Nlift);
    Blift = M(1:Nlift,Nlift+1:end);
    Clift = M(Nlift+1:end,1:Nlift);

    fprintf('Regression done, time = %1.2f s \n', toc);
    %% *********************** Predictor comparison ***************************

    Tmax = 0.1; %tempo di simulazione
    Nsim = Tmax/deltaT;
    RMSE=0;
    randomly_initial_conditions=1;
    %for jbis=1:randomly_initial_conditions %average prediction RMSE over 100 randomly sampled initial conditions
    u_dt = @(i)((-1).^(round(i/30))); % control signal % funzione che restituisce un segnale di controllo che cambia di segno ogni 30 iterazioni
    % Initial condition
    x0 = [0;pi/100;0;0];
    x_true = x0;
    % Lifted initial condition
    xlift = liftFun(x0);

    xbar= [0;0;0;0];
    ubar=0;
    [Ac_xbar, Bc_xbar, c_xbar,Ad_bar,Bd_bar,cd_bar]= linearizza(xbar,ubar,f_u,deltaT);
    X_loc_x0 = x0;

    xbar1=[0;0;0;0];  
    ubar1=0;
    [Ac_xbar1, Bc_xbar1, c_xbar1,Ad_bar1,Bd_bar1,cd_bar1]= linearizza(xbar1,ubar1,f_u,deltaT);
    X_loc_0 = x0;


    % Simulate:
    for i = 0:Nsim-1
        % Koopman predictor
        xlift = [xlift, Alift*xlift(:,end) + Blift*u_dt(i)]; % Lifted dynamics

        % True dynamics
        x_true = [x_true, f_ud(0,x_true(:,end),u_dt(i)) ];

        % Local linearization predictor at x0
        X_loc_x0 = [X_loc_x0, Ad_bar*(X_loc_x0(:,end)-xbar) + Bd_bar*(u_dt(i)-ubar) + cd_bar+xbar];

        % Local linearization predictor at 0
        X_loc_0 = [X_loc_0, Ad_bar1*(X_loc_0(:,end)-xbar1) + Bd_bar1*(u_dt(i)-ubar1) + cd_bar1+xbar1];
    end
    x_koop = Clift * xlift; % Koopman predictions
    %% **************** root mean squared errors (RMSE) *******************
    NUM=0;
    DEN=0;

    for j=1:size(x_koop,2)
        NUM= NUM + norm(x_koop(:,j)-x_true(:,j) )^2   ;
        DEN= DEN + norm(x_true(:,j) )^2  ;
        RMSE=RMSE + 100* ( sqrt(NUM)/sqrt(DEN) );
    end
    %end %end average prediction RMSE over 100 randomly sampled initial
    %conditions SBAGLIATO PROBABILMENTE
    RMSE_AV(Nrbf)= RMSE/randomly_initial_conditions;
end %end Nrbf da 1 a 20 riga 85

%% ****************************  Plots  ***********************************

lw = 4;

figure
plot([0:Nsim-1]*deltaT,u_dt(0:Nsim-1),'linewidth',lw); hold on
title('Control input $u$', 'interpreter','latex'); xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)

figure %x1=teta1 angolo tra braccio e asse x
plot([0:Nsim]*deltaT,x_true(1,:),'linewidth',lw); hold on
plot([0:Nsim]*deltaT,x_koop(1,:), '--r','linewidth',lw)
plot([0:Nsim]*deltaT,X_loc_x0(1,:), '--g','linewidth',lw-1)
plot([0:Nsim]*deltaT,X_loc_0(1,:), '--k','linewidth',lw-1)
%axis([0 Tmax min(x_koop(1,:))-0.1 max(x_koop(1,:))+0.1])
title('Predictor comparison - $x_1$','interpreter','latex'); xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)
LEG = legend('True','Koopman','Local at $x_0$','Local at 0','location','southwest');
set(LEG,'interpreter','latex')

figure %x2=teta2 angolo tra asta pendolo e asse z
plot([0:Nsim]*deltaT,x_true(2,:),'linewidth',lw); hold on
plot([0:Nsim]*deltaT,x_koop(2,:), '--r','linewidth',lw)
plot([0:Nsim]*deltaT,X_loc_x0(2,:), '--g','linewidth',lw-1)
plot([0:Nsim]*deltaT,X_loc_0(2,:), '--k','linewidth',lw-1)
%axis([0 Tmax min(x_koop(2,:))-0.15 max(x_koop(2,:))+0.15])
title('Predictor comparison - $x_2$','interpreter','latex'); xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)
LEG = legend('True','Koopman','Local at $x_0$','Local at 0','location','southwest');
set(LEG,'interpreter','latex')

figure %x3=teta3 variazione (velocità) angolo tra braccio e asse x
plot([0:Nsim]*deltaT,x_true(3,:),'linewidth',lw); hold on
plot([0:Nsim]*deltaT,x_koop(3,:), '--r','linewidth',lw)
plot([0:Nsim]*deltaT,X_loc_x0(3,:), '--g','linewidth',lw-1)
plot([0:Nsim]*deltaT,X_loc_0(3,:), '--k','linewidth',lw-1)
%axis([0 Tmax min(x_koop(3,:))-0.15 max(x_koop(3,:))+0.15])
title('Predictor comparison - $x_3$','interpreter','latex'); xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)
LEG = legend('True','Koopman','Local at $x_0$','Local at 0','location','southwest');
set(LEG,'interpreter','latex')

figure %x4=teta4 variazione (velocità) angolo tra asta pendolo e asse z
plot([0:Nsim]*deltaT,x_true(4,:),'linewidth',lw); hold on
plot([0:Nsim]*deltaT,x_koop(4,:), '--r','linewidth',lw)
plot([0:Nsim]*deltaT,X_loc_x0(4,:), '--g','linewidth',lw-1)
plot([0:Nsim]*deltaT,X_loc_0(4,:), '--k','linewidth',lw-1)
%axis([0 Tmax min(x_koop(4,:))-0.15 max(x_koop(4,:))+0.15])
title('Predictor comparison - $x_4$','interpreter','latex'); xlabel('Time [s]','interpreter','latex');
set(gca,'fontsize',20)
LEG = legend('True','Koopman','Local at $x_0$','Local at 0','location','southwest');
set(LEG,'interpreter','latex')

%{
figure()
plot(X(:,1:1000)');
plot(X(1,1:1000),X(2,1:1000));
%}

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