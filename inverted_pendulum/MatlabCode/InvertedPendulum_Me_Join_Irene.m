%% PENDOLO INVERSO DLQ CONTROL
addpath("Tesi\")
clear
close all
clc
rng(939);
plot_predictions = true;
load IP_0_100.mat
%% Dynamics ***************************************************************
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
% x0=[0;0;0;0];

alfa=@(t,x,u)  6*Lp^3*Rm*mp+8*Lb^2*Lp*Rm*mb+24*Lb^2*Lp*Rm*mp-6*Lp^3*Rm*mp*cos(x(2,:)).^2-18*Lb^2*Lp*Rm*mp*cos(x(2,:)).^2;
beta=@(t,x,u)  4*Lb^2*mb+12*Lb^2*mp-3*Lp^2*mp*cos(x(2,:)).^2+3*Lp^2*mp;
gamma=@(t,x,u) Lp^2*mp-( (9*Lb^2*Lp^2*mp^2*cos(x(2,:)).^2) / (4*Lb^2*mb+12*Lb^2*mp-3*Lp^2*mp*cos(x(2,:)).^2+3*Lp^2*mp) );

f3_1=@(t,x,u) ( 9*Lb*Lp^2*mp*Rm*( sin(x(2,:))-( sin(x(2,:)).^3 ) ) ) .* x(3,:).^2;  % QUI VA x(3) ALLA SECONDA, E NON ALLA TERZA
f3_2=@(t,x,u) ( 6*Lp^3*mp*Rm*sin(2*x(2,:)) )  .*x(3,:).*x(4,:);
f3_3=@(t,x,u) ( -24*b1*Lp*Rm-24*Kgi^2*km*kt*Lp*etag*etam )  *x(3,:);
f3_4=@(t,x,u) ( 12*Lb*Lp^2*mp*Rm*sin(x(2,:)) ) .* x(4,:).^2;
f3_5=@(t,x,u) ( 36*b2*Lb*Rm*cos(x(2,:)) ) .* x(4,:);
f3_6=@(t,x,u) ( 24*etag*etam*Kgi*kt*Lp*u + 9*g*Lb*Lp*mp*Rm*sin(2*x(2,:)) )  ;

f3=@(t,x,u) (f3_1(t,x,u) + f3_2(t,x,u) + f3_3(t,x,u)+ f3_4(t,x,u)+ f3_5(t,x,u)+ f3_6(t,x,u)).*(1./alfa(t,x,u));


f4_1=@(t,x,u) 3*Lp^2*mp*cos(x(2,:)).*sin(x(2,:)).*beta(t,x,u) .*x(3,:).^2;  % QUI VA x(3) ALLA SECONDA, E NON ALLA TERZA
f4_2=@(t,x,u) -36*Lb*Lp^3*mp^2*cos(x(2,:)).^2.*sin(x(2,:)) .*x(3,:).*x(4,:);
f4_3=@(t,x,u) -82*Lb*Lp*mp*cos(x(2,:)) .*(b1 + (etag*etam*Kgi^2*km*kt)/Rm  ) .*x(3,:);
f4_4=@(t,x,u) -36*Lb^2*Lp^2*mp^2*sin(x(2,:)).*cos(x(2,:)) .* x(4,:).^2;
f4_5=@(t,x,u) -12*b2.*beta(t,x,u).*x(4,:);
f4_6=@(t,x,u) ( ( ( ( 6*etag*etam*Kgi*kt*Lb*Lp*mp*cos(x(2,:)).*u ) ./ ( Rm*(4*Lb^2*mb+12*Lb^2*mp-3*Lp^2*mp*cos(x(2,:)).^2+3*Lp^2*mp) ) ) )+0.5*g*Lp*mp*sin(x(2,:)) )*12.*beta(t,x,u);

f4=@(t,x,u) (f4_1(t,x,u) + f4_2(t,x,u) + f4_3(t,x,u)+ f4_4(t,x,u)+ f4_5(t,x,u)+ f4_6(t,x,u)).* (1./(4*beta(t,x,u).*gamma(t,x,u)));
f_u = @(t,x,u)([ x(3,:) ; x(4,:) ; f3(t,x,u) ; f4(t,x,u)] );

n = 4;
m = 1; % number of control inputs.
%% Discretization *********************************************************

deltaT = 0.01;
%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

%% LQ control *************************************************************
% Local linearization predictor at 0
x = sym('x',[n;1]); u = sym('u',[m;1]);
Ac_0 = double(subs(jacobian(f_u(0,x,u),x),[x;u],[zeros(size(x)); zeros(size(u))]));
Bc_0 = double(subs(jacobian(f_u(0,x,u),u),[x;u],[zeros(size(x)); zeros(size(u))]));
% discretize
[Ad_0, Bd_0] = c2d(Ac_0, Bc_0, deltaT);

Q = 100*eye(n);
R = eye(m);
[K_DLQ,S,eigen] = dlqr(Ad_0, Bd_0, Q, R);


%% Collect data ***********************************************************

tic
disp('Starting data collection')
% Numero di secondi per traiettoria = 20 come per VDP
%Nsim/(Ntraj*deltaT)
Nsim = 200;
Ntraj = 1000;

% Random initial conditions
Xcurrent = 0.1*randn(n,Ntraj);

X = []; Y = []; U = [];

for i = 1:Nsim
    u = -K_DLQ*Xcurrent + 0.1*rand(1)-0.05;
    Xnext = f_ud(0,Xcurrent,u);
    X = [X Xcurrent];
    Y = [Y Xnext];
    U = [U u];

    Xcurrent = Xnext; %aggiorna lo stato corrente del sistema per la prossima iterazione del ciclo.
end
fprintf('Data collection DONE, time = %1.2f s \n', toc);


%% ************************** Basis functions *****************************
basisFunction = 'rbf';
rbf_type = 'polyharmonic';
simOut1 = cell(1, 100);
simOut2 = cell(1, 100);
simOut3 = cell(1, 100);
for Nrbf = 1:100
    rng(10);
    cond=10000;
    %while(cond>5000) %ricalcola i centri "cent" fin quando non ho una KliftS buona
        % RBF centers
        cent = rand(n,Nrbf)*6.28 - 3.14;
        % Lifting mapping - RBFs + the state itself
        liftFun = @(xx)( [xx;rbf(xx,cent,rbf_type)-rbf(zeros(n,1),cent,rbf_type)] ); % funzione che prende in input lo stato e lo solleva tramite le rbf
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
        M = WVt * pinv(VVt); % Matrix [A B; C 0]
        Alift = M(1:Nlift,1:Nlift);
        Blift = M(1:Nlift,Nlift+1:end);
        Clift = M(Nlift+1:end,1:Nlift);

        fprintf('Regression done, time = %1.2f s \n', toc);
        disp(Nrbf)

        %% DLQ CONTROL KOOPMAN

        xbar=[0;pi/10;0;0];
        Qlift =100*eye(Nrbf+4);
        Rlift= eye(m);
        QliftS = 100*eye(Nrbf+4+1);
        RliftS= Nrbf;

        %[Klift, ~, Eigenvalues] = dlqr(Alift, Blift, Qlift, Rlift);
        flag1=1;
        flag2=1;
        flag3=1;
        try
            Klift=DLQControl(Alift, Blift, Qlift, Rlift,liftFun(xbar),deltaT);
        catch ME
            flag1=0;
            warning("Error using dlqr: Cannot compute a stabilizing LQR gain (DLQ)-->SALTO ALLA PROSSIMA ITERAZIONE");
            disp(Nrbf);
            break;
        end


    %% Simulink

    if(flag1==1)
        t_sim=20;
        try
            simOut1{Nrbf} = sim('IP_DLQ', 'SimulationMode', 'normal', 'SaveOutput', 'on');
        catch ME
            warning("Warning: Solver is encountering difficulty in simulating model");
            continue;
        end
    end
end
%% ********************** Analisi dei risultati *********************************

%array che contengono errori finali e tempi di assestamento
errori_finali = NaN(1, numel(simOut1));
tempi_assestamento=NaN(1,numel(simOut1));


%cerco il tempo di assestamento: appena esco dalla soglia vuol dire che
%il t precedente era il tempo di assestamento
for i = 1:numel(simOut1)
    %DLQ normale
    if ~isempty(simOut1{i})
        segnale = simOut1{i}.yout{1}.Values.Data(:,2);
        valore_finale = segnale(end);
        valore_iniziale = segnale(1);
        soglia = 0.05 * abs(valore_iniziale);
        tempo_assestamento=0;
        for t = numel(segnale):-1:1
            if abs(segnale(t) - valore_finale) >= soglia
                tempo_assestamento = simOut1{i}.yout{1}.Values.Time(t);
                tempi_assestamento(i)=tempo_assestamento;
                break;
            end
        end
    end
    
    % Estraggo gli errori finali da ciascuna simulazione
    if ~isempty(simOut1{i})
        errore=simOut1{i}.yout{2}.Values.Data(end);
        errori_finali(i) = abs(errore);
        fprintf('Errore del segnale con %3d lift function : %.10f\n', i, errore);
    end
end

%% PLOT DELLE SIMULAZIONI RIUSCITE (QUINDI SENZA INTEGRATORE)
figure()
subplot(2,1,1), plot(errori_finali, 'o', 'MarkerSize', 8), title("Errore finale", 'FontSize', 20), xlabel("n° funzioni sollevamento", 'FontSize', 16);
subplot(2,1,2), plot(tempi_assestamento, 'o', 'MarkerSize', 8), title("Tempo assestamento", 'FontSize', 20), xlabel("n° funzioni sollevamento", 'FontSize', 16);


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
tfinal=40;
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

