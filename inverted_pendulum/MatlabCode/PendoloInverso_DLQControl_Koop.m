%% PENDOLO INVERSO DLQ CONTROL
clear all
close all
addpath("Tesi\")
load Data_for_KoopmanDLQ_InvertedPendolum ;
load IP_0_60.mat
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
x0=[0;0;0;0];
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
simOut1 = cell(1, 100);
simOut2 = cell(1, 100);
simOut3 = cell(1, 100);
for Nrbf = 1:110
    rng(10);
    cond=10000;
    while(cond>5000) %ricalcola i centri "cent" fin quando non ho una KliftS buona
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

        %% DLQ CONTROL KOOPMAN

        xbar=[0;pi/100;0;0];
        Qlift =100*eye(Nrbf+4);
        Rlift= Nrbf;
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
        %title("DLQ CONTROL KOOPMAN SENZA INTEGRATORE")

        %% DLQ CONTROL KOOPMAN con integratore eulero avanti
        %{
        QliftS(Nrbf+2+1,Nrbf+2+1)=10*Nrbf;
        AliftS=[Alift zeros(size(Alift,1),1);-[1,zeros(1,size(Alift,1)-1)] 1];
        BliftS=[Blift;0];
        %[KliftS, ~, EigenvaluesS] = dlqr(AliftS, BliftS, QliftS, RliftS);
        KliftS_EuleroAvanti=DLQControl(AliftS, BliftS, QliftS, RliftS,[liftFun(xbar);0],deltaT);
        %}
        %% DLQ CONTROL KOOPMAN con integratore eulero indietro

        QliftS(Nrbf+4+1,Nrbf+4+1)=100*Nrbf; %stato relativo all'integratore 
        AliftS=[Alift zeros(size(Alift,1),1);-[1,zeros(1,size(Alift,1)-1)]*Alift 1];
        BliftS=[Blift;-[1,zeros(1,size(Alift,1)-1)]*Blift];
        %[KliftS, ~, EigenvaluesS] = dlqr(AliftS, BliftS, QliftS, RliftS);
        try
            KliftS_EuleroIndietro=DLQControl(AliftS, BliftS, QliftS, RliftS,[liftFun(xbar);0],deltaT);
        catch ME
            flag2=0;
            warning("Error using dlqr: Cannot compute a stabilizing LQR gain (eulero indietro)-->SALTO ALLA PROSSIMA ITERAZIONE");
            disp(Nrbf);
            break;
        end

        %title("integratore eulero indietro");


        %% DLQ CONTROL KOOPMAN con integratore (variazioni)

        QliftS(Nrbf+4+1,Nrbf+4+1)=100*Nrbf;
        AliftS=[Alift zeros(size(Alift,1),1);-[1,zeros(1,size(Alift,1)-1)]*Alift 1];
        BliftS=[Blift;-[1,zeros(1,size(Alift,1)-1)]*Blift];
        %[KliftS, ~, EigenvaluesS] = dlqr(AliftS, BliftS, QliftS, RliftS);
        try
            KliftS=DLQControl(AliftS, BliftS, QliftS, RliftS,[liftFun(xbar);0],deltaT);
        catch ME
            flag3=0;
            warning("Error using dlqr: Cannot compute a stabilizing LQR gain integratore variazioni-->SALTO ALLA PROSSIMA ITERAZIONE");
            disp(Nrbf);
            break;
        end
        %title("Variazioni");

        cond=max(abs(KliftS));
        disp(cond);


    end


    %% Simulink

    if(flag1==1)
        t_sim=30;
        try
            simOut1{Nrbf} = sim('IP_DLQ', 'SimulationMode', 'normal', 'SaveOutput', 'on');
        catch ME
            warning("Warning: Solver is encountering difficulty in simulating model");
            continue;
        end
    end
    if(flag2==1)
        t_sim=30;
        try
            simOut2{Nrbf} = sim('IP_DLQ_EI', 'SimulationMode', 'normal', 'SaveOutput', 'on');
        catch ME
            warning("Warning: Solver is encountering difficulty in simulating model");
            continue;
        end
    end
    if(flag3==1)
        t_sim=30;
        try
            simOut3{Nrbf} = sim('IP_DLQ_VAR', 'SimulationMode', 'normal', 'SaveOutput', 'on');
        catch ME
            warning("Warning: Solver is encountering difficulty in simulating model");
            continue;
        end
    end

end
%% ********************** Analisi dei risultati *********************************

%array che contengono errori finali e tempi di assestamento
errori_finali = NaN(1, numel(simOut1));
errori_finali_EI = NaN(1, numel(simOut2));
errori_finali_VAR = NaN(1, numel(simOut3));
tempi_assestamento=NaN(1,numel(simOut1));
tempi_assestamento_EI=NaN(1,numel(simOut2));
tempi_assestamento_VAR=NaN(1,numel(simOut3));

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
    %DLQ+integratoreEI
    if ~isempty(simOut2{i})
        segnale_EI = simOut2{i}.yout{1}.Values.Data(:,2);
        valore_finale_EI = segnale_EI(end);
        valore_iniziale_EI = segnale_EI(1);
        soglia_EI = 0.05 * abs(valore_iniziale_EI);
        tempo_assestamento_EI=0;
        for t = numel(segnale_EI):-1:1
            if abs(segnale_EI(t) - valore_finale_EI) >= soglia_EI
                tempo_assestamento_EI = simOut2{i}.yout{1}.Values.Time(t);
                tempi_assestamento_EI(i)=tempo_assestamento_EI;
                break;
            end
        end
    end
    %DLQ+integratoreVAR
    if ~isempty(simOut3{i})
        segnale_VAR = simOut3{i}.yout{1}.Values.Data(:,2);
        valore_finale_VAR = segnale_VAR(end);
        valore_iniziale_VAR = segnale_VAR(1);
        soglia_VAR = 0.05 * abs(valore_iniziale_VAR);
        tempo_assestamento_VAR=0;
        for t = numel(segnale_VAR):-1:1
            if abs(segnale_VAR(t) - valore_finale_VAR) >= soglia_VAR
                tempo_assestamento_VAR = simOut3{i}.yout{1}.Values.Time(t);
                tempi_assestamento_VAR(i)=tempo_assestamento_VAR;
                break;
            end
        end
    end
    %fprintf('Tempo di assestamento del segnale con %3d lift function : %.2f\n', i, tempo_assestamento);
    %fprintf('Tempo di assestamento del segnale con %3d lift function (DLQ+integratoreEI) : %.2f\n', i, tempo_assestamento_EI);
    %fprintf('Tempo di assestamento del segnale con %3d lift function (DLQ+integratoreVAR) : %.2f\n', i, tempo_assestamento_VAR);
    % Estraggo gli errori finali da ciascuna simulazione
    if ~isempty(simOut1{i})
        errore=simOut1{i}.yout{2}.Values.Data(end);
        errori_finali(i) = abs(errore);
        fprintf('Errore del segnale con %3d lift function : %.10f\n', i, errore);
    end
    if ~isempty(simOut2{i})
        errore_EI=simOut2{i}.yout{2}.Values.Data(end);
        errori_finali_EI(i) = abs(errore_EI);
        fprintf('Errore del segnale con %3d lift function (DLQ+integratoreEI) : %.10f\n', i, errore_EI);
    end
    if ~isempty(simOut3{i})
        errore_VAR=simOut3{i}.yout{2}.Values.Data(end);
        errori_finali_VAR(i) = abs(errore_VAR);
        fprintf('Errore del segnale con %3d lift function (DLQ+integratoreVAR) : %.10f\n', i, errore_VAR);
    end


    % Salvo l'ultimo valore dell'errore nella cella degli ultimi errori


end
%% PLOT DI TUTTO INSIEME
figure()
subplot(6,1,1),plot(errori_finali,'o'),title("errore finale"),xlabel("numero di lift function");
subplot(6,1,2),plot(errori_finali_EI,'o'),title("errore finale con integratore EI"),xlabel("numero di lift function");
subplot(6,1,3),plot(errori_finali_VAR,'o'),title("errore finale con integratore VAR"),xlabel("numero di lift function");
subplot(6,1,4),plot(tempi_assestamento,'o'),title("tempo assestamento"),xlabel("numero di lift function");
subplot(6,1,5),plot(tempi_assestamento_EI,'o'),title("tempo assestamento con integratore EI"),xlabel("numero di lift function");
subplot(6,1,6),plot(tempi_assestamento_VAR,'o'),title("tempo assestamento con integratore VAR"),xlabel("numero di lift function");
%% PLOT DELLE SIMULAZIONI RIUSCITE (QUINDI SENZA INTEGRATORE)
figure()
subplot(2,1,1),plot(errori_finali,'o'),title("errore finale"),xlabel("numero di lift function");
subplot(2,1,2),plot(tempi_assestamento,'o'),title("tempo assestamento"),xlabel("numero di lift function");
%% PLOT INTEGRATORE EI
figure()
subplot(2,1,1),plot(errori_finali_EI,'o'),title("errore finale EI"),xlabel("numero di lift function");
subplot(2,1,2),plot(tempi_assestamento_EI,'o'),title("tempo assestamento EI"),xlabel("numero di lift function");
%% PLOT INTEGRATORE VARIAZIONI
figure()
subplot(2,1,1),plot(errori_finali_VAR,'o'),title("errore finale EI"),xlabel("numero di lift function");
subplot(2,1,2),plot(tempi_assestamento_VAR,'o'),title("tempo assestamento EI"),xlabel("numero di lift function");
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

