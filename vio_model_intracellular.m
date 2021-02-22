function vio_model_intracellular
close all
clc

Tmax = 3600*16; % sec

Ha = @(x,K,n) x.^n./(x.^n+K^n);
Hr = @(x,K,n) K^n./(x.^n+K^n);


%% parameters of balanced nutrient growth
nu_in = 1e5;
P = 6e5; % total intracellular gene concentration (nM) ~600 uM
M_r = 100;
M_v = 0;
M_p = 6200;
tau = 6.4e3; % [tRNA] 5200-6400
L_p = 300; % average protein length
L_r = 300; % ribosome length
L_v = 300; % violacein length

K_tau = 1e5; % dissociation constant of tRNA - a.a.
K_M = 100; % dissociation constant of ribosome - mRNA
K_v = 100;

k_on = 1; % ribosome and mRNA binding rate
k_init = 0.1; % ribosome-mRNA and tRNA binding rate
k_rec = 0.01; % stRNA rescue rate

zeta_r = 0; % chemical protein degradation rate
zeta_v = 0; % chemical vioE degradation rate
k_trans = 20; % translation rate of protein (a.a)
nu_vio = 1/60*1;%*1,8,16; % vioA catalysis rate
k_deg = 0.1;
kappa = 3e9;

M_vs = linspace(0,800,6);
Ms = M_r+M_p+M_vs;
ss = nan(6,length(M_vs));
for i=1:5
    figure('Position',[0 0 1200 600])
end
for i = 1:length(M_vs)
    M_v = M_vs(i);
    M = Ms(i); % 2000-8000 copies per cell
    L = L_p*M_p/M+L_r*M_r/M+L_v*M_v/M;
    init = zeros(5,1);
    init(1) = 4e4; % a, nM
    init(2) = 5e3; % Pr, nM
    init(4) = 1e4; %rho
    init(5) = 0; 
    [T,Y] = ode15s(@ode_fun, [0 Tmax],init);
    T = T/3600;
    ya = Y(:,1);
    ypr = Y(:,2);
    ypv = Y(:,3);
    yrho = Y(:,4);
    yvio = Y(:,5);
    yf1 = k_init./(K_tau+ya).*M*tau.*(L_p/k_trans*ya+K_tau/k_rec);
    gammy = tau*ya./(K_tau+ya)*M./(K_M+M+k_init/k_on*tau+yf1);
    sigmy = k_init*ypr.*gammy;
    muy = sigmy*(M_p+M_r)/M*log(2)/P;

    figure(1)
    subplot(2,3,i)
    plot(T,ya./(K_tau+ya));
    hold all
    axis([0 Tmax/3600 0 1])
    title('L-tryptophan')

    figure(2)
    subplot(2,3,i)
    plot(T,ypr);
    hold all;
    axis([0 Tmax/3600 0 1.5e4])
    title('Pr')

    figure(3)
    subplot(2,3,i)
    plot(T,ypv);
    axis([0 Tmax/3600 0 1e5])
    hold all;
    title('Pv')

    figure(4)
    subplot(2,3,i)
    plot(T,yrho);
    hold all;
    axis([0 Tmax/3600 0 3e9])
    title('\rho')
    
    figure(5)
    subplot(2,3,i)
    plot(T,yvio);
    hold all;
    axis([0 Tmax/3600 0 2.5e8])
    title('vio')
    
    ss(1:5,i) = Y(end,:);
    ss(6,i) = muy(end);
end

%% steady state a
figure;
x = M_vs./Ms*100;
plot(x,ss(1,:)./(ss(1,:)+K_tau),'o','linewidth',2);
p = polyfit(x,ss(1,:)./(ss(1,:)+K_tau),1);
hold on
plot(x,polyval(p,x),'r-','linewidth',2)
xlabel('M_v (%)');
ylabel('a/(K_{\tau}+a)');
axis([0 12 0.2 1.0])
set(gca,'fontsize',12);
title('Mv vs a')

%% steady state P_r
figure;
plot(x,ss(2,:),'o','linewidth',2);
p = polyfit(x,ss(2,:),1);
hold on
plot(x,polyval(p,x),'r-','linewidth',2)
xlabel('M_v (%)');
ylabel('P_r (nM)');
axis([0 12 1.2e4 1.38e4])
set(gca,'fontsize',12);
title('Mv vs Pr');

%% steady state gamma
figure;
aa = ss(1,:);
ss_gamma = Ms.*tau*k_init*L/k_trans.*Ha(aa,K_tau,1)./(K_M+k_init/k_on*tau+Ms+Ms.*tau*k_init*L/k_trans.*Ha(aa,K_tau,1)+Ms.*tau*k_init/k_rec.*Hr(aa,K_tau,1));
plot(x,ss_gamma,'o','linewidth',2);
p = polyfit(x,ss_gamma,1);
hold on
plot(x,polyval(p,x),'r-','linewidth',2)
xlabel('M_v (%)');
ylabel('\gamma_M');
axis([0 12 0.2 1.0])
set(gca,'fontsize',12);
title('Mv vs gamma');

%% steady state M fraction
figure;
plot(x,(Ms-M_vs)./Ms,'o','linewidth',2);
p = polyfit(x,(Ms-M_vs)./Ms,1);
hold on
plot(x,polyval(p,x),'r-','linewidth',2)
xlabel('M_v (%)');
ylabel('(M_r+M_p)/{M}_T');
axis([0 12 0.2 1])
set(gca,'fontsize',12);
title('Mv vs (M_r+M_p)/{M}_T');

%% steady state mu
figure;
ssmu = ss(6,:);
ssmu = ssmu*3600; % convert to /h
plot(x,ssmu,'o','linewidth',2);
p = polyfit(x,ssmu,1);
hold on
plot(x,polyval(p,x),'r-','linewidth',2)
xlabel('M_v (%)');
ylabel('\mu (/hour)');
axis([0 12 0 1.4])
set(gca,'fontsize',12);
title('Mv vs mu');

%% ODE
function output = ode_fun(t,y)
    a = max(1e-1,y(1));
    P_r = max(0,y(2));
    P_v = max(0,y(3));
    rho = max(0,y(4));
    vio = max(0,y(5));
    gamma = M*tau*k_init*L/k_trans*Ha(a,K_tau,1)/(K_M+k_init/k_on*tau+M+M*tau*k_init*L/k_trans*Ha(a,K_tau,1)+M*tau*k_init/k_rec*Hr(a,K_tau,1));
    sigma = k_trans/L*P_r*gamma;
    mu = sigma*(M_p+M_r)/M*log(2)/P;    
    dydt(1) = nu_in-nu_vio*P_v*Ha(a,K_v,1)-sigma*(L+a*log(2)/P*(M_r+M_p)/M); %a
    dydt(2) = sigma*(M_r*L/M/L_r - P_r*log(2)/P*(M_r+M_p)/M) - zeta_r*P_r; %P_r
    dydt(3) = sigma*(M_v*L/M/L_v - P_v*log(2)/P*(M_r+M_p)/M) - zeta_v*P_v;
    dydt(4) = rho*mu*(1-rho/kappa); %rho
    dydt(5) = rho*nu_vio*P_v*1e-5 - k_deg*vio; %vio
    output = dydt';
end
    
end