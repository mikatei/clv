function vio_model_cLV
close all
clc

% simulation time
Tmax = 3600*36; % sec

% parameter definition
mu_E0 = 1.7522;
mu_B0 = 2.16;
p1 = 0.3952;
p2 = 0.32;
Kc = 0.0250;
cmax = 1115; % ng/mL
b0 = 0.97;
Mvs = linspace(0,1,10);

kappa = 1e3;
ep_vio = cmax/kappa;

% seeding ratios
ec = [.01,.01,.2,.2,.01,0];
bs = [.01,.2,.01,.2,0,.01];

figure('Position',[0 0 1200 800],'Name','Simulation result of Ec Bs populations')
hold all
figure('Position',[0 0 1200 800],'Name','Simulation result of Ec Bs ratios')
hold all

for j = 1:6
    Ys = nan(2,length(Mvs));
    for i = 1:length(Mvs)
        Mv = Mvs(i);

        E_init = kappa*ec(j);
        B_init = kappa*bs(j);
        [T,Y] = ode45(@ode_fun,[0 Tmax],[E_init B_init]);
        Ys(:,i) = Y(end,:);
    end
    
    % Ec vs Bs
    figure(1)
    subplot(2,3,j)
    plot(Mvs,Ys,'-o','markers',3,'linewidth',2);
    set(gca,'yscale','log','fontsize',18)
    legend('Ec','Bs')
    axis([-0.1,1.1,kappa/400 kappa*1.1])
    axis square
    xlabel('M_v')
    ylabel('Density (cells/mL)')
    title(sprintf('Ec %.2e, Bs %.2e',ec(j),bs(j)))
    set(gca,'fontsize',18')
    set(gcf,'PaperPositionMode','auto')

    % Ratio
    figure(2)
    subplot(2,3,j)
    plot(Mvs,Ys(1,:)./kappa,'-o','markers',3,'linewidth',2);
    % set(gca,'yscale','log')
    xlabel('M_v')
    ylabel('vioABE fraction')
    axis([-0.1 1.1 0 1.1])
    axis square
    title(sprintf('Ec %.2e, Bs %.2e',ec(j),bs(j)))
    set(gca,'fontsize',18')
    set(gcf,'PaperPositionMode','auto')
end

function output = ode_fun(t,y)
    dydt = zeros(size(y));      
    rho_E = y(1);
    rho_B = y(2);
    c = ep_vio*Mv*rho_E;
    mu_E = mu_E0*(1-p1*Mv);
    mu_B = mu_B0/(1+c/cmax/Kc);
    dydt(1) = mu_E*rho_E*(1- ((b0+p2*Mv)*rho_E + rho_B)/kappa); % rho_E
    dydt(2) = mu_B*rho_B*(1- ((b0+p2*Mv)*rho_E + rho_B)/kappa); % rho_B
    output = dydt;
end
end
