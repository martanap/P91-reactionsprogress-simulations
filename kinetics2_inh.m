
%% Define parameters

% 
% kcat = 0.099 %s^-1
% Km = 26.9*10^-6 %M
% Ks = 577*10^-6 %M
% k_1 = 1; 
% k1 = kcat*60/(Km)+(1/(Km))*k_1;  
% k2 = kcat*60; %min^-1; kcat*60
% k3 = 1;
% k_3 = Ks*k3;



kcat = 0.099 %s^-1
Km = 26.9*10^-6 %M
Ks = 577*10^-6 %M
k1 = (10^6)*60; %min^-1 
k2 = kcat*60; %min^-1; kcat*60
k_1 = Km*k1-k2; %min^-1
k3 = (10^6)*60; %min^-1
k_3 = Ks*k3;


%% Initial conditions
E_init = 0.1*10^-6;  %M
S_init = 82*10^-6; %M
C_init = 0;
P_init = 0;
CI_init = 0;

x_init = [E_init;S_init;C_init;P_init;CI_init];

%% Solve IVP
[t,x] = ode23s(@(t,y) enz_kin_inh(t,y,k1,k_1,k2,k3,k_3),[0,500],x_init);

%% plot
figure
subplot(1,2,1)
plot(t,x(:,4),'-',t,x(:,2),'-','LineWidth',1)
ax = gca;
ax.FontSize = 12;
xlabel('Time (min)')
ylabel('Concentration (M)')
ax.YLim = [0-S_init/100,S_init+S_init/100];
legend('Product', 'Substrate')
grid on

subplot(1,2,2)
plot(t,x(:,1),'-',t,x(:,3),'LineWidth',1)
ax = gca;
ax.FontSize = 12;
ax.YLim = [0-E_init/100,E_init+E_init/100];
xlabel('Time (min)')
ylabel('Concentration (M)')
legend('Enzyme', 'Enzyme:Substrate Complex')
grid on

