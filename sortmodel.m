% %WT P91 parameters on bead
% kcat = 0.076; %s^-1
% Km = 44.4*10^-6; %M
% Ks = 136*10^-6; %M
% 
% %WT P91 parameters in solution
% kcat = 0.099; %s^-1
% Km = 26.9*10^-6; %M
Ks = 577*10^-6; %M


%initial conditions (in M)
E_init = 0.006*10^-6;  
C_init = 0;
P_init = 0;
CI_init = 0;

i=1;
cov=100; % set coverage of screen
time1=780; %set beginning of the screen
time2=1020; %set end of the screen

%simulate variants with different kcat and Km
% for S_init = (1*10^-6):(20*10^-6):(100*10^-6)
%     for kcat_v=(kcat/50):(2*kcat):(kcat*50) 
%         for Km_v=(Km/50):(2*Km):(Km*50)
for S_init = [(26.9*10^-6)]
    for kcat_v=logspace(-3,0,20)
        for Km_v=logspace(-6,-3,20)


            %calculate rate constants (in min^-1)
            k1 = (10^8)*60; %min^-1 
            k2 = kcat_v*60; %min^-1
            k_1 = Km_v*k1-k2; 
            k3 = (10^8)*60; 
            k_3 = Ks*k3;

            %solve IVP
            options = odeset('RelTol',1e-13);
            x_init = [E_init;S_init;C_init;P_init;CI_init];
            [t,x] = ode15s(@(t,y) enz_kin_inh(t,y,k1,k_1,k2,k3,k_3),[0,1020],x_init, options); 
            sol = ode15s(@(t,y) enz_kin_inh(t,y,k1,k_1,k2,k3,k_3),[time1,time2],x_init);  %time interval defines two extreme timepoints of sorting
            
            %retrieve value of P at two extreme points of sort duration,
            %at a randomly drawed timepoint (P_rand), and at set of random
            %timepoints (P_rand_50) with coverage cov
            kcat_o{i} = kcat_v;
            Km_o{i} = Km_v;
            P_t1{i} = (sol.y(4,2));  %
            P_t2{i} = (sol.y(4,end));
            P_mean{i} = ((sol.y(4,2))+(sol.y(4,end)))./2;
            P_rand{i} = ((sol.y(4,end))-(sol.y(4,2))).*rand(1,1)+(sol.y(4,2));
            P_rand_50{i} = ((sol.y(4,end))-(sol.y(4,2))).*rand(cov,1)+(sol.y(4,2));
            S_in{i} = S_init;
            i=i+1;
        end
    end
end

%convert to numerical
P_t1_mat = cell2mat(P_t1(:,:));
P_t2_mat = cell2mat(P_t2(:,:));
S_in_mat = cell2mat(S_in(:,:));
P_mean_mat = cell2mat(P_mean(:,:));
%P_rand =((P_t2_mat-P_t1_mat).*rand(1,1)+P_t1_mat); 
Km_o_mat = cell2mat(Km_o(:,:));
kcat_o_mat = cell2mat(kcat_o(:,:));
P_rand_mat = cell2mat(P_rand(:,:));
P_rand_50_mat = cell2mat(P_rand_50(:,:));

P_rand_avg = mean(P_rand_50_mat); %calculate average of each variant

Out_summ = [Km_o_mat; kcat_o_mat; P_rand_mat; S_in_mat]; %results array for 1 droplet/variant
Out_summ_avg = [Km_o_mat; kcat_o_mat; P_rand_avg; S_in_mat]; %results array for set coverage
Out_summ_mean = [Km_o_mat; kcat_o_mat; P_mean_mat; S_in_mat]; %results array for a single (mean) incubation time 


% %split into separate cells of 3D array based on S_init
Out_summ(5,:) = categorical(Out_summ(4,:));
n=unique(Out_summ(5,:));
for j = n
    ind(j,:) = Out_summ(5,:) == j;
    Out_summ_grouped{j} = Out_summ(:,ind(j,:));
end

Out_summ_avg(5,:) = categorical(Out_summ_avg(4,:));
n=unique(Out_summ_avg(5,:));
for j = n
    ind(j,:) = Out_summ_avg(5,:) == j;
    Out_summ_grouped_avg{j} = Out_summ_avg(:,ind(j,:));
end

Out_summ_mean(5,:) = categorical(Out_summ_mean(4,:));
n=unique(Out_summ_mean(5,:));
for j = n
    ind(j,:) = Out_summ_mean(5,:) == j;
    Out_summ_grouped_mean{j} = Out_summ_mean(:,ind(j,:));
end


%plot each cell of Out_summ_grouped array
figure();
set(gca,'color','none')
for plotID=1:j
    subplot(3,j,plotID);
    loglog(Out_summ_grouped{1,plotID}(2,:), Out_summ_grouped{1,plotID}(3,:), 'dw'); title(Out_summ_grouped{1,plotID}(4,1));
    ylim([0, 100*10^-6])
    hold on
    pointsize=8;
    scatter(Out_summ_grouped{1,plotID}(2,:), Out_summ_grouped{1,plotID}(3,:), pointsize, Out_summ_grouped{1,plotID}(1,:));
    colorbar
    set(gca,'ColorScale','log')
    xlabel('kcat (1/s)')
    ylabel('[Fluorescein] (M)')
    xline(0.099)
    yline(50*10^-9, 'r')
    grid on
    
    subplot(3,j,j+plotID);
    loglog(Out_summ_grouped{1,plotID}(1,:), Out_summ_grouped{1,plotID}(3,:), 'dw'); title(Out_summ_grouped{1,plotID}(4,1));
    ylim([0, 100*10^-6])
    hold on
    pointsize=8;
    scatter(Out_summ_grouped{1,plotID}(1,:), Out_summ_grouped{1,plotID}(3,:), pointsize, Out_summ_grouped{1,plotID}(2,:));
    colorbar
    set(gca,'ColorScale','log')
    xlabel('Km (M)')
    ylabel('[Fluorescein] (M)')
    xline(26.9*10^-6)
    yline(50*10^-9, 'r')
    grid on
    
    subplot(3,j,2*j+plotID);
    scatter3(Out_summ_grouped{1,plotID}(2,:), Out_summ_grouped{1,plotID}(1,:), Out_summ_grouped{1,plotID}(3,:), 10, 'filled')
    zlim([0, 100*10^-6])
    set(gca,'Xscale','log','Zscale','log','Yscale','log')
    title(Out_summ_grouped{1,plotID}(4,1));
    xlabel('kcat (1/s)')
    zlabel('[Fluorescein] (M)')
    ylabel('Km (M)')
    grid on 
    hold on
    %add reference point
   
end   
sgtitle("Coverage = 1," + " Sort time: " + time1 + "-" + time2 + " min, [Enzyme]: " + E_init + " M" )


%plot each cell of Out_summ_grouped_avg array
figure();
set(gca,'color','none')
for plotID=1:j
    subplot(3,j,plotID);
    loglog(Out_summ_grouped_avg{1,plotID}(2,:), Out_summ_grouped_avg{1,plotID}(3,:), 'dw'); title(Out_summ_grouped_avg{1,plotID}(4,1));
    ylim([0, 100*10^-6])
    hold on
    pointsize=8;
    scatter(Out_summ_grouped_avg{1,plotID}(2,:), Out_summ_grouped_avg{1,plotID}(3,:), pointsize, Out_summ_grouped_avg{1,plotID}(1,:));
    colorbar
    set(gca,'ColorScale','log')
    xlabel('kcat (1/s)')
    ylabel('[Fluorescein] (M)')
    xline(0.099)
    yline(50*10^-9, 'r')
    grid on
    
    subplot(3,j,j+plotID);
    loglog(Out_summ_grouped_avg{1,plotID}(1,:), Out_summ_grouped_avg{1,plotID}(3,:), 'dw'); title(Out_summ_grouped_avg{1,plotID}(4,1));
    ylim([0, 100*10^-6])
    hold on
    pointsize=8;
    scatter(Out_summ_grouped_avg{1,plotID}(1,:), Out_summ_grouped_avg{1,plotID}(3,:), pointsize, Out_summ_grouped_avg{1,plotID}(2,:));
    colorbar
    set(gca,'ColorScale','log')
    xlabel('Km (M)')
    ylabel('[Fluorescein] (M)')
    xline(26.9*10^-6)
    yline(50*10^-9, 'r')
    grid on
    
    subplot(3,j,2*j+plotID);
    scatter3(Out_summ_grouped_avg{1,plotID}(2,:), Out_summ_grouped_avg{1,plotID}(1,:), Out_summ_grouped_avg{1,plotID}(3,:), 10, 'filled')
    zlim([0, 100*10^-6])
    set(gca,'Xscale','log','Zscale','log','Yscale','log')
    title(Out_summ_grouped_avg{1,plotID}(4,1));
    xlabel('kcat (1/s)')
    zlabel('[Fluorescein] (M)')
    ylabel('Km (M)')
    grid on 
    hold on
    %add reference point
end   
sgtitle("Coverage = " + cov +  " Sort time: " + time1 + "-" + time2 + " min, [Enzyme]: " + E_init + " M" )




%plot each cell of Out_summ_grouped_mean array
figure();
set(gca,'color','none')
for plotID=1:j
    subplot(3,j,plotID);
    loglog(Out_summ_grouped_mean{1,plotID}(2,:), Out_summ_grouped_mean{1,plotID}(3,:), 'dw'); title(Out_summ_grouped_mean{1,plotID}(4,1));
    ylim([0, 100*10^-6])
    hold on
    pointsize=8;
    scatter(Out_summ_grouped_mean{1,plotID}(2,:), Out_summ_grouped_mean{1,plotID}(3,:), pointsize, Out_summ_grouped_mean{1,plotID}(1,:));
    colorbar
    set(gca,'ColorScale','log')
    xlabel('kcat (1/s)')
    ylabel('[Fluorescein] (M)')
    xline(0.099)
    yline(50*10^-9, 'r')
    grid on
    
    subplot(3,j,j+plotID);
    loglog(Out_summ_grouped_mean{1,plotID}(1,:), Out_summ_grouped_mean{1,plotID}(3,:), 'dw'); title(Out_summ_grouped_mean{1,plotID}(4,1));
    ylim([0, 100*10^-6])
    hold on
    pointsize=8;
    scatter(Out_summ_grouped_mean{1,plotID}(1,:), Out_summ_grouped_mean{1,plotID}(3,:), pointsize, Out_summ_grouped_mean{1,plotID}(2,:));
    colorbar
    set(gca,'ColorScale','log')
    xlabel('Km (M)')
    ylabel('[Fluorescein] (M)')
    xline(26.9*10^-6)
    yline(50*10^-9, 'r')
    grid on
    
    subplot(3,j,2*j+plotID);
    scatter3(Out_summ_grouped_mean{1,plotID}(2,:), Out_summ_grouped_mean{1,plotID}(1,:), Out_summ_grouped_mean{1,plotID}(3,:), 10, 'filled')
    zlim([0, 100*10^-6])
    set(gca,'Xscale','log','Zscale','log','Yscale','log')
    title(Out_summ_grouped_avg{1,plotID}(4,1));
    xlabel('kcat (1/s)')
    zlabel('[Fluorescein] (M)')
    ylabel('Km (M)')
    grid on 
    hold on
    %add reference point
end   
sgtitle("Coverage = " + 1 +  " Sort time: " + (time1+time2)/2 + " min, [Enzyme]: " + E_init + " M" )





%plot each cell of Out_summ_grouped_avg array
figure();
set(gca,'color','none')
for plotID=1:j
 subplot(1,j,plotID);
    loglog(Out_summ_grouped_avg{1,plotID}(2,:)./Out_summ_grouped_avg{1,plotID}(1,:), Out_summ_grouped_avg{1,plotID}(3,:), 'dw'); title(Out_summ_grouped_avg{1,plotID}(4,1));
    ylim([0, 100*10^-6])
    hold on
    pointsize=8;
    scatter(Out_summ_grouped_avg{1,plotID}(2,:)./Out_summ_grouped_avg{1,plotID}(1,:), Out_summ_grouped_avg{1,plotID}(3,:), pointsize, Out_summ_grouped_avg{1,plotID}(2,:));
    colorbar
    c.Label.String = 'kcat (1/s)';
    set(gca,'ColorScale','log')
    xlabel('kcat/Km(1/(Ms))')
    ylabel('[Fluorescein] (M)')
    xline(3668)
    yline(50*10^-9, 'r')
    grid on
   
end   
sgtitle("Coverage = " + cov +  " Sort time: " + time1 + "-" + time2 + " min, [Enzyme]: " + E_init + " M" )

