close all
clear all
clc

%% Parameters

    Ts = 0.01;            % sample time
    Ttotal = 30;          % trajectory time
    t = 0:Ts:Ttotal;      % time vector
    
    N = numel(t);         % horizonte
    parms = 1;            % parâmetro que caracteriza o valor de mi
                          % parms = 0 -> mu = 0 / parms = 1 -> mu ~= 0
    mi = 1e+8;            % parâmetro de penalidade ( > 0)
    alfa = 0.01;          % pertence ao parâmetro de minimização ( > 0)
                          % lambda = 1 + alfa
                          
%% System Data 

    model; % system's parameters and model

% Variaveis globais
    warning off
    global m n q F G Q R H Ef Eg P0 F_inc G_inc un_rate
    warning  on

%% Reference Calculation

    t_lc = 0:Ts:10; % total time for lane change maneuvers
    
    u_ref = zeros(1,numel(t)); %creating u_ref vector 
    u_ref(:,1000:1500) = 0.01*sin(0.4*pi*t_lc(1:501)); %turn left
    u_ref(:,2000:2500) = -0.01*sin(0.4*pi*t_lc(501:1001)); %turn righ
    
    u_ref = [u_ref/2;u_ref/2];

    x_ref(:,1)=[0;0;0;0;0;0];
%--------------------------------------------------------------------------
    psi_ref(1) = x_ref(6,1);
    posey_ref(1) = x_ref(5,1);
    posex_ref(1) = 0;

for k=1:numel(t)-1
    
    x_ref(:,k+1)=F*x_ref(:,k)+G*u_ref(:,k);
    
    psi_ref(k+1) = psi_ref(k)+x_ref(2,k)*Ts;
    velX_ref(k) = v - x_ref(1,k)*psi_ref(k+1);
    velY_ref(k) = x_ref(1,k) + v*psi_ref(k+1);
    posey_ref(k+1) = posey_ref(k) + velY_ref(k)*Ts;
    posex_ref(k+1) = posex_ref(k) + velX_ref(k)*Ts;
    x_ref(1:4,k) = [0;0;0;0];
end
 
%input disturbance

    t_dis = 0:Ts:2;         %disturbance time 
    u_dis = zeros(1,numel(t)); %creating u_dis vector
    u_dis(:,1600:1800) = (exp(1)).^(-2*t_dis)*0.2.*sin(2*pi*t_dis);
    u_dis = [u_dis/2;u_dis/2];

% Iinitial conditions

     x(:,1) = x0;
     
% minimum and maximum steering
lim_inf = -0.22; % maximum steering angle ... remember, it is the half of the actual value
lim_sup = 0.22;  % minimum steering angle
max_rate = 100.926; % maximum steering rate
          
%% Robust Recursive Regulator
 
   [P,K,L] = SL_Regulador_Robusto_alternativo(N,parms,mi,alfa);

%% H_infinity
 
   [u_h,x_h,erro_h] = h_infinity(x0,N,x_ref,u_dis,lim_sup,lim_inf,Ts,max_rate);
 
%% Real System Simulation
 
%--------------------------------------------------------------------------
    psi(1) = x(6,1);
    posey(1) = x(5,1);
    posex(1) = 0;
   
%--------------------------------------------------------------------------

    psi_h(1) = x_h(6,1);
    posey_h(1) = x_h(5,1);
    posex_h(1) = 0;
    
%--------------------------------------------------------------------------

% Robust Linear Quadratic Regulator - with uncertainties 
for k = 1:N-1
    
    erro(:,k) = x(:,k) - x_ref(:,k);
    u(:,k) = K(:,:,k)*erro(:,k);
    u(1,k) = min(max(lim_inf,u(1,k)),lim_sup);
    u(2,k) = min(max(lim_inf,u(2,k)),lim_sup);
    
    if(k>1 && ((u(1,k)-u(1,k-1))/Ts)>max_rate)
        u(1,k) = Ts*max_rate+u(1,k-1);
        u(2,k) = Ts*max_rate*u(2,k-1);
    elseif(k>1 && (((u(1,k)-u(1,k-1))/Ts)<-max_rate))
        u(1,k) = Ts*(-max_rate)+u(1,k-1);
        u(2,k) = Ts*(-max_rate)+u(2,k-1);
    end
        
    x(:,k+1) = (F_inc)*x(:,k) + (G_inc)*u(:,k);
   
        psi(k+1) = psi(k)+x(2,k)*Ts;
        velX(k) = v - x(1,k)*psi(k+1);
        velY(k) = x(1,k) + v*psi(k+1);
        posey(k+1) = posey(k) + velY(k)*Ts;
        posex(k+1) = posex(k) + velX(k)*Ts;
        
% findig the position fot H_infinity controller

        psi_h(k+1) = psi_h(k)+x_h(2,k)*Ts;
        velX_h(k) = v - x_h(1,k)*psi_h(k+1);
        velY_h(k) = x_h(1,k) + v*psi_h(k+1);
%         posey_h(k+1) = posey_h(k) + velY_h(k)*Ts;
        posey_h(k+1) = x_h(5,k+1);
        posex_h(k+1) = posex_h(k) + velX_h(k)*Ts;

end

%% Performance metrics

%maximum variation in control input 
max_st_rate_u = max(abs(diff(u(1,:))/Ts))
max_st_rate_h = max(abs(diff(u_h(1,:))/Ts))

%norm L2 for Robust Linear Quadratic Regulator

sum_y = 0;
sum_yaw = 0;
sum_articulation_rate = 0;
sum_articulation = 0;
sum_displacement = 0;
sum_orientation = 0;

for k = 1:N-1
sum_y = sum_y + (norm(erro(1,k)))^2;
sum_yaw = sum_yaw + (norm(erro(2,k)))^2;
sum_articulation_rate = sum_articulation_rate + (norm(erro(3,k)))^2;
sum_articulation = sum_articulation + (norm(erro(4,k)))^2;
sum_displacement = sum_displacement + (norm(erro(5,k)))^2;
sum_orientation = sum_orientation + (norm(erro(6,k)))^2;
end

L2_y = sqrt(sum_y/Ttotal);
L2_yaw = sqrt(sum_yaw/Ttotal);
L2_articulation_rate = sqrt(sum_articulation_rate/Ttotal);
L2_articulation = sqrt(sum_articulation/Ttotal);
L2_displacement = sqrt(sum_displacement/Ttotal);
L2_orientation = sqrt(sum_orientation/Ttotal);

% norm L2 for H_infinity controller

sum_y_h = 0;
sum_yaw_h = 0;
sum_articulation_rate_h = 0;
sum_articulation_h = 0;
sum_displacement_h = 0;
sum_orientation_h = 0;

for k = 1:N-1
    sum_y_h = sum_y_h + (norm(erro_h(1,k)))^2;
    sum_yaw_h = sum_yaw_h + (norm(erro_h(2,k)))^2;
    sum_articulation_rate_h = sum_articulation_rate_h + (norm(erro_h(3,k)))^2;
    sum_articulation_h = sum_articulation_h + (norm(erro_h(4,k)))^2;
    sum_displacement_h = sum_displacement_h + (norm(erro_h(5,k)))^2;
    sum_orientation_h = sum_orientation_h + (norm(erro_h(6,k)))^2;
end

L2_y_h = sqrt(sum_y_h/Ttotal);
L2_yaw_h = sqrt(sum_yaw_h/Ttotal);
L2_articulation_rate_h = sqrt(sum_articulation_rate_h/Ttotal);
L2_articulation_h = sqrt(sum_articulation_h/Ttotal);
L2_displacement_h = sqrt(sum_displacement_h/Ttotal);
L2_orientation_h = sqrt(sum_orientation_h/Ttotal);

%% Plots
 
% % Control input
%     figure(1)
%     hold on
%     plot(t,2*u_h(1,1:end-1),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
%     plot(t(1:end-1),2*u(1,:),'-.k','LineWidth',1.5)
%     legend('H_{\infty}','RLQR')
%     xlabel('Time (s)','Interpreter','latex')
%     ylabel('$\alpha$ (rad)','Interpreter','latex')
%     title('\textbf{Control input $\alpha$}','Interpreter','latex')
%     grid on
% 
% % Lateral velocity
%     figure(2)
%     hold on
%     plot(t,x_h(1,1:numel(t)),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
%     plot(t,x(1,:),'-.k','LineWidth',1.5)
%     legend('H_{\infty}','RLQR')
%     xlabel('Time (s)','Interpreter','latex')
%     ylabel('$\dot{y}_{1}$ (m/s)','Interpreter','latex')
%     title('\textbf{Lateral Velocity ($\dot{y}_{1}$)}','Interpreter','latex')
%     grid on
% 
% % Yaw rate
%     figure(3)
%     hold on
%     plot(t,x_h(2,1:numel(t)),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
%     plot(t,x(2,:),'-.k','LineWidth',1.5)
%     legend('H_{\infty}','RLQR')
%     xlabel('Time (s)','Interpreter','latex')
%     ylabel('$\dot{\psi}$ (rad/s)','Interpreter','latex')
%     title('\textbf{Yaw rate ($\dot{\psi}$)}','Interpreter','latex')
%     grid on
% 
% % Articulation rate
%     figure(4)
%     hold on
%     plot(t,x_h(3,1:numel(t)),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
%     plot(t,x(3,:),'-.k','LineWidth',1.5)
%     legend('H_{\infty}','RLQR')
%     xlabel('Time (s)','Interpreter','latex')
%     ylabel('$\dot{\phi}$ (rad/s)','Interpreter','latex')
%     title('\textbf{Articulation angle rate ($\dot{\phi}$)}','Interpreter','latex')
%     grid on
% 
% % Articulation angle
%     figure(5)
%     hold on
%     plot(t,x_h(4,1:numel(t)),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
%     plot(t,x(4,:),'-.k','LineWidth',1.5)
%     legend('H_{\infty}','RLQR')
%     xlabel('Time (s)','Interpreter','latex')
%     ylabel('$\phi$ (rad)','Interpreter','latex')
%     title('\textbf{Articulation Angle ($\phi$)}','Interpreter','latex')
%     grid on
% 
% % Displacement error
%     figure(6)
%     hold on
%     plot(t(1:end-1),erro_h(5,1:3000),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
%     plot(t(1:end-1),erro(5,:),'-.k','LineWidth',1.5)
%     legend('H_{\infty}','RLQR')
%     xlabel('Time (s)','Interpreter','latex')
%     ylabel('$\rho$ (m)','Interpreter','latex')
%     title('\textbf{Displacement error ($\rho$)}','Interpreter','latex')
%     grid on
% 
% % Orientation error
%     figure(7)
%     hold on
%     plot(t(1:end-1),erro_h(6,1:3000),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
%     plot(t(1:end-1),erro(6,:),'-.k','LineWidth',1.5)
%     legend('H_{\infty}','RLQR')
%     xlabel('Time (s)','Interpreter','latex')
%     ylabel('$\theta$ (rad)','Interpreter','latex')
%     title('\textbf{Orientation error ($\theta$)}','Interpreter','latex')
%     grid on
% 
% Global position
    figure(8)
    hold on
    plot(posex_ref,posey_ref,':','Color',[0.33,0.33,0.33],'LineWidth',1.5)
    plot(posex_h,posey_h,'Color',[0.5,0.5,0.5],'LineWidth',1.5)
    plot(posex,posey,'-.k','LineWidth',1.5)
    axis([posex(1) posex(end) -0.5 6]);
    annotation('arrow',[0.15 0.21],[0.2 0.325]) %put an arrow
%     annotation('arrow',[0.5 0.45],[0.5 0.6]) %put an arrow
%     annotation('arrow',[0.69 0.76],[0.5 0.58]) %put an arrow
    legend({'Reference','H_{\infty}','RLQR'},'FontSize',13)
    xlabel('X (m)','Interpreter','latex','FontSize',14.5)
    ylabel('Y (m)','Interpreter','latex','FontSize',14.5)
    title('\textbf{Global Position}','Interpreter','latex','FontSize',14.5)
    grid on
    %begining
    axes('position',[.16 .35 .1 .15])%[positionX positionY height width]
    box on
    indexOfInterest = (posex_ref >0) & (posex_ref<40); %interval
    hold on
    plot(posex_ref(indexOfInterest),posey_ref(indexOfInterest),':','Color',[0.33,0.33,0.33],'LineWidth',1.5)
    plot(posex_h(indexOfInterest),posey_h(indexOfInterest),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
    plot(posex(indexOfInterest),posey(indexOfInterest),'-.k','LineWidth',1.5)
    grid on
    axis tight
%     %first lane change
%     axes('position',[.31 .61 .15 .15]) %[positionX positionY width height]
%     box on
%     indexOfInterest2 = (posex_ref >215) & (posex_ref<275);%interval
%     hold on
%     plot(posex_ref(indexOfInterest2),posey_ref(indexOfInterest2),':','Color',[0.33,0.33,0.33],'LineWidth',1.5)
%     plot(posex_h(indexOfInterest2),posey_h(indexOfInterest2),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
%     plot(posex(indexOfInterest2),posey(indexOfInterest2),'-.k','LineWidth',1.5)
%     grid on
%     axis tight
%     %second lane change
%     axes('position',[.74 .61 .15 .15]) %[positionX positionY width height]
%     box on
%     indexOfInterest3 = (posex_ref >380) & (posex_ref<460);%interval
%     hold on
%     plot(posex_ref(indexOfInterest3),posey_ref(indexOfInterest3),':','Color',[0.33,0.33,0.33],'LineWidth',1.5)
%     plot(posex_h(indexOfInterest3),posey_h(indexOfInterest3),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
%     plot(posex(indexOfInterest3),posey(indexOfInterest3),'-.k','LineWidth',1.5)
%     grid on
%     axis tight
    
    
% Subplot state vector

    figure(10)
    subplot(3,2,1) 
    hold on
    plot(t,x_h(1,1:numel(t)),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
    plot(t,x(1,:),'-.k','LineWidth',1.5)
%     axis([t(1) t(end) -3 3])
    legend({'H_{\infty}','RLQR'},'Location','NorthEast','FontSize',13)
    ylabel('$\dot{y}_{1}$ (m/s)','Interpreter','latex','FontSize',14.5)
    title('\textbf{Lateral Velocity ($\dot{y}_{1}$)}','Interpreter','latex','FontSize',14.5)
    grid on
    
    subplot(3,2,2)
    hold on
    plot(t,x_h(2,1:numel(t)),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
    plot(t,x(2,:),'-.k','LineWidth',1.5)
%     axis([t(1) t(end) -2 2])
    legend({'H_{\infty}','RLQR'},'Location','NorthEast','FontSize',13)
    ylabel('$\dot{\psi}$ (rad/s)','Interpreter','latex','FontSize',14.5)
    title('\textbf{Yaw rate ($\dot{\psi}$)}','Interpreter','latex','FontSize',14.5)
    grid on
    
    subplot(3,2,3)
    hold on
    plot(t,x_h(3,1:numel(t)),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
    plot(t,x(3,:),'-.k','LineWidth',1.5)
%     axis([t(1) t(end) -3.5 3.5])
    legend({'H_{\infty}','RLQR'},'Location','NorthEast','FontSize',13)
    ylabel('$\dot{\phi}$ (rad/s)','Interpreter','latex','FontSize',14.5)
    title('\textbf{Articulation angle rate ($\dot{\phi}$)}','Interpreter','latex','FontSize',14.5)
    grid on
   
    subplot(3,2,4)
    hold on
    plot(t,x_h(4,1:numel(t)),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
    plot(t,x(4,:),'-.k','LineWidth',1.5)
    legend({'H_{\infty}','RLQR'},'Location','NorthEast','FontSize',13)
    ylabel('$\phi$ (rad)','Interpreter','latex','FontSize',14.5)
    title('\textbf{Articulation Angle ($\phi$)}','Interpreter','latex','FontSize',14.5)
    grid on
   
    subplot(3,2,5)
    hold on
    plot(t(1:end-1),erro_h(5,1:3000),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
    plot(t(1:end-1),erro(5,:),'-.k','LineWidth',1.5)
    legend({'H_{\infty}','RLQR'},'Location','NorthEast','FontSize',13)
    xlabel('Time (s)','Interpreter','latex','FontSize',14.5)
    ylabel('$\rho$ (m)','Interpreter','latex','FontSize',14.5)
    title('\textbf{Displacement error ($\rho$)}','Interpreter','latex','FontSize',14.5)
    grid on
    
    subplot(3,2,6)
    hold on
    plot(t(1:end-1),erro_h(6,1:3000),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
    plot(t(1:end-1),erro(6,:),'-.k','LineWidth',1.5)
    legend({'H_{\infty}','RLQR'},'Location','NorthEast','FontSize',13)
    xlabel('Time (s)','Interpreter','latex','FontSize',14.5)
    ylabel('$\theta$ (rad)','Interpreter','latex','FontSize',14.5)
    title('\textbf{Orientation error ($\theta$)}','Interpreter','latex','FontSize',14.5)
    grid on
    
    %subplot control input and global position
    figure(11)
    subplot(2,1,1)
    hold on
    plot(posex_ref,posey_ref,':','Color',[0.33,0.33,0.33],'LineWidth',1.5)
    plot(posex_h,posey_h,'Color',[0.5,0.5,0.5],'LineWidth',1.5)
    plot(posex,posey,'-.k','LineWidth',1.5)
    annotation('arrow',[0.15 0.185],[0.625 0.678])
    axis([posex(1) posex(end) -0.5 6]);
    legend({'Reference','H_{\infty}','RLQR'},'Location','NorthEast','FontSize',14)
    xlabel('X (m)','Interpreter','latex','FontSize',16)
    ylabel('Y (m)','Interpreter','latex','FontSize',16)
    title('\textbf{Global Position}','Interpreter','latex','FontSize',16)
    grid on
    %begining
    axes('position',[.16 .7 .08 .1])%[positionX positionY height width]
    box on
    indexOfInterest = (posex_ref >0) & (posex_ref<40); %interval
    hold on
    plot(posex_ref(indexOfInterest),posey_ref(indexOfInterest),':','Color',[0.33,0.33,0.33],'LineWidth',1.5)
    plot(posex_h(indexOfInterest),posey_h(indexOfInterest),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
    plot(posex(indexOfInterest),posey(indexOfInterest),'-.k','LineWidth',1.5)
    grid on
    axis tight
%     %first lane change
%     axes('position',[.35 .77 .08 .15]) %[positionX positionY width height]
%     box on
%     indexOfInterest2 = (posex_ref >215) & (posex_ref<245);%interval
%     hold on
%     plot(posex_ref(indexOfInterest2),posey_ref(indexOfInterest2),'Color',[0.2,0.8,0.2],'LineWidth',1.5)
%     plot(posex(indexOfInterest2),posey(indexOfInterest2),'r','LineWidth',1.5)
%     plot(posex_h(indexOfInterest2),posey_h(indexOfInterest2),':b','LineWidth',1.5)
%     grid on
%     axis tight
%     %second lane change
%     axes('position',[.73 .77 .08 .15]) %[positionX positionY width height]
%     box on
%     indexOfInterest3 = (posex_ref >394) & (posex_ref<445);%interval
%     hold on
%     plot(posex_ref(indexOfInterest3),posey_ref(indexOfInterest3),'Color',[0.2,0.8,0.2],'LineWidth',1.5)
%     plot(posex(indexOfInterest3),posey(indexOfInterest3),'r','LineWidth',1.5)
%     plot(posex_h(indexOfInterest3),posey_h(indexOfInterest3),':b','LineWidth',1.5)
%     grid on
%     axis ([394 445 -0.1 0.6498])
    
    subplot(2,1,2)
    hold on
    plot(t,2*u_h(1,1:end-1),'Color',[0.5,0.5,0.5],'LineWidth',1.5)
    plot(t(1:end-1),2*u(1,:),'-.k','LineWidth',1.5)
%     axis([t(1) t(end) -1 1])
    legend({'H_{\infty}','RLQR'},'Location','NorthEast','FontSize',14)
    xlabel('Time (s)','Interpreter','latex','FontSize',16)
    ylabel('$\alpha$ (rad)','Interpreter','latex','FontSize',16)
    title('\textbf{Control input $\alpha$}','Interpreter','latex','FontSize',16)
    grid on