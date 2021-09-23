%% Plots
clear all;close all;

%% KCS
Lpp = 4.367;
rho = 1025;

%% Loop over all depths
Depths = [499 416 333]; %mm

for i = 1:3
    % Load .mat data file
    load(sprintf('Data_KVLCC2_FHR_PMM_Shallow_%3.0fmm.mat',Depths(i)));
    
    %% SWAY
    
    U_s = PMMY2(:,4); % Velocity for harmonic sway
    v_sdot = PMMY2(:,8);
    Y_s = PMMY2(:,11);% Sway force
    N_s = PMMY2(:,12);% Yaw force
    K_s = PMMY2(:,24);% Drift force
    Drift_sangle = PMMY2(:,25).*pi/180;
    
    % Non-dim
    Acc_s = U_s.^2/Lpp;
    Force_s = 0.5*rho*Lpp^2.*U_s.^2;
    Mom_s = 0.5*rho*Lpp^3.*U_s.^2;
    
    v_sdot_ND = v_sdot./Acc_s;
    Y_sND = Y_s./Force_s;
    N_sND = N_s./Mom_s;
    K_sND = K_s./Mom_s;
    
    %% YAW
    
    U_y = PMMPSI2(:,4); % Velocity for harmonic sway
    k = 1:25:length(U_y);
    for j = 1:length(k)
        Vel_y(j) = U_y(k(j));
    end
    
    Y_y = PMMPSI2(:,11);
    N_y = PMMPSI2(:,12);
    K_y = PMMPSI2(:,24);
    Drift_yangle = PMMPSI2(:,25).*pi/180;
    r_y = PMMPSI2(:,6).*pi/180;
    r_ydot = PMMPSI2(:,9).*pi/180;
    
    % Non-dim
    Ang_yvel = U_y./Lpp;
    Ang_yacc = U_y.^2./Lpp^2;
    Force_y = 0.5*rho*Lpp^2.*U_y.^2;
    Mom_y = 0.5*rho*Lpp^3.*U_y.^2;
    
    r_yND = r_y./Ang_yvel;
    r_ydot_ND = r_ydot./Ang_yacc;
    Y_yND = Y_y./Force_y;
    N_yND = N_y./Mom_y;
    K_yND = K_y./Mom_y;
    
    
    %% Plot SWAY
    figure()
    plot(v_sdot_ND,N_sND*1e3,'.')
    grid on
    title(sprintf('KVLCC2 - Yaw moment as function of acceleration - Depth = 0.%3.0f m',Depths(i)))
    legend(sprintf('Velocity = %1.3f m/s \n Drift angles = %1.3f^o, %1.3f^o',U_s(1), max(Drift_sangle(1:24)),max(Drift_sangle(25:48))'),'location','southoutside' )
    xlabel('vdot´ [-]')
    ylabel('N´*10^3 [-]')
    
    
    figure()
    plot(v_sdot_ND,Y_sND*1e3,'.')
    grid on
    title(sprintf('KVLCC2 - Sway force as function of acceleration - Depth = 0.%3.0f m',Depths(i)))
    legend(sprintf('Velocity = %1.3f m/s \n Drift angles = %1.3f^o, %1.3f^o',U_s(1), max(Drift_sangle(1:24)),max(Drift_sangle(25:48))'),'location','southoutside' )
    xlabel('vdot´ [-]')
    ylabel('Y´*10^3 [-]')
    
    figure()
    plot(v_sdot_ND,K_sND*1e3,'.')
    grid on
    title(sprintf('KVLCC2 - Drift moment as function of acceleration - Depth = 0.%3.0f m',Depths(i)))
    legend(sprintf('Velocity = %1.3f m/s \n Drift angles = %1.3f^o, %1.3f^o',U_s(1), max(Drift_sangle(1:24)),max(Drift_sangle(25:48))'),'location','southoutside' )
    xlabel('vdot´ [-]')
    ylabel('K´*10^3 [-]')
    
    %% Plot YAW
    figure()
    plot(r_yND,N_yND*1e3,'o')
    grid on
    title(sprintf('KVLCC2 - Yaw moment as function of angular velocity - Depth = 0.%3.0f m',Depths(i)))
    legend(sprintf('U = %1.3f m/s \n U = %1.3f m/s',Vel_y(1:5),Vel_y(6:11)),'location','southoutside' )
    xlabel('r´ [-]')
    ylabel('N´*10^3 [-]')
    
    figure()
    plot(r_yND,Y_yND*1e3,'x')
    grid on
    title(sprintf('KVLCC2 - Sway force as function of angular velocity - Depth = 0.%3.0f m',Depths(i)))
    legend(sprintf('U = %1.3f m/s \n U = %1.3f m/s',Vel_y(1:5),Vel_y(6:11)),'location','southoutside' )
    xlabel('r´ [-]')
    ylabel('Y´*10^3 [-]')
    
    figure()
    plot(r_yND,K_yND*1e3,'.')
    grid on
    title(sprintf('KVLCC2 - Drift moment as function of angular velocity - Depth = 0.%3.0f m',Depths(i)))
    legend(sprintf('U = %1.3f m/s \n U = %1.3f m/s',Vel_y(1:5),Vel_y(6:11)),'location','southoutside' )
    xlabel('r´ [-]')
    ylabel('K´*10^3 [-]')
    
    figure()
    plot(r_ydot_ND,N_yND*1e3,'o')
    grid on
    title(sprintf('KVLCC2 - Yaw moment as function of angular acceleration - Depth = 0.%3.0f m',Depths(i)))
    legend(sprintf('U = %1.3f m/s \n U = %1.3f m/s',Vel_y(1:5),Vel_y(6:11)),'location','southoutside' )
    xlabel('rdot´ [-]')
    ylabel('N´*10^3 [-]')
    
    figure()
    plot(r_ydot_ND,Y_yND*1e3,'x')
    grid on
    title(sprintf('KVLCC2 - Sway force as function of angular acceleration - Depth = 0.%3.0f m',Depths(i)))
    legend(sprintf('U = %1.3f m/s \n U = %1.3f m/s',Vel_y(1:5),Vel_y(6:11)),'location','southoutside' )
    xlabel('rdot´ [-]')
    ylabel('Y´*10^3 [-]')
    
    figure()
    plot(r_ydot_ND,K_yND*1e3,'.')
    grid on
    title(sprintf('KVLCC2 - Drift moment as function of angular acceleration - Depth = 0.%3.0f m',Depths(i)))
    legend(sprintf('U = %1.3f m/s \n U = %1.3f m/s',Vel_y(1:5),Vel_y(6:11)),'location','southoutside' )
    xlabel('rdot´ [-]')
    ylabel('K´*10^3 [-]')
    
    
end