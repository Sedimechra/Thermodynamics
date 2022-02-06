%% ASEN 3113 Lab 2

% Group Members:
% Tristan Workman
% Alicia Wu
% Autumn Martinez
% Gerardo Romero
% Linus Shmitz
% Ty Banach

clc; clear all;
%% Question 2

% Loading in data files
Aluminum_25V = readmatrix("Aluminum_25V_240ma");
Aluminum_28V = readmatrix("Aluminum_28V_269ma");
Brass_26V = readmatrix("Brass_26V_245ma");
Brass_29V = readmatrix("Brass_29V_273ma");
Steel_21V = readmatrix("Steel_21V_192ma");

% Pre defining the total length of the rod
L = 2*0.127; % [m]

% Defining a vector for the locations of the thermocouples
thermocL = linspace(0.034925,0.123825,8);

% k values for each material (Aluminum, brass, steel)
k_AL = 130; % [W/(m*K)]
k_BR = 115; % [W/(m*K)]
k_ST = 16.2; % [W/(m*K)]

% Density for each material
density_AL = 2810; % [kg/m^3]
density_BR = 8500; % [kg/m^3]
density_ST = 8000; % [kg/m^3]

% Cp values for each material
cp_AL = 960; % [J/kg*K]
cp_BR = 380; % [J/kg*K]
cp_ST = 500; % [J/kg*K]

% calculating thermal diffusivity contant
alpha_AL = k_AL / (cp_AL * density_AL);
alpha_BR = k_BR / (cp_BR * density_BR);
alpha_ST = k_ST / (cp_ST * density_ST);

% calculating power of each voltage / current combination
AL_25V_P = 25 * 0.24;
AL_28V_P = 28 * 0.269;
BR_26V_P = 26 * 0.245;
BR_29V_P = 29 * 0.273;
ST_21V_P = 21 * 0.192;

Area = 0.000507; % [m^2]

% preallocating steady state thermocouple temperature vectors
thermoc_AL_25V = zeros(1,8);
thermoc_AL_28V = zeros(1,8);
thermoc_BR_26V = zeros(1,8);
thermoc_BR_29V = zeros(1,8);
thermoc_ST_21V = zeros(1,8);

% finding the last temperature value at each thermocouple
for i = 1:8
    thermoc_AL_25V(1,i) = Aluminum_25V(end,i+2);
    thermoc_AL_28V(1,i) = Aluminum_28V(end,i+2);
    thermoc_BR_26V(1,i) = Brass_26V(end,i+2);
    thermoc_BR_29V(1,i) = Brass_29V(end,i+2);
    thermoc_ST_21V(1,i) = Steel_21V(end,i+2);
end

% polyfitting the steady-state temperature values versus distance from
% chiller
polyAL_25V = polyfit(thermocL,thermoc_AL_25V,1);
polyAL_28V = polyfit(thermocL,thermoc_AL_28V,1);
polyBR_26V = polyfit(thermocL,thermoc_BR_26V,1);
polyBR_29V = polyfit(thermocL,thermoc_BR_29V,1);
polyST_21V = polyfit(thermocL,thermoc_ST_21V,1);

% assigning a variable to the initial temperatures (according to polyfit)
AL_25V_T0 = polyAL_25V(1,2);
AL_28V_T0 = polyAL_28V(1,2);
BR_26V_T0 = polyBR_26V(1,2);
BR_29V_T0 = polyBR_29V(1,2);
ST_21V_T0 = polyST_21V(1,2);

% calculating anyltical H values using power (calculated earlier)
H_any_AL_25V = AL_25V_P / (k_AL * Area);
H_any_AL_28V = AL_28V_P / (k_AL * Area);
H_any_BR_26V = BR_26V_P / (k_BR * Area);
H_any_BR_29V = BR_29V_P / (k_BR * Area);
H_any_ST_21V = ST_21V_P / (k_ST * Area);

% vectors containing anylitical vs experimental H values and T0 values
H_exp = [polyAL_25V(1,1) polyAL_28V(1,1) polyBR_26V(1,1) polyBR_29V(1,1) polyST_21V(1,1)];
H_any = [H_any_AL_25V H_any_AL_28V H_any_BR_26V H_any_BR_29V H_any_ST_21V];
T0 = [AL_25V_T0 AL_28V_T0 BR_26V_T0 BR_29V_T0 ST_21V_T0];

time = [1:1:12000];

% define a vector of H * x + T0 (15 values long, from x=0 to x=0.15m
% there will be 5 of these vectors, one for each material / voltage
% combination

% preallocating u(x,t) matricies
% u_AL_25V_exp = zeros(3000,8);
u_AL_25V_any = zeros(12000,8);
% u_AL_28V_exp = zeros(3000,8);
u_AL_28V_any = zeros(12000,8);
% u_BR_26V_exp = zeros(3000,8);
u_BR_26V_any = zeros(12000,8);
% u_BR_29V_exp = zeros(3000,8);
u_BR_29V_any = zeros(12000,8);
% u_ST_21V_exp = zeros(3000,8);
u_ST_21V_any = zeros(12000,8);

% defining the non summation portion of u(x,t)
% w_AL_25V_exp = thermocL*polyAL_25V(1,1) + AL_25V_T0;
w_AL_25V_any = thermocL*H_any(1,1) + AL_25V_T0;
% w_AL_28V_exp = thermocL*polyAL_28V(1,1) + AL_28V_T0;
w_AL_28V_any = thermocL*H_any(1,2) + AL_28V_T0;
% w_BR_26V_exp = thermocL*polyBR_26V(1,1) + BR_26V_T0;
w_BR_26V_any = thermocL*H_any(1,3) + BR_26V_T0;
% w_BR_29V_exp = thermocL*polyBR_29V(1,1) + BR_29V_T0;
w_BR_29V_any = thermocL*H_any(1,4) + BR_29V_T0;
% w_ST_21V_exp = thermocL*polyST_21V(1,1) + ST_21V_T0;
w_ST_21V_any = thermocL*H_any(1,5) + ST_21V_T0;


for j = 1:8
    % in the loop, for each vector at each x value, add 10 values of the
    % summation
    for i = 1:12000
%         u_AL_25V_exp(i,j) = w_AL_25V_exp(1,j);
        u_AL_25V_any(i,j) = w_AL_25V_any(1,j);
%         u_AL_28V_exp(i,j) = w_AL_28V_exp(1,j);
        u_AL_28V_any(i,j) = w_AL_28V_any(1,j);
%         u_BR_26V_exp(i,j) = w_BR_26V_exp(1,j);
        u_BR_26V_any(i,j) = w_BR_26V_any(1,j);
%         u_BR_29V_exp(i,j) = w_BR_29V_exp(1,j);
        u_BR_29V_any(i,j) = w_BR_29V_any(1,j);
%         u_ST_21V_exp(i,j) = w_ST_21V_exp(1,j);
        u_ST_21V_any(i,j) = w_ST_21V_any(1,j);
%         sum_AL_25V_exp = 0;
        sum_AL_25V_any = 0;
%         sum_AL_28V_exp = 0;
        sum_AL_28V_any = 0;
%         sum_BR_26V_exp = 0;
        sum_BR_26V_any = 0;
%         sum_BR_29V_exp = 0;
        sum_BR_29V_any = 0;
%         sum_ST_21V_exp = 0;
        sum_ST_21V_any = 0;
        for n = 1:10
%             sum_AL_25V_exp = sum_AL_25V_exp + -(8 .* H_any(1,1) .* L .* (-1).^n)./(pi.^2 .* (2 .* n - 1).^2) .* sin(j .* ((2 .* n - 1) .* pi)./(2 .* L)) .* exp(-(((2 .* n - 1) .* pi)./(2 .* L)).^2 .* alpha_AL .* i);
            sum_AL_25V_any = sum_AL_25V_any + -(8 .* H_any(1,1) .* L .* (-1).^n)./(pi.^2 .* (2 .* n - 1).^2) .* sin(j .* ((2 .* n - 1) .* pi)./(2 .* L)) .* exp(-(((2 .* n - 1) .* pi)./(2 .* L)).^2 .* alpha_AL .* i);
%             sum_AL_28V_exp = sum_AL_28V_exp + -(8 .* H_any(1,2) .* L .* (-1).^n)./(pi.^2 .* (2 .* n - 1).^2) .* sin(j .* ((2 .* n - 1) .* pi)./(2 .* L)) .* exp(-(((2 .* n - 1) .* pi)./(2 .* L)).^2 .* alpha_AL .* i);
            sum_AL_28V_any = sum_AL_28V_any + -(8 .* H_any(1,2) .* L .* (-1).^n)./(pi.^2 .* (2 .* n - 1).^2) .* sin(j .* ((2 .* n - 1) .* pi)./(2 .* L)) .* exp(-(((2 .* n - 1) .* pi)./(2 .* L)).^2 .* alpha_AL .* i);
%             sum_BR_26V_exp = sum_BR_26V_exp + -(8 .* H_any(1,3) .* L .* (-1).^n)./(pi.^2 .* (2 .* n - 1).^2) .* sin(j .* ((2 .* n - 1) .* pi)./(2 .* L)) .* exp(-(((2 .* n - 1) .* pi)./(2 .* L)).^2 .* alpha_BR .* i);
            sum_BR_26V_any = sum_BR_26V_any + -(8 .* H_any(1,3) .* L .* (-1).^n)./(pi.^2 .* (2 .* n - 1).^2) .* sin(j .* ((2 .* n - 1) .* pi)./(2 .* L)) .* exp(-(((2 .* n - 1) .* pi)./(2 .* L)).^2 .* alpha_BR .* i);
%             sum_BR_29V_exp = sum_BR_29V_exp + -(8 .* H_any(1,4) .* L .* (-1).^n)./(pi.^2 .* (2 .* n - 1).^2) .* sin(j .* ((2 .* n - 1) .* pi)./(2 .* L)) .* exp(-(((2 .* n - 1) .* pi)./(2 .* L)).^2 .* alpha_BR .* i);
            sum_BR_29V_any = sum_BR_29V_any + -(8 .* H_any(1,4) .* L .* (-1).^n)./(pi.^2 .* (2 .* n - 1).^2) .* sin(j .* ((2 .* n - 1) .* pi)./(2 .* L)) .* exp(-(((2 .* n - 1) .* pi)./(2 .* L)).^2 .* alpha_BR .* i);
%             sum_ST_21V_exp = sum_ST_21V_exp + -(8 .* H_any(1,5) .* L .* (-1).^n)./(pi.^2 .* (2 .* n - 1).^2) .* sin(j .* ((2 .* n - 1) .* pi)./(2 .* L)) .* exp(-(((2 .* n - 1) .* pi)./(2 .* L)).^2 .* alpha_ST .* i);
            % test = -(8 .* H_any(1,5) .* L .* (-1).^n)./(pi.^2 .* (2 .* n - 1).^2) .* sin(j .* ((2 .* n - 1) .* pi)./(2 .* L)) .* exp(-(((2 .* n - 1) .* pi)./(2 .* L)).^2 .* alpha_ST .* i);
            sum_ST_21V_any = sum_ST_21V_any + -(8 .* H_any(1,5) .* L .* (-1).^n)./(pi.^2 .* (2 .* n - 1).^2) .* sin(j .* ((2 .* n - 1) .* pi)./(2 .* L)) .* exp(-(((2 .* n - 1) .* pi)./(2 .* L)).^2 .* alpha_ST .* i);
        end
%         test = sum_ST_21V_exp;
%         u_AL_25V_exp(i,j) = u_AL_25V_exp(i,j) + sum_AL_25V_exp;
        u_AL_25V_any(i,j) = u_AL_25V_any(i,j) + sum_AL_25V_any;
%         u_AL_28V_exp(i,j) = u_AL_28V_exp(i,j) + sum_AL_28V_exp;
        u_AL_28V_any(i,j) = u_AL_28V_any(i,j) + sum_AL_28V_any;
%         u_BR_26V_exp(i,j) = u_BR_26V_exp(i,j) + sum_BR_26V_exp;
        u_BR_26V_any(i,j) = u_BR_26V_any(i,j) + sum_BR_26V_any;
%         u_BR_29V_exp(i,j) = u_BR_29V_exp(i,j) + sum_BR_29V_exp;
        u_BR_29V_any(i,j) = u_BR_29V_any(i,j) + sum_BR_29V_any;
%         u_ST_21V_exp(i,j) = u_ST_21V_exp(i,j) + sum_ST_21V_exp;
        u_ST_21V_any(i,j) = u_ST_21V_any(i,j) + sum_ST_21V_any;
    end
end

figure
AL28Plotany = plot(time,u_AL_28V_any,'--',"LineWidth",1.5)
hold on
AL28Plotexp = plot(Aluminum_28V(:,1),Aluminum_28V(:,3:end),"LineWidth",1.5)
xlim([0,2600])
title("Aluminum at 28V")
xlabel("Time [s]")
ylabel("Temperature [C]")
legend([AL28Plotany(1) AL28Plotexp(1)],{'Analytical Data','Experimental Data'},"Location","Best")
grid on

figure
AL25VPlotany = plot(time,u_AL_25V_any,'--',"LineWidth",1.5)
hold on
AL25VPlotexp = plot(Aluminum_25V(:,1),Aluminum_25V(:,3:end),"LineWidth",1.5)
xlim([0,2600])
title("Aluminum at 25V")
xlabel("Time [s]")
ylabel("Temperature [C]")
legend([AL25VPlotany(1) AL25VPlotexp(1)],{'Analytical Data','Experimental Data'},"Location","Best")
grid on 

figure
BR26VPlotany = plot(time,u_BR_26V_any,'--',"LineWidth",1.5)
hold on
BR26VPlotexp = plot(Brass_26V(:,1),Brass_26V(:,3:end),"LineWidth",1.5)
xlim([0,5000])
title("Brass at 26V")
xlabel("Time [s]")
ylabel("Temperature [C]")
legend([BR26VPlotany(1) BR26VPlotexp(1)],{'Analytical Data','Experimental Data'},"Location","Best")
grid on 

figure
BR29VPlotany = plot(time,u_BR_29V_any,'--',"LineWidth",1.5)
hold on
BR29VPlotexp = plot(Brass_29V(:,1),Brass_29V(:,3:end),"LineWidth",1.5)
xlim([0,5000])
title("Brass at 29V")
xlabel("Time [s]")
ylabel("Temperature [C]")
legend([BR29VPlotany(1) BR29VPlotexp(1)],{'Analytical Data','Experimental Data'},"Location","Best")
grid on 

figure
ST21VPlotany = plot(time,u_ST_21V_any,'--',"LineWidth",1.5)
hold on
ST21VPlotexp = plot(Steel_21V(:,1),Steel_21V(:,3:end),"LineWidth",1.5)
xlim([0,12000])
title("Steel at 21V")
xlabel("Time [s]")
ylabel("Temperature [C]")
legend([ST21VPlotany(1) ST21VPlotexp(1)],{'Analytical Data','Experimental Data'},"Location","Best")
grid on 

%% Question 5
AL25VSS = Aluminum_25V(end,end) * 0.9999;
AL28VSS = Aluminum_28V(end,end) * 0.9999;
BR26VSS = Brass_26V(end,end) * 0.9999;
BR29VSS = Brass_29V(end,end) * 0.9999;

indAL25V = find(min(abs(Aluminum_25V(:,end) - AL25VSS)) == abs(Aluminum_25V(:,end) - AL25VSS));
timeAL25V = min(Aluminum_25V(indAL25V,1));

indAL28V = find(min(abs(Aluminum_28V(:,end) - AL28VSS)) == abs(Aluminum_28V(:,end) - AL28VSS));
timeAL28V = min(Aluminum_28V(indAL28V,1));

indBR26V = find(min(abs(Brass_26V(:,end) - BR26VSS)) == abs(Brass_26V(:,end) - BR26VSS));
timeBR26V = min(Brass_26V(indBR26V,1));

indBR29V = find(min(abs(Brass_29V(:,end) - BR29VSS)) == abs(Brass_29V(:,end) - BR29VSS));
timeBR29V = min(Brass_29V(indBR29V,1));

Length = 0.1905; % [m]

AL25VQ5 = alpha_AL * timeAL25V / Length^2;
AL28VQ5 = alpha_AL * timeAL28V / Length^2;
BR26VQ5 = alpha_BR * timeBR26V / Length^2;
BR29VQ5 = alpha_BR * timeBR29V / Length^2;

%% plots for question 1
% 
% plot(thermocL,thermoc_AL_25V,'o','markerfacecolor','k')
% hold on
% plot([0,0.01,.15],[0,0.01,.15]*polyAL_25V(1,1) + AL_25V_T0,'LineWidth',1.5)
% plot([0,0.01,.15],[0,0.01,.15]*H_any(1,1) + AL_25V_T0,'LineWidth',1.5)
% xlabel("Distance from Chiller [cm]")
% xticklabels({'0','5','10','15'})
% ylabel("Temperature of Rod [C]")
% legend("Thermocouple Temperature","Experimental H Best Fit","Analytical H Best Fit","Location","NorthWest")
% title("Temperature vs. Distance from Chiller for Aluminum at 25V")
% grid on
% hold off
% 
% figure
% plot(thermocL,thermoc_AL_28V,'o','markerfacecolor','k')
% hold on
% plot([0,0.01,.15],[0,0.01,.15]*polyAL_28V(1,1) + AL_28V_T0,'LineWidth',1.5)
% plot([0,0.01,.15],[0,0.01,.15]*H_any(1,2) + AL_28V_T0,'LineWidth',1.5)
% xlabel("Distance from Chiller [cm]")
% xticklabels({'0','5','10','15'})
% ylabel("Temperature of Rod [C]")
% legend("Thermocouple Temperature","Experimental H Best Fit","Analytical H Best Fit","Location","NorthWest")
% title("Temperature vs. Distance from Chiller for Aluminum at 28V")
% grid on
% hold off
% 
% figure
% plot(thermocL,thermoc_BR_26V,'o','markerfacecolor','k')
% hold on
% plot([0,0.01,.15],[0,0.01,.15]*polyBR_26V(1,1) + BR_26V_T0,'LineWidth',1.5)
% plot([0,0.01,.15],[0,0.01,.15]*H_any(1,3) + BR_26V_T0,'LineWidth',1.5)
% xlabel("Distance from Chiller [cm]")
% xticklabels({'0','5','10','15'})
% ylabel("Temperature of Rod [C]")
% legend("Thermocouple Temperature","Experimental H Best Fit","Analytical H Best Fit","Location","NorthWest")
% title("Temperature vs. Distance from Chiller for Brass at 26V")
% grid on
% hold off
% 
% figure
% plot(thermocL,thermoc_BR_29V,'o','markerfacecolor','k')
% hold on
% plot([0,0.01,.15],[0,0.01,.15]*polyBR_29V(1,1) + BR_29V_T0,'LineWidth',1.5)
% plot([0,0.01,.15],[0,0.01,.15]*H_any(1,4) + BR_29V_T0,'LineWidth',1.5)
% xlabel("Distance from Chiller [cm]")
% xticklabels({'0','5','10','15'})
% ylabel("Temperature of Rod [C]")
% legend("Thermocouple Temperature","Experimental H Best Fit","Analytical H Best Fit","Location","NorthWest")
% title("Temperature vs. Distance from Chiller for Brass at 29V")
% grid on
% hold off
% 
% figure
% plot(thermocL,thermoc_ST_21V,'o','markerfacecolor','k')
% hold on
% plot([0,0.01,.15],[0,0.01,.15]*polyST_21V(1,1) + ST_21V_T0,'LineWidth',1.5)
% plot([0,0.01,.15],[0,0.01,.15]*H_any(1,5) + ST_21V_T0,'LineWidth',1.5)
% xlabel("Distance from Chiller [cm]")
% xticklabels({'0','5','10','15'})
% ylabel("Temperature of Rod [C]")
% legend("Thermocouple Temperature","Experimental H Best Fit","Analytical H Best Fit","Location","NorthWest")
% title("Temperature vs. Distance from Chiller for Steel at 21V")
% grid on
% hold off

% plots: x position vs temperature, with 2 lines: experimental H and
% anylitical H with same T0 y intercept
