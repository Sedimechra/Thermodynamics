%% ASEN 3113 - Lab 2 - Heat Conduction Lab
% Section 011 - Group 14
% 
% Authors:
%     1. Luca Bonarrigo
%     2. Pete Dillman
%     3. Mikaela Felix
%     4. Nathaniel Shiba
%     5. Ryan Sievers
% 
% Created: 10/12/2021 
% Last edited: 10/12/2021 
%

clc; clear; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load in experimental data
filenames = {'Aluminum_25V_240mA','Aluminum_28V_269mA',...
    'Brass_26V_245mA','Brass_29V_273mA','Steel_21V_192mA'};

data = cell(1,length(filenames));

for i=1:length(filenames)
   data{i} = load_data(filenames{i}); 
end

%% QUESTION 1
% Assuming X0 lies 1 3/8 in to the left of the first thermocouple, 
% determine the temperature at the cold end of the rod (T0). This can be
% done using the same linear extrapolation that was performed in the 
% prelab. Also, determine H_exp and plot the experimental and analytical 
% steady states. Create a plot with the experimental and analytical 
% slopes as well as the experimental steady state temperature distribution.
% Additionally, create a table containing T0 in °C, H_exp in °C/m, and
% H_analytical in°C/m for each data set. Compare these results.

% Hint: T0 is defined as the initial temperature of the rod. 
% H_analytical has an equation to determine its value.

H_exp = zeros(1,length(filenames));
T_0 = zeros(1,length(filenames));

for i=1:length(filenames)
    [H_exp(i),T_0(i)] = find_T0(data{i},filenames{i});
end

table1 = table(filenames',T_0',H_exp','VariableNames',...
    {'Trial','T_0 [°C]','H_exp [°C'});

%% QUESTION 2
% Calculating Diffusivities
% Properties read in order of p, cp, k

x = [11/8 15/8 19/8 23/8 27/8 21/8 35/8 39/8];
L = 5;
%t = 0:10:12600;

al_prop = [2810 960 130];
brass_prop = [8500 380 115];
steel_prop = [8000 500 16.2];

% alpha = k / (p * cp)
alpha_al = al_prop(3) / (al_prop(1) * al_prop(2)) * 100^2 / 2.54^2;
alpha_brass = brass_prop(3) / (brass_prop(1) * brass_prop(2)) * 100^2 / 2.54^2;
alpha_steel = steel_prop(3) / (steel_prop(1) * steel_prop(2)) * 100^2 / 2.54^2;
alpha = [alpha_al alpha_al alpha_brass alpha_brass alpha_steel];

% predefine u(t) arrays for 8 thermocouples
u_al25 = cell(1,8);
u_al28 = cell(1,8);
u_brass26 = cell(1,8);
u_brass29 = cell(1,8);
u_steel = cell(1,8);
u = {u_al25,u_al28,u_brass26,u_brass29,u_steel};

for k=1:5
    t = 0:10:data{1,k}.Time(end);
    for i=1:8
        % predefine array for each thermocouple
        u{k}{i} = zeros(1,length(t));
        
        % use find_sum to calculate u(t) using Fourier series
        for j=1:length(t)
            u{k}{1,i}(1,j) = find_sum(t(j),x(i),L,H_exp(k),alpha(k),T_0(k));
        end
    end
    % plot results
    temp_plot(data,u{k},filenames,k);
end

% diff_temp = abs(temp - u_al); % comparing experimental temp (data) against analytical temp (u)
% plot(t, diff_temp)
%% question 3

%% question 5

% u(x,t) = T_0 + Hx + sum(n=1->infinity) b_n*sin(lambda_n*x) *
% exp(-lambda_n^2 * alpha * t)

x = 4+7/8;
L = 5;
alpha = [0.381:0.05:0.581];
alpha = alpha/2.54^2;

t = 1:1:1000;
figure()
hold on
for m = 1:length(alpha)
    u = zeros(1,length(t)); 
    for k = 1:length(t)
        u(k) = find_sum(t(k),x,L,H,alpha(m),T_0);
    end
    plot(t,u);
end
xlabel('Time to heat [s]');
ylabel('Temperature [C]');
title('Convergent temperature at 4.875 inches for n=inf');
lgd = legend(string(alpha));
lgd.Title.String = 'Thermal Diffusivity';
%fprintf("Temperature at x = %.3f inches converges to %.2f deg Celsius after %.0f iterations\n",x,u(end),num);

function u = find_sum(t,x,L,H,alpha,T_0)
num = 10;
n = 1:num;
lambda_n = (2*n - 1)*pi/(2*L);
%lambda_n = pi/(2*L);
n_odd = [1:2:num-1];
b_n_odd = -8*H*L./((2*n_odd - 1)*pi).^2;

n_even = [2:2:num];
b_n_even = 8*H*L./((2*n_even - 1)*pi).^2;

b_n = [];
for i=1:length(n_even)
   b_n(end+1) = b_n_odd(i);
   b_n(end+1) = b_n_even(i);
end

SUM = zeros(1,num);
    for i=1:num-1
        A = b_n(i)*sin(lambda_n(i)*x)*exp(-lambda_n(i)^2 * alpha * t);
        %A = b_n(i)*sin(lambda_n*x)*exp(-lambda_n^2 * alpha * t);
        SUM(i+1) = SUM(i)+A; 
    end

u = SUM + T_0 + H*x;
u = u(end);
end

function T = load_data(filename)
T = readtable(filename);
T = renamevars(T,["Var1","Var2","Var3","Var4","Var5","Var6","Var7",...
    "Var8","Var9","Var10",],["Time","TC0","TC1","TC2","TC3","TC4","TC5",...
    "TC6","TC7","TC8",]);
end

function [H_exp,T_0] = find_T0(data,name)
    % locations of thermocouples 1-8, from X0 [in]
    d = [1+3/8 1+7/8 2+3/8 2+7/8 3+3/8 3+7/8 4+3/8 4+7/8];
    
    % extract T1-8 temps for very first data point (t = 0)
    T = table2array(data(end,3:end)); 
    
    figure()
    hold on
    plot(d,T);
    ylabel('Temperature [C]');
    xlabel('Linear distance [in]');
    title(strcat('Temperature vs distance along rod: ',name));

    % slope and y offset
    coef = polyfit(d,T,1);

    % slope
    H_exp = coef(1);

    % extrapolate; y intercept = T_0, since T = Hd + b and d_0 = 0
    T_0 = coef(2);

    % compare with experimental data
    fit = polyval(coef,d);
    plot(d,fit,'--r');
    legend('Experimental Data','Best Fit');
    hold off

end

function [] = temp_plot(data,u,filenames,idx)
    colors = {'#0072BD','#D95319','#EDB120','#7E2F8E',...
    '#77AC30','#4DBEEE','#A2142F','#0000FF'};

    expdata = {data{1,idx}.TC1,data{1,idx}.TC2,data{1,idx}.TC3,data{1,idx}.TC4,...
        data{1,idx}.TC5,data{1,idx}.TC6,data{1,idx}.TC7,data{1,idx}.TC8};
    t = data{1,idx}.Time;
    
    figure()
    hold on
    for i=1:8
       plot(t,u{1,i}(1,:),'Color',colors{i});
       plot(t,expdata{i},'Color',colors{i},'LineStyle','--');
    end
    
    xlabel('Time (s)');
    ylabel('Temperature (°C)');
    title(strcat('Temperature vs Time for ',filenames(idx)));
    legend('u(t) - TC1','Exp T - TC1','u(t) - TC2','Exp T - TC2',...
        'u(t) - TC3','Exp T - TC3','u(t) - TC4','Exp T - TC4',...
        'u(t) - TC5','Exp T - TC5','u(t) - TC6','Exp T - TC6',...
        'u(t) - TC7','Exp T - TC7','u(t) - TC8','Exp T - TC8');
    hold off;

end