temperatures = [18.53, 22.47, 26.87, 30.05, 35.87, 38.56, 41.50, 46.26];
thermoc = linspace(0.034925*100,0.123825*100,8);
x = polyfit(thermoc,temperatures,1);

% plot(thermoc,temperatures,'o',"LineWidth",2)
% xlim([0,14])
% ylim([0,50])
% grid on
% hold on
% [0:15].*x(1,1)+x(1,2)
% plot([0:15],([0:15].*x(1,1)+x(1,2)),"LineWidth",2)
% title("Thermocouple Location vs. Temperature")
% xlabel("Location [cm]")
% ylabel("Temperature [C]")
% legend("Experimental Values","Best Fit Line","Location","best")


T_0 = 7.949 + 273.15; % [K]
H = 309.48; % [K/m]
x = thermoc(end)/100;
L = 0.127; % [m]
alpha = 4.82 * 10^-5;
u_xt = zeros(1,11);
u_xt_2 = zeros(1,11);
u_xt(1,1) = (T_0 + H*x);
u_xt_2(1,1) = (T_0 + H*x);
t = 1; % [s]
t_2 = 1000;

% BELOW IS FOR PART 2

% for i=2:11
%     u_xt(1,i) = u_xt(1,i-1) + (-((8 * H * L * (-1)^(i-1))/(pi^2 * (2*(i-1)-1)^2)) * sin((((2*(i-1) - 1) * pi)/(2*L)) * x) * exp(-(((2*(i-1) - 1) * pi)/(2*L))^2 * alpha * t));
%     u_xt_2(1,i) = u_xt_2(1,i-1) + (-((8 * H * L * (-1)^(i-1))/(pi^2 * (2*(i-1)-1)^2)) * sin((((2*(i-1) - 1) * pi)/(2*L)) * x) * exp(-(((2*(i-1) - 1) * pi)/(2*L))^2 * alpha * t_2));
% end
% plot([0:10],u_xt,"LineWidth",2)
% title("Number of Terms in Summation vs. Convergence")
% hold on
% plot([0:10],u_xt_2,"LineWidth",2)
% grid on
% xlabel("n Number of Terms")
% ylabel("Temperature [K]")
% legend("t = 1 s","t = 1000 s","Location","best")

% BELOW IS FOR PART 3


time = [0:1000];
temps_3 = (T_0 + H * x) - (8 * H * L)/pi^2 * sin((pi/(2*L))*x) * exp(-(pi/(2*L))^2 * (2*alpha) .* time);
temps_3_2 = (T_0 + H * x) - (8 * H * L)/pi^2 * sin((pi/(2*L))*x) * exp(-(pi/(2*L))^2 * alpha .* time);
temps_3_3 = (T_0 + H * x) - (8 * H * L)/pi^2 * sin((pi/(2*L))*x) * exp(-(pi/(2*L))^2 * (.5*alpha) .* time);
plot(time,temps_3,"Linewidth",2)
grid on
hold on
plot(time,temps_3_2,"Linewidth",2)
plot(time,temps_3_3,"Linewidth",2)
xlabel("Time [s]")
ylabel("Temperature [K]")
title("Temperature vs Time for Differing Thermal Diffusivity")
legend("\alpha = 2\alpha_{Al}","\alpha = \alpha_{Al}","\alpha = 0.5\alpha_{Al}","Location","best")