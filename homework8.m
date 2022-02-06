clc, clear ALL;
wavelengths = [0.01:.01:1000]; % [micro meters]
C1 = 3.74177 * 10^8; % [W * micro m^4 /m^2]
C2 = 1.43878 * 10^4; % [micro meter  * K]
T = 5780; % [K]
Emissiv = C1 ./ (wavelengths.^5 .* (exp(C2 ./ (wavelengths .* T)) - 1));

plot(log10(wavelengths),log10(Emissiv),"LineWidth",1.5)
% xticks([-1:0.5:1])
xticklabels(["10^{-2}"," ","10^{-1}"," ","10^0"," ","10^1"," ","10^2"," ","10^3"])
xlim([log10(0.01),log10(1000)])
ylim([-6,8.5])
xline(log10(0.7),"--","Linewidth",1.5)
xline(log10(0.4),"--","Linewidth",1.5)
yticklabels(["10^{-6}","10^{-4}","10^{-2}","10^{0}","10^{2}","10^{4}","10^{6}","10^{8}"])
grid on
ylabel("Emissive Power [W/m^2\mum]")
xlabel("Wavelength [\mum]")
title("Wavelength vs Emissive Power at T = 5780 K")