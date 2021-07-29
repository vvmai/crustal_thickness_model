clear
run_id = 'rutter04ab';

% Rutter & Brodie's parameters:
%Q = 242000;
%A = 10^-4.93;
%n = 2.97;
% Rutter & Brodie's parameters:
%Q = 220000; 
%A = 0.4;
%m = 2;

output = load(['./rutter04ab_X.out']);
out = output(100:10:end,1:end);
nout = length(out);
%1 = id, 2 = simplified chi2, 3 = original chi2, 4 = grain size exponent, 5 = Q diffusion 
%6 = stress exp., 7 = Q dislocation, 8 = A diffusion, 9 = A dislocation
m = out(:,4); Q_dif = out(:,5); n = out(:,6); Q_dis = out(:,7);
A_dif = log10(out(:,10)); A_dis = log10(out(:,11));

% Geotherm
% Depth from 0 to 40 km, in increments of 1 km
z = linspace(0,120000,100);
% Heat flow equation
k = 2.7; % Thermal conductivity [W/mK] (assume constant)
D = 10000; % Depth of upper crust producting radioactive heat [m]
T_0 = 0; % surface temperature [C]
Q_0 = 0.046; % Surface heat flow [W/m2]
Q_r = 0.025; % Nonradiogenic heat flow [W/m2]
A_0 = (Q_0 - Q_r)/D; % Surface heat production [W/m3]
T = zeros(1, length(z));
Q = zeros(1, length(z));
Q(1) = Q_r;
T(1) = T_0;
for i=2:length(z)
    % Exponentially decreasing rate of radiogenic heat production
    A(i) = A_0*exp(-z(i)/D);
    % Heat flow as a function of depth
    Q(i) = Q(1) - A(i)*D;
    % Temperature as a function of depth
    T(i) = T(1) + (Q(1)/k)*z(i) + (A_0*D^2)/(k) - (A(i).*D^2)/(k); 
end
%{
% Plot geotherm
figure(1);
set(gca,'xaxislocation','top'); hold on;
plot(T,-z/1000);
plot(400*ones(length(z)),-z/1000,'k--');
plot(200*ones(length(z)),-z/1000,'k--');
xlabel('Temperature (C)')
ylabel('Depth (km)')
%}
figure(3);
R = 8.3145;
refd = 1000; % Reference grain size [um]
refE = 1e-15; % Reference strain rate [s-1]
T = T + 273.15; % Switch to K

% RB04 flow law
s_dis = (refE./(10^-4.93*exp(-242000./(R.*T)))).^(1/2.97);
s_dif = ((refE*refd^2)./(0.4*exp(-220000./(R.*T))));
% Inversion flow law
s_dis_MCMC = (refE./(10^mean(A_dis)*exp(-mean(Q_dis)./(R.*T)))).^(1/mean(n));
s_dif_MCMC = ((refE*refd^mean(m))./(10^mean(A_dif)*exp(-mean(Q_dif)./(R.*T))));

T = T - 273.15; % Switch to C
semilogx(s_dis,-z/1000,'b');
set(gca,'xaxislocation','top'); hold on;
semilogx(s_dif,-z/1000,'r');

semilogx(s_dif_MCMC,-z/1000,'r--');
semilogx(s_dis_MCMC,-z/1000,'b--');

xlim([1e-1 1e5])
ylim([-120 0])
xlabel('Stress [MPa]')
ylabel('Depth [km]')
title(['Stress vs. Depth [e = ' num2str(refE) ' s^{-1}]'])
%}