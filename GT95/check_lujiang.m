clear

file = load(['./GT95.dat']);
T_all = file(:,1);
e_all = file(:,5); % Strain rate (s^-1)
s_all = file(:,7); % Stress (MPa)
%flow law
A_0 = exp(-9.12);
n = 4.0;
Q = 223000;
R = 8.3145;
normT = 1373;

e_dot = A_0.*s_all.^n.*exp(-Q./(R.*T_all));
s_corr = 0.73*s_all;
A = e_dot./(s_corr.^n.*exp(-Q./(R.*T_all)));
fprintf("calculated A: %d\n", mean(A))
fprintf("Lu Jiang A: %d\n", exp(-7.58))