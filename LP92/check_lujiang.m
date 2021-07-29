clear

file = load(['./LP92.dat']);
T_all = file(:,1);
e_all = file(:,5); % Strain rate (s^-1)
s_all = file(:,7); % Stress (MPa)

%flow law
n = 4.0;
Q = 152000;
R = 8.3145;
normT = 1300;

e_0 = s_all.^n.*exp(-Q./(R.*T_all));
A = e_all./e_0;
fprintf("calculated A: %d\n", mean(A));
fprintf("Lu Jiang A: %d\n", exp(-18.24));
norm_factor = (exp(-Q./(R.*normT)))./(exp(-Q./(R.*T_all)));
loglog(s_all, e_all.*norm_factor, 'ro');
hold on;
loglog(s_all, exp(-18.24).*s_all.^n.*exp(-Q./(R.*normT)), 'b-');
loglog(s_all, mean(A).*s_all.^n.*exp(-Q./(R.*normT)), 'g-');