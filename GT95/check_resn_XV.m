% check_resn.m
clear
run_id = 'GT95_V';

file = load(['./' run_id '.dat']);
T_all = file(:,1); % Temperature (K)
P_all = file(:,3)*1e9; % Pressure (Pa)
e_all = file(:,5); % Strain rate (s^-1)
s_all = file(:,7); % Stress (MPa)

% Separate samples into arrays nested in a cell array
for i=1:5
    data = file(file(:,9) == i,:);
    T{i} = data(:,1); % Temperature [K]
    dT{i} = data(:,2); 
    P{i} = data(:,3)*1e9; % Pressure [Pa]
    dP{i} = data(:,4)*1e9;
    e{i} = data(:,5); % Strain rate [s^-1]
    de{i} = data(:,6);
    s{i} = data(:,7); % Stress [MPa]
    ds{i} = data(:,8);
end

output = load(['./' run_id '.out']);
% skip first 1000 runs, then select results from every 100 runs
out = output(100:10:end,1:end);
nout = length(out);
% 4 = Stress exponent, 5 = Activation energy (J/mol), 21 = A
chi2 = out(:,3); n = out(:,4); Q = out(:,5); V = out(:,6); A = out(:,12);
% 6-10 = Inter-run bias relative to sample 1
for i=7:11
    X{i-6} = out(:,i);
end

% Print summary data
disp(['id=' run_id]);
disp(['n = ' num2str(mean(n)) ' +/- ' num2str(1*std(n))]);
disp(['Q = ' num2str(mean(Q)/1e3) ' +/- ' num2str(1*std(Q)/1e3) ' kJ/mol']);
disp(['V = ' num2str(mean(V)*1e6) ' +/- ' num2str(1*std(V)*1e6) ' m^3/mol']);
disp(['log10(A) = ' num2str(mean(log10(A))) ' +/- ' num2str(1*std(log10(A)))]);
disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(1*std(chi2))]);
for i=2:5
    disp(['X_' num2str(i) ' = ' num2str(mean(X{i})) ' +/- ' num2str(1*std(X{i}))]);
end

R = 8.3145;
nx = 20;
%stress data
xs = linspace(min(s_all),max(s_all)+150,nx);
%temperature data
xt = linspace(min(T_all)-100,max(T_all)+300,nx);
%pressure data
xP = linspace(min(P_all)-0.25*1e9,max(P_all)+0.5*1e9,nx);
%normalize data
normT = 1373; % Temperature [K]
norms = 220; % Stress [MPa]
normP = 1.5*1e9; % Pressure [Pa]

figure(1);
% Strain vs. stress
for i=1:nout
    %model strain rate
    e_pred = A(i)*xs.^n(i).*exp(-(Q(i)+normP*V(i))/(normT*R));
    loglog(xs,e_pred,'r:');
    if i==1
        hold on; box on; axis tight;
    end
end

for i=1:5
    norm_factor{i} = exp(-(mean(Q)+normP*mean(V))/(normT*R))./exp(-(mean(Q)+P{i}*mean(V))./(T{i}.*R));
    loglog(s{i},e{i}.*norm_factor{i}*exp(-mean(X{i})),'bo');
    % error bars
    [ex{i},ey{i}] = calc_err(s{i},ds{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    loglog(ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
title('Strain rate vs. Stress')
xlim([30,300])
xticks([0, 100, 200, 300])

figure(2); hold off;
% Strain vs. temp
for i=1:nout
    %model strain rate
    e_pred = A(i)*norms^n(i).*exp(-(Q(i)+normP*V(i))./(xt.*R));
    semilogy(1e4./xt,e_pred,'r:');
    if i==1
        hold on; box on; axis tight;
    end
end

for i=1:5
    e_pred_P_s{i} = norms^mean(n).*exp(-(mean(Q)+normP*mean(V))./(T{i}.*R));
    e_pred0{i} = s{i}.^mean(n).*exp(-(mean(Q)+P{i}.*mean(V))./(T{i}.*R));
    semilogy(1e4./T{i},e{i}.*e_pred_P_s{i}./e_pred0{i}.*exp(-mean(X{i})),'bo');
    % error bars
    [ex{i},ey{i}] = calc_err(T{i},dT{i},e{i}.*e_pred_P_s{i}./e_pred0{i},(de{i}./e{i}).*e{i}.*e_pred_P_s{i}./e_pred0{i});
    semilogy(1e4./ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end
ylabel('Strain rate [s^{-1}]')
xlabel('10^4/T [K^{-1}]')
title('Strain rate vs. Temperature')
xlim([7,8.5])

figure(3); hold off;
% Strain vs. pressure
for i=1:nout
    %model strain rate
    e_pred = A(i)*norms^n(i).*exp(-(Q(i)+xP.*V(i))./(normT*R));
    semilogy(xP./1e9,e_pred,'r:');
    if i==1
        hold on; box on; axis tight;
    end
end
for i=1:5
    e_pred_T_s{i} = norms^mean(n).*exp(-(mean(Q)+P{i}.*mean(V))./(normT*R));
    e_pred0{i} = s{i}.^mean(n).*exp(-(mean(Q)+P{i}.*mean(V))./(T{i}.*R));
    semilogy(P{i}./1e9,e{i}.*e_pred_T_s{i}./e_pred0{i}.*exp(-mean(X{i})),'bo');
    % error bars
    [ex{i},ey{i}] = calc_err(P{i},dP{i},e{i}.*e_pred_T_s{i}./e_pred0{i},(de{i}./e{i}).*e{i}.*e_pred_T_s{i}./e_pred0{i});
    semilogy(ex{i}./1e9,ey{i}.*exp(-mean(X{i})),'b-');
end
ylabel('Strain rate [s^{-1}]')
xlabel('Pressure [GPa]')
title('Strain rate. vs. Pressure')

figure(4); hold off;
subplot(1,2,1);
for i=1:5
    norm_factor{i} = exp(-(mean(Q)+normP*mean(V))/(normT*R))./exp(-(mean(Q)+P{i}*mean(V))./(T{i}.*R));
    loglog(s{i},e{i}.*norm_factor{i},'-o');
    hold on;
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([40,250])
title("Without inter-run bias")
subplot(1,2,2);
for i=1:5
    norm_factor{i} = exp(-(mean(Q)+normP*mean(V))/(normT*R))./exp(-(mean(Q)+P{i}*mean(V))./(T{i}.*R));
    loglog(s{i},e{i}.*norm_factor{i}*exp(-mean(X{i})),'-o');
    hold on;
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([40,250])
title("With inter-run bias")

figure(5); hold off;
for i=1:5
    plot(i,mean(X{i}),'bo');
    errorbar(i,mean(X{i}),std((X{i})),'b');
    if i==1
        ylabel('X_m')
        xlabel('m')
        xlim([1,6])
        hold on;
    end
end

 % Histograms
figure(6);
subplot(3,2,1);
hist(n,20); xlabel('n');
subplot(3,2,2);
hist(Q/1e3,20); xlabel('Q [kJ/mol]');
subplot(3,2,3);
hist(log10(A),20); xlabel('log_{10}(A)');
subplot(3,2,4);
hist(V*1e6,20); xlabel('V [m^3/mol]');
subplot(3,2,5);
hist(chi2,20); xlabel('\chi^2');
%}