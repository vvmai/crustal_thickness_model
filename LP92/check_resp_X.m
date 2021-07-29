% check_resp.m
clear
run_id = 'LP92';

file = load(['./' run_id '.dat']);
T_all = file(:,1);
e_all = file(:,5); % Strain rate (s^-1)
s_all = file(:,7); % Stress (MPa)
d_all = file(:,9); % Grain size (um)
COH_all = file(:,11);  % Water content (ppm H/Si)

% Separate samples into arrays nested in a cell array
for i=1:14
    data = file(file(:,13) == i,:);
    T{i} = data(:,1);
    dT{i} = data(:,2);
    e{i} = data(:,5);
    de{i} = data(:,6);
    s{i} = data(:,7);
    ds{i} = data(:,8);
    d{i} = data(:,9);
    dd{i} = data(:,10);
    COH{i} = data(:,11);
    dCOH{i} = data(:,12);
end

output = load(['./' 'LP92_pX' '.out']);
% skip first 1000 runs, then select results from every 100 runs
out = output(100:10:end,1:end);
nout = length(out);
% 4 = Water exponent, 5 = Stress exponent, 6 = Activation energy (J/mol), 21 = A
chi2 = out(:,3); r_dif = out(:,4); m = out(:,5); Q_dif = out(:,6);
r_dis = out(:,7); n = out(:,8); Q_dis = out(:,9);
A_dif = out(:,24); A_dis = out(:,25);
% 10-23 = Inter-run bias relative to sample 1
for i=10:23
    X{i-9} = out(:,i);
end

% Print summary data
disp(['id=' run_id]);
disp(['n = ' num2str(mean(n)) ' +/- ' num2str(std(n))]);
disp(['m = ' num2str(mean(m)) ' +/- ' num2str(std(m))]);
disp(['r_dis = ' num2str(mean(r_dis)) ' +/- ' num2str(std(r_dis))]);
disp(['r_dif = ' num2str(mean(r_dif)) ' +/- ' num2str(std(r_dif))]);
disp(['Q_dis = ' num2str(mean(Q_dis)/1e3) ' +/- ' num2str(std(Q_dis)/1e3) ' kJ/mol']);
disp(['Q_dif = ' num2str(mean(Q_dif)/1e3) ' +/- ' num2str(std(Q_dif)/1e3) ' kJ/mol']);
disp(['log10(A_dis) = ' num2str(mean(log10(A_dis))) ' +/- ' num2str(std(log10(A_dis)))]);
disp(['log10(A_dif) = ' num2str(mean(log10(A_dif))) ' +/- ' num2str(std(log10(A_dif)))]);
disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(std(chi2))]);
for i=2:14
    disp(['X_' num2str(i) ' = ' num2str(mean(X{i})) ' +/- ' num2str(std(X{i}))]);
end

R = 8.3145;
nx = 50;
%stress data
xs = linspace(min(s_all)-50,max(s_all)+150,nx);
%temperature data
xt = linspace(min(T_all)-100,max(T_all)+300,nx);
%COH data
xCOH = linspace(min(COH_all)-4e-4,max(COH_all)+5e-3,nx);
%normalize data
normT = 1373;
norms = 220;
normd = 1.3; 
normCOH = 4000e-6;

figure(1);
% Strain vs. stress
for i=1:nout
    %model strain rate
    e_dis = A_dis(i)*xs.^n(i).*normCOH^r_dis(i).*exp(-Q_dis(i)/(normT*R));
    e_dif = A_dif(i)*xs.*normd^-m(i).*normCOH^r_dif(i).*exp(-Q_dif(i)/(normT*R));
    loglog(xs,e_dis+e_dif,'r:');
    if i==1
        hold on; box on; axis tight;
    end
end

for i=1:14
    e_pred_T_d_COH{i} = mean(log10(A_dis))*s{i}.^mean(n).*normCOH.^mean(r_dis).*exp(-mean(Q_dis)./(normT.*R)) + mean(log10(A_dif))*s{i}.*normd^-mean(m).*normCOH^mean(r_dif).*exp(-mean(Q_dif)/(normT*R));
    e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*COH{i}.^mean(r_dis).*exp(-mean(Q_dis)./(T{i}.*R)) + mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*COH{i}.^mean(r_dif).*exp(-mean(Q_dif)./(T{i}.*R));
    loglog(s{i},e{i}.*(e_pred_T_d_COH{i}./e_pred0{i}).*exp(-mean(X{i})),'bo');
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([100,700])

 % Histograms
figure(4);
subplot(3,3,1);
hist(n,20); xlabel('n');
subplot(3,3,2);
hist(r_dis,20); xlabel('r_{dis}');
subplot(3,3,3);
hist(Q_dis/1e3,20); xlabel('Q_{dis} [kJ/mol]');
subplot(3,3,4);
hist(log10(A_dis),20); xlabel('log_{10}(A_dis)');
subplot(3,3,5);
hist(m,20); xlabel('m');
subplot(3,3,6);
hist(r_dif,20); xlabel('r_{dif}');
subplot(3,3,7);
hist(Q_dif/1e3,20); xlabel('Q_{dif} [kJ/mol]');
subplot(3,3,8);
hist(log10(A_dif),20); xlabel('log_{10}(A_{dif})');
subplot(3,3,9);
hist(chi2,20); xlabel('\chi^2');
%}