% check_resn.m
clear
run_id = 'GT95';

file = load(['./' run_id '.dat']);
T_GT_all = file(:,1);
e_GT_all = file(:,5); % Strain rate (s^-1)
s_GT_all = file(:,7); % Stress (MPa)
file2 = load(['./fukuda18.dat']);
T_f_all = file2(:,1); % Temperature (K)
e_f_all = file(:,5); % Strain rate (s^-1)
s_f_all = file(:,7); % Stress (MPa)

% Separate samples into arrays nested in a cell array
for i=1:5
    data = file(file(:,9) == i,:);
    T_GT{i} = data(:,1);
    dT_GT{i} = data(:,2);
    e_GT{i} = data(:,5);
    de_GT{i} = data(:,6);
    s_GT{i} = data(:,7);
    ds_GT{i} = data(:,8);
end

for i=1:file2(end,end)
    data2 = file2(file2(:,end) == i,:);
    T_f{i} = data2(:,1); % Temperature [K]
    dT_f{i} = data2(:,2); 
    e_f{i} = data2(:,5); % Strain rate [s^-1]
    de_f{i} = data2(:,6);
    s_f{i} = data2(:,7); % Stress [MPa]
    ds_f{i} = data2(:,8);
end

output = load(['./' run_id '_X.out']);
% skip first 1000 runs, then select results from every 100 runs
out = output(100:10:end,1:end);
nout = length(out);
% 4 = Stress exponent, 5 = Activation energy (J/mol), 21 = A
chi2 = out(:,3); n = out(:,4); Q = out(:,5); A = out(:,11);
% 6-10 = Inter-run bias relative to sample 1
for i=6:10
    X{i-5} = out(:,i);
end

% Print summary data
disp(['id=' run_id]);
disp(['n = ' num2str(mean(n)) ' +/- ' num2str(1*std(n))]);
disp(['Q = ' num2str(mean(Q)/1e3) ' +/- ' num2str(1*std(Q)/1e3) ' kJ/mol']);
disp(['log10(A) = ' num2str(mean(log10(A))) ' +/- ' num2str(1*std(log10(A)))]);
disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(1*std(chi2))]);
for i=2:5
    disp(['X_' num2str(i) ' = ' num2str(mean(X{i})) ' +/- ' num2str(1*std(X{i}))]);
end

R = 8.3145;
nx = 20;
%stress data
xs = logspace(log10(min(s_GT_all))-1,log10(max(s_GT_all))+3,nx);
%temperature data
xt = linspace(min(T_GT_all)-400,max(T_GT_all)+300,nx);
%normalize data
normT = 1373;
norms = 220;

figure(1);
% Strain vs. stress
for i=1:nout
    %model strain rate
    e_pred = A(i)*xs.^n(i).*exp(-Q(i)/(normT*R));
    loglog(xs,e_pred,'r:');
    if i==1
        hold on; box on; axis tight;
    end
end

for i=1:5
    norm_factor_GT{i} = exp(-mean(Q)/(normT*R))./exp(-mean(Q)./(T_GT{i}.*R));
    % GT data
    loglog(s_GT{i},e_GT{i}.*norm_factor_GT{i}*exp(-mean(X{i})),'bo');
    % error bars
    [ex{i},ey{i}] = calc_err(s_GT{i},ds_GT{i},e_GT{i}.*norm_factor_GT{i},(de_GT{i}./e_GT{i}).*e_GT{i}.*norm_factor_GT{i});
    loglog(ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end
for i=1:file2(end:end)
    % fukuda data
    norm_factor_f{i} = exp(-mean(Q)/(normT*R))./exp(-mean(Q)./(T_f{i}.*R));
    loglog(s_f{i},e_f{i}.*norm_factor_f{i},'go');
    % error bars
    [ex{i},ey{i}] = calc_err(s_f{i},ds_f{i},e_f{i}.*norm_factor_f{i},(de_f{i}./e_f{i}).*e_f{i}.*norm_factor_f{i});
    loglog(ex{i},ey{i},'g-');
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([20,500])

figure(2);
% Strain vs. temp
for i=1:nout
    %model strain rate
    e_pred = A(i)*norms^n(i).*exp(-Q(i)./(xt.*R));
    semilogy(1e4./xt,e_pred,'r:');
    if i==1
        hold on; box on; axis tight;
    end
end

for i=1:5
    norm_factor{i} = norms^mean(n)./s_GT{i}.^mean(n);
    semilogy(1e4./T_GT{i}+(0.01*i),e_GT{i}.*norm_factor{i}*exp(-mean(X{i})),'bo');
    % error bars
    [ex{i},ey{i}] = calc_err(T_GT{i},dT_GT{i},e_GT{i}.*norm_factor{i},(de_GT{i}./e_GT{i}).*e_GT{i}.*norm_factor{i});
    semilogy(1e4./ex{i}+(0.01*i),ey{i}.*exp(-mean(X{i})),'b-');
end
for i=1:file2(end:end)
    disp("i: " + i)
    % fukuda data
    norm_factor_f{i} = norms^mean(n)./s_f{i}.^mean(n);
    semilogy(1e4./T_f{i}+(0.01*i),e_f{i}.*norm_factor_f{i},'go');
    % error bars
    [ex{i},ey{i}] = calc_err(T_f{i},dT_f{i},e_f{i}.*norm_factor_f{i},(de_f{i}./e_f{i}).*e_f{i}.*norm_factor_f{i});
    semilogy(1e4./ex{i}+(0.01*i),ey{i},'g-');
end
ylabel('Strain rate [s^{-1}]')
xlabel('10^4/T [K^{-1}]')
%}
%{
figure(5); hold off;
subplot(2,1,1);
for i=1:5
    norm_factor{i} = exp(-mean(Q)/(normT*R))./exp(-mean(Q)./(T_GT{i}.*R));
    loglog(s_GT{i},e_GT{i}.*norm_factor{i},'-o');
    hold on;
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([40,250])
title("Without inter-run bias")
subplot(2,1,2);
for i=1:5
    norm_factor{i} = exp(-mean(Q)/(normT*R))./exp(-mean(Q)./(T_GT{i}.*R));
    loglog(s_GT{i},e_GT{i}.*norm_factor{i}.*exp(-mean(X{i})),'-o');
    hold on;
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([40,250])
title("With inter-run bias")

figure(3); hold off;
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

figure(4);
subplot(2,2,1);
hist(n,20); xlabel('n');
subplot(2,2,2);
hist(Q/1e3,20); xlabel('Q [kJ/mol]');
subplot(2,2,3);
hist(log10(A),20); xlabel('log_{10}(A)');
subplot(2,2,4);
hist(chi2,20); xlabel('\chi^2');
%}