% check_resn.m
clear
run_id = 'GT95';

file = load(['./' run_id '.dat']);
T_all = file(:,1);
e_all = file(:,5); % Strain rate (s^-1)
s_all = file(:,7); % Stress (MPa)

% Separate samples into arrays nested in a cell array
for i=1:5
    data = file(file(:,9) == i,:);
    T{i} = data(:,1);
    dT{i} = data(:,2);
    e{i} = data(:,5);
    de{i} = data(:,6);
    s{i} = data(:,7);
    ds{i} = data(:,8);
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
xs = linspace(min(s_all)-15,max(s_all)+150,nx);
%temperature data
xt = linspace(min(T_all)-100,max(T_all)+300,nx);
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
    norm_factor{i} = exp(-mean(Q)/(normT*R))./exp(-mean(Q)./(T{i}.*R));
    loglog(s{i},e{i}.*norm_factor{i}*exp(-mean(X{i})),'bo');
    % error bars
    [ex{i},ey{i}] = calc_err(s{i},ds{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    loglog(ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end
% inversion NOT using stress correction
loglog(xs, exp(-7.58).*xs.^4.*exp(-223000/(R*normT)),'c', 'Linewidth', 1)
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([30,300])

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
    norm_factor{i} = norms^mean(n)./s{i}.^mean(n);
    semilogy(1e4./T{i},e{i}.*norm_factor{i}*exp(-mean(X{i})),'bo');
    % error bars
    [ex{i},ey{i}] = calc_err(T{i},dT{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    semilogy(1e4./ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end
ylabel('Strain rate [s^{-1}]')
xlabel('10^4/T [K^{-1}]')
xlim([7,8.5])
%}
%{
figure(5); hold off;
subplot(2,1,1);
for i=1:5
    norm_factor{i} = exp(-mean(Q)/(normT*R))./exp(-mean(Q)./(T{i}.*R));
    loglog(s{i},e{i}.*norm_factor{i},'-o');
    hold on;
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([40,250])
title("Without inter-run bias")
subplot(2,1,2);
for i=1:5
    norm_factor{i} = exp(-mean(Q)/(normT*R))./exp(-mean(Q)./(T{i}.*R));
    loglog(s{i},e{i}.*norm_factor{i}.*exp(-mean(X{i})),'-o');
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
%}
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