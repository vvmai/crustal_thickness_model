% check_resp.m
clear
run_id = 'GT95_d';

% load data file
    file = load(['./' run_id '_fixed.dat']);
    T_all = file(:,1); % Temperature (K)
    e_all = file(:,5); % Strain rate (s^-1)
    s_all = file(:,7); % Stress (MPa)
    d_all = file(:,9); % grain size (um)

% Separate samples into arrays nested in a cell array
for i=1:5
    data = file(file(:,11) == i,:);
    T{i} = data(:,1); % Temperature [K]
    dT{i} = data(:,2); 
    e{i} = data(:,5); % Strain rate [s^-1]
    de{i} = data(:,6);
    s{i} = data(:,7); % Stress [MPa]
    ds{i} = data(:,8);
    d{i} = data(:,9); % Grain size [um]
    dd{i} = data(:,10);
end
output = load(['./' run_id '_fixed.out']);
% skip first 1000 runs, then select results from every 100 runs
out = output(100:10:end,1:end);
nout = length(out);
% 4 = Stress exponent, 5 = Activation energy (J/mol), 6 = A
chi2 = out(:,3); m = out(:,4); Q_dif = out(:,5); n = out(:,6); Q_dis = out(:,7);
A_dif = out(:,end-1); A_dis = out(:,end);
% 10-23 = Inter-run bias relative to sample 1
for i=8:12
    X{i-7} = out(:,i);
end

% Print summary data
    disp(['id=' run_id]);
    disp(['n = ' num2str(mean(n)) ' +/- ' num2str(1*std(n))]);
    disp(['Q_dis = ' num2str(mean(Q_dis)/1e3) ' +/- ' num2str(1*std(Q_dis)/1e3) ' kJ/mol']);
    disp(['log10(A_dis) = ' num2str(mean(log10(A_dis))) ' +/- ' num2str(1*std(log10(A_dis)))]);
    disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(1*std(chi2))]);

    disp(['m = ' num2str(mean(m)) ' +/- ' num2str(1*std(m))]);
    disp(['Q_dif = ' num2str(mean(Q_dif)/1e3) ' +/- ' num2str(1*std(Q_dif)/1e3) ' kJ/mol']);
    disp(['log10(A_dif) = ' num2str(mean(log10(A_dif))) ' +/- ' num2str(1*std(log10(A_dif)))]);
    for i=2:5
        disp(['X_' num2str(i) ' = ' num2str(mean(X{i})) ' +/- ' num2str(1*std(X{i}))]);
    end

% set constants
    R = 8.3145;
    nx = 20;
    %stress data
    xs = logspace(log10(min(s_all))-1,log10(max(s_all))+2,nx);
    %temperature data
    xt = linspace(min(T_all)-100,max(T_all)+100,nx);
    %grain size data
    xd = logspace(log10(min(d_all))-1,log10(max(d_all))+2,nx);
    %normalize data
    normT = 1373; % Temperature [K]
    norms = 250; % Stress [MPa]
    normd = 100; % Grain size [um]

mat = [log10(A_dis) log10(A_dif) n m Q_dis Q_dif];
[cc] = corrcoef(mat)


% Strain vs. stress
figure(1); hold off;
for i=1:nout
    %model strain rate
    e_dis = A_dis(i).*xs.^n(i).*exp(-Q_dis(i)./(normT*R));
    e_dif = A_dif(i).*xs.*normd^-m(i).*exp(-Q_dif(i)./(normT*R));
    loglog(xs,e_dis+e_dif,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('strain rate [s^{-1}]')
        xlabel('stress [MPa]')
        title('GT95')
        set(gca, 'LineWidth', 1)
        set(gca,'FontSize',20)
        xlim([10^(1) 10^(3)])
        text(12, 10^(-0.9), ['T = ' num2str(normT) ' K' newline 'd = ' num2str(normd) ' \mu m'], 'FontSize', 20)
    end
end

for i=1:5
    e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*exp(-mean(Q_dis)./(T{i}.*R)) + ...
    mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*exp(-mean(Q_dif)./(T{i}.*R));
    e_pred_T_d_P{i} = mean(log10(A_dis))*s{i}.^mean(n).*exp(-mean(Q_dis)./(normT.*R)) + ...
    mean(log10(A_dif))*s{i}.*normd^-mean(m).*exp(-mean(Q_dif)/(normT*R));
    norm_factor{i} = e_pred_T_d_P{i}./e_pred0{i};
    loglog(s{i},e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
    %loglog(s{i},e{i},'g*');
    % error bars
    [ex{i},ey{i}] = calc_err(s{i},ds{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    loglog(ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end

% GT95's flow law -- difference due to activation volume?
loglog(xs, exp(-7.58).*xs.^4.*exp(-223000/(R*normT)),'c', 'Linewidth', 1.5)
mean_dis = 10^mean(log10(A_dis))*xs.^mean(n).*exp(-(mean(Q_dis))/(normT*R));
mean_dif = 10^mean(log10(A_dif))*xs.*normd^-mean(m).*exp(-(mean(Q_dif))/(normT*R));
loglog(xs,mean_dif,'k-','LineWidth', 1.5);
loglog(xs,mean_dis,'k--','LineWidth',1.5);

figure(2); hold off;
% Strain vs. temp
for i=1:nout
    %model strain rate
    e_dis = A_dis(i)*norms.^n(i).*exp(-Q_dis(i)./(xt.*R));
    e_dif = A_dif(i)*norms.*normd^-m(i).*exp(-Q_dif(i)./(xt.*R));
    semilogy(1e4./xt,e_dis + e_dif,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('strain rate [s^{-1}]')
        xlabel('10^4/T [K^{-1}]')
        xlim([7.5 8.5])
        ylim([1e-6 1e-2])
        set(gca, 'LineWidth', 1)
        set(gca,'FontSize',20)
        text(8.08, 10^(-2.5),['\sigma = ' num2str(norms) ' MPa' newline 'd = ' num2str(normd) ' \mu m'], 'FontSize', 20)
    end
end

for i=1:5
    e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*exp(-mean(Q_dis)./(T{i}.*R)) + ...
    mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*exp(-mean(Q_dif)./(T{i}.*R));
    e_pred_s_d_P{i} = mean(log10(A_dis))*norms.^mean(n).*exp(-mean(Q_dis)./(T{i}.*R)) + ...
    mean(log10(A_dif))*norms.*normd^-mean(m).*exp(-mean(Q_dif)./(T{i}.*R));
    norm_factor{i} = e_pred_s_d_P{i}./e_pred0{i};
    semilogy(1e4./T{i},e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
    %semilogy(1e4./T{i},e{i},'g*');
    % error bars
    [ex{i},ey{i}] = calc_err(T{i},dT{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    semilogy(1e4./ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end
loglog(1e4./xt, exp(-7.58)*norms^4.*exp(-223000./(R.*xt)),'c', 'Linewidth', 1.5)
mean_dis = 10^mean(log10(A_dis))*norms.^mean(n).*exp(-(mean(Q_dis))./(xt.*R));
mean_dif = 10^mean(log10(A_dif))*norms.*normd^-mean(m).*exp(-(mean(Q_dif))./(xt.*R));
loglog(1e4./xt,mean_dif,'k-','LineWidth', 1.5);
loglog(1e4./xt,mean_dis,'k--','LineWidth',1.5);

figure(3); hold off;
% Strain vs. grain size
for i=1:nout
    %model strain rate
    e_dis = A_dis(i)*norms.^n(i).*exp(-Q_dis(i)./(normT.*R));
    e_dif = A_dif(i)*norms*xd.^-m(i).*exp(-Q_dif(i)./(normT.*R));
    loglog(xd,e_dis + e_dif,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('strain rate [s^{-1}]')
        xlabel('grain size [\mum]')
        xlim([10^(1.6) 10^(2.3)])
        ylim([10^(-5) 10^(0)])
        set(gca, 'LineWidth', 1)
        set(gca,'FontSize',20)
        text(10^(2), 10^(-0.6), ['\sigma = ' num2str(norms) ' MPa' newline 'T = ' num2str(normT) ' K'], 'FontSize', 20)
    end
end

for i=1:5
    e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*exp(-mean(Q_dis)./(T{i}.*R)) + ...
    mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*exp(-mean(Q_dif)./(T{i}.*R));
    e_pred_s_T_P{i} = mean(log10(A_dis))*norms.^mean(n).*exp(-mean(Q_dis)./(normT.*R)) + ...
    mean(log10(A_dif))*norms.*d{i}.^-mean(m).*exp(-mean(Q_dif)./(normT.*R));
    norm_factor{i} = e_pred_s_T_P{i}./e_pred0{i};
    loglog(d{i}+i,e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
    %loglog(d{i},e{i},'g*');
    % error bars
    [ex{i},ey{i}] = calc_err(d{i},dd{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    loglog(ex{i}+i,ey{i}.*exp(-mean(X{i})),'b-');
end
loglog(xd, ones(size(xd))*exp(-7.58)*norms^4*exp(-223000/(R*normT)),'c', 'Linewidth', 1.5)
mean_dis = 10^mean(log10(A_dis))*norms.^mean(n).*exp(-(mean(Q_dis))./(normT.*R));
mean_dif = 10^mean(log10(A_dif))*norms.*xd.^-mean(m).*exp(-(mean(Q_dif))./(normT.*R));
loglog(xd,mean_dif,'k-','LineWidth', 1.5);
loglog(xd,ones(size(xd)).*mean_dis,'k--','LineWidth',1.5);


% Histograms
 figure(6);
 subplot(3,3,1);
 hist(n,20); xlabel('n');
 subplot(3,3,2);
 hist(m,20); xlabel('m');
 subplot(3,3,3);
 hist(Q_dis/1e3,20); xlabel('Q_{dis} [kJ/mol]');
 subplot(3,3,4);
 hist(Q_dif/1e3,20); xlabel('Q_{dif} [kJ/mol]');
 subplot(3,3,5);
 hist(log10(A_dis),20); xlabel('log_{10}(A_{dis})');
 subplot(3,3,6);
 hist(log10(A_dif),20); xlabel('log_{10}(A_{dif})');
 subplot(3,3,9);
 hist(chi2,20); xlabel('\chi^2');
%}

%{
figure(5); hold off;
subplot(2,1,1);
for i=1:5
    e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*exp(-mean(Q_dis)./(T{i}.*R)) + ...
    mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*exp(-mean(Q_dif)./(T{i}.*R));
    e_pred_T_d_P{i} = mean(log10(A_dis))*s{i}.^mean(n).*exp(-mean(Q_dis)./(normT.*R)) + ...
    mean(log10(A_dif))*s{i}.*normd^-mean(m).*exp(-mean(Q_dif)/(normT*R));
    norm_factor{i} = e_pred_T_d_P{i}./e_pred0{i};
    loglog(s{i},e{i}.*norm_factor{i},'-o');
    hold on;
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([0,250])
title("Without inter-run bias")
set(gca, 'LineWidth', 1)
set(gca,'FontSize',20)

subplot(2,1,2);
for i=1:5
    e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*exp(-mean(Q_dis)./(T{i}.*R)) + ...
    mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*exp(-mean(Q_dif)./(T{i}.*R));
    e_pred_T_d_P{i} = mean(log10(A_dis))*s{i}.^mean(n).*exp(-mean(Q_dis)./(normT.*R)) + ...
    mean(log10(A_dif))*s{i}.*normd^-mean(m).*exp(-mean(Q_dif)/(normT*R));
    norm_factor{i} = e_pred_T_d_P{i}./e_pred0{i};
    loglog(s{i},e{i}.*norm_factor{i}*exp(-mean(X{i})),'-o');
    hold on;
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([0,250])
title("With inter-run bias")
set(gca, 'LineWidth', 1)
set(gca,'FontSize',20)
%}
%}