% check_resp.m
clear
run_id = 'LP92';

file = load(['./' run_id '.dat']);
T_all = file(:,1);
P_all = file(:,3);
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

output = load(['./' run_id '_pXC.out']);
out = output(1000:2:end,1:end);
nout = length(out);
% 4 = Water exponent, 5 = Stress exponent, 6 = Activation energy (J/mol), 21 = A
chi2 = out(:,3); r_dif = out(:,4); m = out(:,5); Q_dif = out(:,6);
r_dis = out(:,7); n = out(:,8); Q_dis = out(:,9);
A_dif = out(:,end-1); A_dis = out(:,end);
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

% set constants
    R = 8.3145;
    nx = 50;
    %stress data
    xs = logspace(log10(min(s_all))-2,log10(max(s_all))+3,nx);
    %temperature data
    xt = linspace(min(T_all)-50,max(T_all)+50,nx);
    %COH data
    xCOH = linspace(min(COH_all)-1e-4,max(COH_all)+5e-3,nx);
    %grain size data
    xd = logspace(log10(min(d_all))-1,log10(max(d_all))+2,nx);
    %normalize data
    normT = 1273; % Temperature [K]
    norms = 300; % Stress [MPa]
    normd = 100; % Grain size [um]
    normCOH = 0.005; % Water content (ppm H/Si)

figure(1); hold off;
% Strain vs. stress
for i=1:nout
    %model strain rate
    e_dis = A_dis(i)*xs.^n(i).*normCOH^r_dis(i).*exp(-Q_dis(i)/(normT*R));
    e_dif = A_dif(i)*xs.*normd^-m(i).*normCOH^r_dif(i).*exp(-Q_dif(i)/(normT*R));
    loglog(xs,e_dif + e_dis,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('strain rate [s^{-1}]')
        xlabel('stress [MPa]')
        title('LP92')
        xlim([10^(1),10^(3.5)])
        ylim([1e-9,1e0])
        set(gca, 'LineWidth', 1)
        set(gca,'FontSize', 20)
        text(10^(1.05), 10^(-1.9), ['T = ' num2str(normT) ' K' newline 'd = ' num2str(normd) ' \mu m' newline 'C_{OH} = ' num2str(normCOH) ' H/10^6Si'], 'FontSize', 20)
    end
end

for i=1:14
    e_pred_T_d_COH{i} = mean(log10(A_dis))*s{i}.^mean(n).*normCOH.^mean(r_dis).*exp(-mean(Q_dis)./(normT.*R)) + mean(log10(A_dif))*s{i}.*normd^-mean(m).*normCOH^mean(r_dif).*exp(-mean(Q_dif)/(normT*R));
    e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*COH{i}.^mean(r_dis).*exp(-mean(Q_dis)./(T{i}.*R)) + mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*COH{i}.^mean(r_dif).*exp(-mean(Q_dif)./(T{i}.*R));
    loglog(s{i},e{i}.*(e_pred_T_d_COH{i}./e_pred0{i}).*exp(-mean(X{i})),'bo');
    % errorbars
    [ex, ey] = calc_err(s{i}, ds{i}, e{i}.*(e_pred_T_d_COH{i}./e_pred0{i}).*exp(-mean(X{i})), de{i}.*(e_pred_T_d_COH{i}./e_pred0{i}).*exp(-mean(X{i})));
    loglog(ex, ey, 'b-')
end
loglog(xs, 1.2e-8.*xs.^4.*exp(-152000/(normT*R)),'c', 'Linewidth', 1.5);
loglog(xs, 10^mean(log10(A_dis)).*xs.^mean(n).*normCOH^mean(r_dis).*exp(-mean(Q_dis)./(normT.*R)), 'k--','Linewidth', 1.5);
loglog(xs, 10^mean(log10(A_dif)).*xs.*normd^-mean(m).*normCOH^mean(r_dif).*exp(-mean(Q_dif)./(normT.*R)), 'k-','Linewidth', 1.5);
fprintf("for e = 1e-15, T = %d, K\nCalculated: %d\nLP: %d\n", normT, (1e-15./(10^mean(log10(A_dis))*exp(-mean(Q_dis)./(R.*normT)))).^(1/mean(n)), (1e-15./(1.2e-8.*exp(-152000./(R.*normT)))).^(1/4))

figure(2); hold off;
% Strain vs. temp
for i=1:nout
    %model strain rate
    e_dis = A_dis(i)*norms.^n(i).*normCOH^r_dis(i).*exp(-Q_dis(i)./(xt.*R));
    e_dif = A_dif(i)*norms.*normd^-m(i).*normCOH^r_dif(i).*exp(-Q_dif(i)./(xt.*R));
    semilogy(1e4./xt, e_dis + e_dif,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('strain rate [s^{-1}]')
        xlabel('10^4/T [K^{-1}]')
        ylim([1e-9 1e0])
        set(gca, 'LineWidth', 1)
        set(gca,'FontSize',20)
        text(7.45, 10^(-1.8), ['\sigma = ' num2str(norms) ' MPa' newline 'd = ' num2str(normd) ' \mu m' newline 'C_{OH} = ' num2str(normCOH) ' H/10^6Si'], 'FontSize', 20)
    end
end

for i=1:14
    e_pred_s_d_COH{i} = mean(log10(A_dis))*norms.^mean(n).*normCOH.^mean(r_dis).*exp(-mean(Q_dis)./(T{i}.*R)) + mean(log10(A_dif))*norms.*normd^-mean(m).*normCOH^mean(r_dif).*exp(-mean(Q_dif)./(T{i}.*R));
    e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*COH{i}.^mean(r_dis).*exp(-mean(Q_dis)./(T{i}.*R)) + mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*COH{i}.^mean(r_dif).*exp(-mean(Q_dif)./(T{i}.*R));
    semilogy(1e4./T{i},e{i}.*(e_pred_s_d_COH{i}./e_pred0{i}).*exp(-mean(X{i})),'bo');
    % errorbars
    [ex, ey] = calc_err(T{i}, dT{i}, e{i}.*(e_pred_s_d_COH{i}./e_pred0{i}).*exp(-mean(X{i})), de{i}.*(e_pred_s_d_COH{i}./e_pred0{i}).*exp(-mean(X{i})));
    semilogy(1e4./ex, ey, 'b-')
end
semilogy(1e4./xt, 1.2e-8*norms^4.*exp(-152000./(xt.*R)),'c', 'Linewidth', 1.5);
semilogy(1e4./xt, 10^mean(log10(A_dis))*norms^mean(n).*normCOH^mean(r_dis).*exp(-mean(Q_dis)./(xt.*R)), 'k--','Linewidth', 1.5);
semilogy(1e4./xt, 10^mean(log10(A_dif))*norms.*normd^-mean(m).*normCOH^mean(r_dif).*exp(-mean(Q_dif)./(xt.*R)), 'k','Linewidth', 1.5);

figure(3); hold off;
% Strain vs. grain size
for i=1:nout
    %model strain rate
    e_dis = A_dis(i)*norms.^n(i).*normCOH^r_dis(i).*exp(-Q_dis(i)./(normT.*R));
    e_dif = A_dif(i)*norms.*xd.^-m(i).*normCOH^r_dif(i).*exp(-Q_dif(i)./(normT.*R));
    loglog(xd, e_dif + e_dis,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('strain rate [s^{-1}]')
        xlabel('grain size [\mu m]')
        ylim([1e-6 1e0])
        xlim([10^(0.5) 1e3])
        set(gca, 'LineWidth', 1)
        set(gca,'FontSize',20)
        text(10^(1.4), 10^(-1.2), ['T = ' num2str(normT) ' K' newline '\sigma = ' num2str(norms) ' MPa' newline 'C_{OH} = ' num2str(normCOH) ' H/10^6Si'], 'FontSize', 20)
    end
end

for i=1:14
    e_pred_s_T_COH{i} = mean(log10(A_dis))*norms.^mean(n).*normCOH.^mean(r_dis).*exp(-mean(Q_dis)./(normT.*R)) + mean(log10(A_dif))*norms.*d{i}.^-mean(m).*normCOH^mean(r_dif).*exp(-mean(Q_dif)./(normT.*R));
    e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*COH{i}.^mean(r_dis).*exp(-mean(Q_dis)./(T{i}.*R)) + mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*COH{i}.^mean(r_dif).*exp(-mean(Q_dif)./(T{i}.*R));
    loglog(d{i},e{i}.*(e_pred_s_T_COH{i}./e_pred0{i}).*exp(-mean(X{i})),'bo');
    % errorbars
    [ex, ey] = calc_err(d{i}, dd{i}, e{i}.*(e_pred_s_T_COH{i}./e_pred0{i}).*exp(-mean(X{i})), de{i}.*(e_pred_s_T_COH{i}./e_pred0{i}).*exp(-mean(X{i})));
    loglog(ex, ey, 'b-')
end
loglog(xd, ones(size(xd)).*1.2e-8*norms^4.*exp(-152000./(normT.*R)),'c', 'Linewidth', 1.5);
loglog(xd, ones(size(xd)).*10^mean(log10(A_dis))*norms^mean(n).*normCOH^mean(r_dis)*exp(-mean(Q_dis)./(normT*R)), 'k--','Linewidth', 1.5);
loglog(xd, ones(size(xd)).*10^mean(log10(A_dif))*norms.*xd.^-mean(m).*normCOH^mean(r_dif)*exp(-mean(Q_dif)./(normT*R)), 'k','Linewidth', 1.5);

%{
figure(4); hold off;
% Strain vs. water content
for i=1:nout
    %model strain rate
    e_dis = A_dis(i)*norms.^n(i).*xCOH.^r_dis(i).*exp(-Q_dis(i)./(normT.*R));
    e_dif = A_dif(i)*norms.*normd.^-m(i).*xCOH.^r_dif(i).*exp(-Q_dif(i)./(normT.*R));
    loglog(xCOH./1e-6,e_dis+e_dif,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('Strain rate [s^{-1}]')
        xlabel('C_{OH} [ppm H/Si]')
        title('Strain rate vs. Water content')
    end
end

for i=1:14
    e_pred_s_T_d{i} = mean(log10(A_dis))*norms.^mean(n).*COH{i}.^mean(r_dis).*exp(-mean(Q_dis)./(normT.*R)) + mean(log10(A_dif))*norms.*normd.^-mean(m).*COH{i}.^mean(r_dif).*exp(-mean(Q_dif)./(normT.*R));
    e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*COH{i}.^mean(r_dis).*exp(-mean(Q_dis)./(T{i}.*R)) + mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*COH{i}.^mean(r_dif).*exp(-mean(Q_dif)./(T{i}.*R));
    loglog(COH{i}./1e-6,e{i}.*(e_pred_s_T_d{i}./e_pred0{i}).*exp(-mean(X{i})),'bo');
end
%}
% Inter-run bias
%{
figure(5); hold off;
for i=1:14
    plot(i,mean(X{i}),'bo');
    errorbar(i,mean(X{i}),std((X{i})),'b');
    if i==1
        ylabel('X_m')
        xlabel('m')
        xlim([1,15])
        hold on;
    end
end
%}

% Histograms
figure(6);
subplot(3,3,1);
hist(n,20); xlabel('n');
subplot(3,3,2);
hist(r_dis,20); xlabel('r_{dis}');
subplot(3,3,3);
hist(Q_dis/1e3,20); xlabel('Q_{dis} [kJ/mol]');
subplot(3,3,4);
hist(log10(A_dis),20); xlabel('log_{10}(A_{dis})');
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