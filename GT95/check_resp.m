% check_resp.m
clear
run_id = 'GT95_d';
% load data file
    file = load(['./' run_id '.dat']);
    T_all = file(:,1); % Temperature (K)
    P_all = file(:,3)*1e9; % Pressure (Pa)
    e_all = file(:,5); % Strain rate (s^-1)
    s_all = file(:,7); % Stress (MPa)
    d_all = file(:,9); % grain size (um)

% Separate samples into arrays nested in a cell array
for i=1:5
    data = file(file(:,11) == i,:);
    T{i} = data(:,1); % Temperature [K]
    dT{i} = data(:,2); 
    P{i} = data(:,3)*1e9; % Pressure [Pa]
    dP{i} = data(:,4)*1e9;
    e{i} = data(:,5); % Strain rate [s^-1]
    de{i} = data(:,6);
    s{i} = data(:,7); % Stress [MPa]
    ds{i} = data(:,8);
    d{i} = data(:,9); % Grain size [um]
    dd{i} = data(:,10);
end
output = load(['./' run_id 'C.out']);
% skip first 1000 runs, then select results from every 20 runs
out = output(100:2:end,1:end);
nout = length(out);
% 4 = Stress exponent, 5 = Activation energy (J/mol), 6 = A
chi2 = out(:,3); m = out(:,4); Q_dif = out(:,5); V_dif = out(:,6); n = out(:,7); Q_dis = out(:,8); V_dis = out(:,9);
A_dif = out(:,15); A_dis = out(:,16);
% 10-23 = Inter-run bias relative to sample 1
for i=10:14
    X{i-9} = out(:,i);
end

% Print summary data
    disp(['id=' run_id]);
    disp(['n = ' num2str(mean(n)) ' +/- ' num2str(1*std(n))]);
    disp(['Q_dis = ' num2str(mean(Q_dis)/1e3) ' +/- ' num2str(1*std(Q_dis)/1e3) ' kJ/mol']);
    disp(['V_dis = ' num2str(mean(V_dis)*1e6) ' +/- ' num2str(1*std(V_dis)*1e6) ' m^3/mol']);
    disp(['log10(A_dis) = ' num2str(mean(log10(A_dis))) ' +/- ' num2str(1*std(log10(A_dis)))]);
    disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(1*std(chi2))]);

    disp(['m = ' num2str(mean(m)) ' +/- ' num2str(1*std(m))]);
    disp(['Q_dif = ' num2str(mean(Q_dif)/1e3) ' +/- ' num2str(1*std(Q_dif)/1e3) ' kJ/mol']);
    disp(['log10(A_dif) = ' num2str(mean(log10(A_dif))) ' +/- ' num2str(1*std(log10(A_dif)))]);
    disp(['V_dif = ' num2str(mean(V_dif)*1e6) ' +/- ' num2str(1*std(V_dif)*1e6) ' m^3/mol']);
    for i=2:5
        disp(['X_' num2str(i) ' = ' num2str(mean(X{i})) ' +/- ' num2str(1*std(X{i}))]);
    end

% set constants
    R = 8.3145;
    nx = 20;
    %stress data
    xs = linspace(min(s_all)-20,max(s_all)+150,nx);
    %temperature data
    xt = linspace(min(T_all)-100,max(T_all)+100,nx);
    %grain size data
    xd = linspace(min(d_all)-50,max(d_all)+100,nx);
    %pressure data
    xP = linspace(min(P_all)-0.25*1e9,max(P_all)+0.5*1e9,nx);
    %normalize data
    normT = mean(T_all); % Temperature [K]
    norms = mean(s_all); % Stress [MPa]
    normP = mean(P_all); % Pressure [Pa]
    normd = mean(d_all); % Grain size [um]
    norm_params = [normT, normP, norms, normd];

% Strain vs. stress
figure(1); hold off;
for i=1:nout
    %model strain rate
    e_dis = A_dis(i)*xs.^n(i).*exp(-(Q_dis(i)+normP*V_dis(i))/(normT*R));
    e_dif = A_dif(i)*xs.*normd^-m(i).*exp(-(Q_dif(i)+normP*V_dif(i))/(normT*R));
    loglog(xs,e_dis+e_dif,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('Strain rate [s^{-1}]')
        xlabel('Stress [MPa]')
        title('Strain rate vs. Stress')
        set(gca, 'XTick', [50:50:350])
    end
end

for i=1:5
    norm_factor = calc_norm("s",out,file,norm_params);
    loglog(s{i},e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
    %loglog(s{i},e{i},'g*');
    % error bars
    [ex{i},ey{i}] = calc_err(s{i},ds{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    loglog(ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end
% GT95's flow law -- difference due to activation volume?
loglog(xs, exp(-7.58).*xs.^4.*exp(-223000/(R*normT)),'c', 'Linewidth', 1)

figure(2); hold off;
% Strain vs. temp
for i=1:nout
    %model strain rate
    e_dis = A_dis(i)*norms.^n(i).*exp(-(Q_dis(i)+normP*V_dis(i))./(xt.*R));
    e_dif = A_dif(i)*norms.*normd^-m(i).*exp(-(Q_dif(i)+normP*V_dif(i))./(xt.*R));
    semilogy(1e4./xt,e_dis + e_dif,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('Strain rate [s^{-1}]')
        xlabel('10^4/T [K^{-1}]')
        title('Strain rate. vs. Temperature')
    end
end

for i=1:5
    norm_factor = calc_norm("T",out,file,norm_params);
    semilogy(1e4./T{i},e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
    %semilogy(1e4./T{i},e{i},'g*');
    % error bars
    [ex{i},ey{i}] = calc_err(T{i},dT{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    semilogy(1e4./ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end
loglog(1e4./xt, exp(-7.58)*norms^4.*exp(-223000./(R.*xt)),'c', 'Linewidth', 1)

figure(3); hold off;
% Strain vs. grain size
for i=1:nout
    %model strain rate
    e_dis = A_dis(i)*norms.^n(i).*exp(-(Q_dis(i)+normP*V_dis(i))./(normT.*R));
    e_dif = A_dif(i)*norms*xd.^-m(i).*exp(-(Q_dif(i)+normP*V_dif(i))./(normT.*R));
    loglog(xd,e_dis + e_dif,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('Strain rate [s^{-1}]')
        xlabel('Grain size [um]')
    end
end

for i=1:5
    norm_factor = calc_norm("d",out,file,norm_params);
    loglog(d{i},e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
    %loglog(d{i},e{i},'g*');
    % error bars
    [ex{i},ey{i}] = calc_err(d{i},dd{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    loglog(ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end
loglog(xd, ones(size(xd))*exp(-7.58)*norms^4*exp(-223000/(R*normT)),'c', 'Linewidth', 1)
%{
figure(4); hold off;
% Strain vs. pressure
for i=1:nout
    %model strain rate
    e_dis = A_dis(i)*norms.^n(i).*exp(-(Q_dis(i)+xP.*V_dis(i))./(normT.*R));
    e_dif = A_dif(i)*norms*normd.^-m(i).*exp(-(Q_dif(i)+xP.*V_dif(i))./(normT.*R));
    semilogy(xP./1e9,e_dis + e_dif,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('Strain rate [s^{-1}]')
        xlabel('Pressure [GPa]')
        title('Strain rate. vs. Pressure')
    end
end
for i=1:5
    norm_factor = calc_norm("P",out,file,norm_params);
    semilogy(P{i}./1e9,e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
    % error bars
    [ex{i},ey{i}] = calc_err(P{i},dP{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    semilogy(ex{i}./1e9,ey{i}.*exp(-mean(X{i})),'b-');
end

% Inter-run bias
%{
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
%}
%}
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
 subplot(3,3,7);
 hist(V_dis*1e6,20); xlabel('V_{dis} [m^3/mol]');
 subplot(3,3,8);
 hist(V_dif*1e6,20); xlabel('V_{dif} [m^3/mol]');
 subplot(3,3,9);
 hist(chi2,20); xlabel('\chi^2');

 %{
figure(5); hold off;
subplot(1,2,1);
for i=1:5
    norm_factor = calc_norm("s",out,file,norm_params);
    loglog(s{i},e{i}.*norm_factor{i},'-o');
    hold on;
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([0,250])
title("Without inter-run bias")
subplot(1,2,2);
for i=1:5
    norm_factor = calc_norm("s",out,file,norm_params);
    loglog(s{i},e{i}.*norm_factor{i}*exp(-mean(X{i})),'-o');
    hold on;
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([0,250])
title("With inter-run bias")
%}