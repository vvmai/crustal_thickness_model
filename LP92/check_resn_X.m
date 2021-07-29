% check_resp.m
clear
run_id = 'LP92';

file = load(['./' run_id '_COH.dat']);
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
    COH{i} = data(:,11);
    dCOH{i} = data(:,12);
end

output = load(['./' 'LP92' '_COH.out']);
% skip first 1000 runs, then select results from every 100 runs
out = output(100:10:end,1:end);
nout = length(out);
% 4 = Water exponent, 5 = Stress exponent, 6 = Activation energy (J/mol), 21 = A
chi2 = out(:,3); r = out(:,4); n = out(:,5); Q = out(:,6); A = out(:,21);
% 7-20 = Inter-run bias relative to sample 1
for i=7:20
    X{i-6} = out(:,i);
end

% Print summary data
disp(['id=' run_id]);
disp(['n = ' num2str(mean(n)) ' +/- ' num2str(2*std(n))]);
disp(['r = ' num2str(mean(r)) ' +/- ' num2str(2*std(r))]);
disp(['Q = ' num2str(mean(Q)/1e3) ' +/- ' num2str(2*std(Q)/1e3) ' kJ/mol']);
disp(['log10(A) = ' num2str(mean(log10(A))) ' +/- ' num2str(2*std(log10(A)))]);
disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(2*std(chi2))]);
for i=2:14
    disp(['X_' num2str(i) ' = ' num2str((mean(X{i})))]);
end

R = 8.3145;
nx = 100;
%stress data
xs = linspace(min(s_all)-100,max(s_all)+1000,nx);
%temperature data
xt = linspace(min(T_all)-100,max(T_all)+300,nx);
%COH data
xCOH = logspace(log10(min(COH_all))-2,log10(max(COH_all))+3,nx);
%normalize data
normT = 1373;
norms = 250;
normCOH = 4500;

mat = [log10(A) n r Q];
[cc] = corrcoef(mat)

figure(1);
% Strain vs. stress
for i=1:nout
    %model strain rate
    e_pred = A(i)*xs.^n(i).*normCOH^r(i).*exp(-Q(i)/(normT*R));
    loglog(xs,e_pred,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('strain rate [s^{-1}]')
        xlabel('stress [MPa]')
        title('LP92')
        xlim([10^(1.8),10^(3.2)])
        ylim([10^(-6.5),1e-1])
        set(gca, 'LineWidth', 1)
        set(gca,'FontSize', 20)
        text(10^(1.85), 10^(-1.9), ['T = ' num2str(normT) ' K' newline 'C_{OH} = ' num2str(normCOH) ' H/10^6Si'], 'FontSize', 20)
    end
end

for i=1:14
    norm_factor{i} = (normCOH^mean(r)*exp(-mean(Q)/(normT*R)))./(COH{i}.^mean(r).*exp(-mean(Q)./(T{i}.*R)));
    loglog(s{i},e{i}.*norm_factor{i}*exp(-mean(X{i})),'bo');
    % error bars
    [ex{i},ey{i}] = calc_err(s{i},ds{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    loglog(ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end

loglog(xs, 1.2e-8.*xs.^4.*exp(-152000/(normT*R)),'c', 'Linewidth', 1.5);
loglog(xs, 10^mean(log10(A)).*xs.^mean(n).*normCOH^mean(r).*exp(-mean(Q)./(normT.*R)), 'k--','Linewidth', 1.5);

%xlim([100,700])

figure(2); hold off;
% Strain vs. temp
for i=1:nout
    %model strain rate
    e_pred = A(i)*norms^n(i)*normCOH^r(i).*exp(-Q(i)./(xt.*R));
    loglog(10^4./xt,e_pred,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('strain rate [s^{-1}]')
        xlabel('10^4/T [K^{-1}]')
        ylim([1e-7 1e-2])
        xlim([7.5,9.5])
        set(gca, 'LineWidth', 1)
        set(gca,'FontSize',20)
        text(7.55, 10^(-2.7), ['\sigma = ' num2str(norms) ' MPa' newline 'C_{OH} = ' num2str(normCOH) ' H/10^6Si'], 'FontSize', 20)
    end
end

for i=1:14
    norm_factor{i} = (normCOH^mean(r)*norms^mean(n))./(COH{i}.^mean(r).*s{i}.^mean(n));
    semilogy((0.999+i*0.0005)*10^4./T{i},e{i}.*norm_factor{i}*exp(-mean(X{i})),'bo');
    % error bars
    [ex{i},ey{i}] = calc_err(T{i},dT{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    semilogy((0.999+i*0.0005)*10^4./ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end

semilogy(1e4./xt, 1.2e-8*norms^4.*exp(-152000./(xt.*R)),'c', 'Linewidth', 1.5);
semilogy(1e4./xt, 10^mean(log10(A))*norms^mean(n).*normCOH^mean(r).*exp(-mean(Q)./(xt.*R)), 'k--','Linewidth', 1.5);


figure(3); hold off;
% Strain vs. COH
for i=1:nout
    %model strain rate
    e_pred = A(i)*norms^n(i).*xCOH.^r(i).*exp(-Q(i)/(normT*R));
    loglog(xCOH,e_pred,'r:');
    if i==1
        hold on; box on; axis tight;
        set(gca, 'LineWidth', 1)
        set(gca,'FontSize',20)
        ylim([1e-6 1e-2])
        xlim([1e2 1e5])
        text(10^(2.1), 10^(-2.5), ['T = ' num2str(normT) ' K' newline '\sigma = ' num2str(norms) ' MPa'], 'FontSize', 20)
    end
end

for i=1:14
    norm_factor{i} = (norms^mean(n)*exp(-mean(Q)/(normT*R)))./(s{i}.^mean(n).*exp(-mean(Q)./(T{i}.*R)));
    loglog(COH{i},e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
    % error bars
    [ex{i},ey{i}] = calc_err(COH{i},dCOH{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    loglog(ex{i},ey{i}.*exp(-mean(X{i})),'b-');
end

loglog(xCOH, ones(size(xCOH)).*1.2e-8*norms^4.*exp(-152000./(normT.*R)),'c', 'Linewidth', 1.5);
loglog(xCOH, 10^mean(log10(A))*norms^mean(n).*xCOH.^mean(r).*exp(-mean(Q)/(normT*R)), 'k--','Linewidth', 1.5);

ylabel('Strain rate [s^{-1}]')
xlabel('C_{OH} [ppm H/Si]')
%}
%{
figure(4); hold off;
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
 
figure(4);
subplot(3,2,1);
hist(n,20); xlabel('n');
subplot(3,2,2);
hist(r,20); xlabel('r');
subplot(3,2,3);
hist(Q/1e3,20); xlabel('Q [kJ/mol]');
subplot(3,2,4);
hist(log10(A),20); xlabel('log_{10}(A)');
subplot(3,2,5);
hist(chi2,20); xlabel('\chi^2');
%}