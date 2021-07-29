% check_resp.m
clear
run_id = 'fukuda18';
% load data file
    file = load(['./' run_id '_d_fixed.dat']);
    T_all = file(:,1); % Temperature (K)
    P_all = file(:,3)*1e9; % Pressure (Pa)
    e_all = file(:,5); % Strain rate (s^-1)
    s_all = file(:,7); % Stress (MPa)
    d_all = file(:,9); % grain size (um)
    f_H2O_all = file(:,11); % water fugacity [MPa]

% Separate samples into arrays nested in a cell array
for i=1:file(end,end)
    data = file(file(:,end) == i,:);
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
    f_H2O{i} = data(:,11); % water fugacity [MPa]
    df_H2O{i} = data(:,12);
end
% get output data
    output = load(['./' run_id '_f_d_fixed.out']);
    % skip first 1000 runs, then select results from every 10 runs
    out = output(100:10:end,1:end);
    %out = output;
    nout = length(out);
    chi2 = out(:,3);
    m = out(:,4); n = out(:,5); r = out(:,6); Q = out(:,7);
    A = out(:,end);
    % 8-21 = Inter-run bias relative to sample 1

    mat = [log10(A) n m r Q];
    [cc] = corrcoef(mat)

% Print summary data
    disp(['id=' run_id]);
    disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(1*std(chi2))]);
    disp(['m = ' num2str(mean(m)) ' +/- ' num2str(1*std(m))]);
    disp(['n = ' num2str(mean(n)) ' +/- ' num2str(1*std(n))]);
    disp(['r = ' num2str(mean(r)) ' +/- ' num2str(1*std(r))]);
    disp(['Q = ' num2str(mean(Q)/1e3) ' +/- ' num2str(1*std(Q)/1e3) ' kJ/mol']);
    disp(['log10(A) = ' num2str(mean(log10(A))) ' +/- ' num2str(1*std(log10(A)))]);

% set constants
    R = 8.3145;
    nx = 20;
    %stress data
    xs = logspace(log10(min(s_all))-0.5,log10(max(s_all))+1,nx);
    %temperature data
    xt = linspace(min(T_all)-100,max(T_all)+600,nx);
    %grain size data
    xd = logspace(log10(min(d_all))-2,log10(max(d_all))+2,nx);
    %pressure data
    xP = linspace(min(P_all)-0.25*1e9,max(P_all)+0.5*1e9,nx);
    %water fugacity data
    xf = logspace(log10(min(f_H2O_all))-1,log10(max(f_H2O_all))+1,nx);
    %normalize data
    normT = 1073; % Temperature [K]
    norms = 250; % Stress [MPa]
    normP = 1.5*1e9; % Pressure [Pa]
    normd = 10; % Grain size [um]
    normf_H2O = 4.5*1e3; % Water fugacity [MPa]
    norm_params = [normT, normP, norms, normd, normf_H2O];

% Strain vs. stress
figure(1); hold off;
    for i=1:nout
        %model strain rate
        e_pred = A(i)*xs.^n(i).*normd^-m(i)*normf_H2O^r(i)*exp(-(Q(i))/(normT*R));
        loglog(xs, e_pred,'r:');
        if i==1
            hold on; box on;
            ylabel('strain rate [s^{-1}]')
            xlabel('stress [MPa]')
            title('F18 (GBS)')
            set(gca, 'LineWidth', 1)
            set(gca,'FontSize',20)
            ylim([10^(-7.5) 1e0])
            xlim([1e1,0.5*1e4])
            text(12, 10^(-1.8),['T = ' num2str(normT) ' K' newline 'd = ' num2str(normd) ' \mu m' newline 'f_{H_2O} = ' num2str(normf_H2O) ' MPa'], 'FontSize', 20)
        end
    end
    for i=1:file(end,end)
        % original data
        %loglog(s{i},e{i},'go');
        %[ex_0{i},ey_0{i}] = calc_err(s{i},ds{i},e{i},(de{i}./e{i}).*e{i});
        %loglog(ex_0{i},ey_0{i},'g-');

        % normalized data
        norm_factor{i} = (s{i}.^mean(n).*normd^-mean(m)*normf_H2O^mean(r)*exp(-(mean(Q))/(normT*R)))./...
        (s{i}.^mean(n).*d{i}.^-mean(m).*f_H2O{i}.^mean(r).*exp(-(mean(Q))./(T{i}.*R)));
        loglog(s{i},e{i}.*norm_factor{i},'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(s{i},ds{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        loglog(ex{i},ey{i},'b-');
    end
loglog(xs, 10^-2.97.*xs.^(1.7).*normd^-0.51.*normf_H2O.*exp(-183000/(R*normT)),'c', 'Linewidth', 1.5);
loglog(xs, 10^mean(log10(A)).*xs.^mean(n).*normd^-mean(m).*normf_H2O^mean(r).*exp(-mean(Q)./(normT.*R)), 'k','Linewidth',1.5);    
%}

% Strain vs. temp
figure(2); hold off;
    for i=1:nout
        %model strain rate
        e_pred = A(i)*norms.^n(i).*normd^-m(i)*normf_H2O^r(i).*exp(-(Q(i))./(xt.*R));
        semilogy(1e4./xt, e_pred,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('strain rate [s^{-1}]')
            xlabel('10^4/T [K^{-1}]')
            set(gca, 'LineWidth', 1)
            set(gca,'FontSize', 20)
            xlim([7.5 11.5])
            ylim([1e-7,1e0])
            text(9.4, 10^(-1.65),['\sigma = ' num2str(norms) ' MPa' newline 'd = ' num2str(normd) ' \mu m' newline 'f_{H_2O} = ' num2str(normf_H2O) ' MPa'], 'FontSize', 20)
        end
    end

    for i=1:13
        % original data
        %semilogy(1e4./T{i}+(0.01*i),e{i},'go');
        %[ex_0{i},ey_0{i}] = calc_err(T{i},dT{i},e{i},(de{i}./e{i}).*e{i});
        %semilogy(1e4./ex_0{i}+(0.01*i),ey_0{i},'g-');

        % normalized data
        norm_factor{i} = (norms^mean(n).*normd^-mean(m).*normf_H2O^mean(r).*exp(-(mean(Q))./(T{i}.*R)))./(s{i}.^mean(n).*d{i}.^-mean(m).*f_H2O{i}.^mean(r).*exp(-(mean(Q))./(T{i}.*R)));
        semilogy(1e4./T{i}+(0.01*i),e{i}.*norm_factor{i},'bo');
        % error bars
        [ex{i}, ey{i}] = calc_err(T{i}, dT{i}, e{i}.*norm_factor{i}, de{i}.*norm_factor{i});
        semilogy(1e4./ex{i}+(0.01*i),ey{i},'b-');
    end
semilogy(1e4./xt, 10^-2.97.*norms^(1.7).*normd^-0.51.*normf_H2O.*exp(-183000./(R.*xt)),'c', 'Linewidth', 1.5);
semilogy(1e4./xt, 10^mean(log10(A)).*norms.^mean(n).*normd^-mean(m).*normf_H2O^mean(r).*exp(-mean(Q)./(xt.*R)), 'k','Linewidth',1.5);     

% Strain vs. grain size
figure(3); hold off;
    for i=1:nout
        %model strain rate
        e_pred = A(i)*norms.^n(i).*xd.^-m(i)*normf_H2O^r(i).*exp(-(Q(i))./(normT.*R));
        loglog(xd,e_pred,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('strain rate [s^{-1}]')
            xlabel('grain size [\mum]')
            set(gca, 'LineWidth', 1)
            set(gca,'FontSize',20)
            ylim([1e-6,10^(-1.5)])
            xlim([1e0 1e2])
            text(1.1, 10^(-2.5), ['\sigma = ' num2str(norms) ' MPa' newline 'T = ' num2str(normT) ' K' newline 'f_{H_2O} = ' num2str(normf_H2O) ' MPa'], 'FontSize', 20)
        end
    end
    
    for i=1:file(end,end)
        % original data
        %loglog(d{i},e{i},'go');
        % error bars
        %[ex_0{i},ey_0{i}] = calc_err(d{i},dd{i},e{i},(de{i}./e{i}).*e{i});
        %loglog(ex_0{i},ey_0{i},'g-');

        % normalized data
        norm_factor{i} = (norms^mean(n)*normf_H2O^mean(r).*exp(-(mean(Q))./(normT*R)))./...
        (s{i}.^mean(n).*f_H2O{i}.^mean(r).*exp(-(mean(Q))./(T{i}.*R)));
        loglog(d{i},e{i}.*norm_factor{i},'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(d{i},dd{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        loglog(ex{i},ey{i},'b-');
    end
loglog(xd, 10^-2.97.*norms.^(1.7).*xd.^-0.51.*normf_H2O.*exp(-183000/(R*normT)),'c', 'Linewidth', 1.5);
loglog(xd, 10^mean(log10(A)).*norms.^mean(n).*xd.^-mean(m).*normf_H2O^mean(r).*exp(-mean(Q)./(normT.*R)), 'k','Linewidth',1.5);       

% Strain vs. water fugacity
figure(4); hold off;
for i=1:nout
    %model strain rate
    e_pred = A(i)*norms.^n(i).*normd.^-m(i).*xf.^r(i).*exp(-(Q(i))./(normT.*R));
    loglog(xf,e_pred,'r:');
    if i==1
        hold on; box on; axis tight;
        ylabel('strain rate [s^{-1}]')
        xlabel('water fugacity [MPa]')
        ylim([1e-7,1e-2])
        xlim([1e3 1e4])
        set(gca, 'LineWidth', 1)
        set(gca,'FontSize',20)
        text(10^(3.05), 10^(-3), ['\sigma = ' num2str(norms) ' MPa' newline 'T = ' num2str(normT) ' K' newline 'd = ' num2str(normd) ' \mu m'], 'FontSize', 20)
    end
end

for i=1:file(end,end)
    % original data
    %loglog(f_H2O{i}+(10*i),e{i},'go');
    % error bars
    %[ex_0{i},ey_0{i}] = calc_err(f_H2O{i},df_H2O{i},e{i},(de{i}./e{i}).*e{i});
    %loglog(ex_0{i}+(10*i),ey_0{i},'g-');

    % normalized data
    norm_factor{i} = (norms^mean(n)*normd^-mean(m).*exp(-(mean(Q))./(normT*R)))./...
    (s{i}.^mean(n).*d{i}.^-mean(m).*exp(-(mean(Q))./(T{i}.*R)));
    loglog(f_H2O{i}+(10*i),e{i}.*norm_factor{i},'bo');
    % error bars
    [ex{i},ey{i}] = calc_err(f_H2O{i},df_H2O{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
    loglog(ex{i}+(10*i),ey{i},'b-');
end
loglog(xf, 10^-2.97.*norms.^(1.7).*normd.^-0.51.*xf.*exp(-183000/(R*normT)),'c', 'Linewidth', 1);
loglog(xf, 10^mean(log10(A)).*norms.^mean(n).*normd.^-mean(m).*xf.^mean(r).*exp(-mean(Q)./(normT.*R)), 'k');       
%}
%{
% Histograms
figure(7);
 subplot(3,3,1);
 hist(n,20); xlabel('n');
 subplot(3,3,2);
 hist(r,20); xlabel('r');
 subplot(3,3,3);
 hist(m,20); xlabel('m');
 subplot(3,3,4);
 hist(Q/1e3,20); xlabel('Q [kJ/mol]');
 subplot(3,3,5);
 hist(log10(A),20); xlabel('log_{10}(A)');
 subplot(3,3,7);
 hist(chi2,20); xlabel('\chi^2');
 %}
 %}