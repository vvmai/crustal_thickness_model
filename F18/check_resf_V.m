% check_resp.m
clear
run_id = 'fukuda18';
% load data file
    file = load(['./' run_id '_V.dat']);
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
    output = load(['./' run_id '_f_V.out']);
    % skip first 1000 runs, then select results from every 10 runs
    out = output(100:10:end,1:end);
    %out = output;
    nout = length(out);
    chi2 = out(:,3);
    m = out(:,4); n = out(:,5); r = out(:,6); Q = out(:,7); V = out(:,8);
    A = out(:,end);

% Print summary data
    disp(['id=' run_id]);
    disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(1*std(chi2))]);
    disp(['m = ' num2str(mean(m)) ' +/- ' num2str(1*std(m))]);
    disp(['n = ' num2str(mean(n)) ' +/- ' num2str(1*std(n))]);
    disp(['r = ' num2str(mean(r)) ' +/- ' num2str(1*std(r))]);
    disp(['Q = ' num2str(mean(Q)/1e3) ' +/- ' num2str(1*std(Q)/1e3) ' kJ/mol']);
    disp(['V = ' num2str(mean(V)*1e6) ' +/- ' num2str(1*std(V)*1e6) ' m^3/mol']);
    disp(['log10(A) = ' num2str(mean(log10(A))) ' +/- ' num2str(1*std(log10(A)))]);

% set constants
    R = 8.3145;
    nx = 20;
    %stress data
    xs = logspace(min(s_all)-30,max(s_all)+600,nx);
    %temperature data
    xt = linspace(min(T_all)-100,max(T_all)+200,nx);
    %grain size data
    xd = linspace(min(d_all)-2,max(d_all)+100,nx);
    %pressure data
    xP = linspace(min(P_all)-0.25*1e9,max(P_all)+0.25*1e9,nx);
    %water fugacity data
    xf = linspace(min(f_H2O_all)-300,max(f_H2O_all)+1000,nx);
    %normalize data
    normT = 1073; % Temperature [K]
    norms = 380; % Stress [MPa]
    normP = 1.5*1e9; % Pressure [Pa]
    normd = 18; % Grain size [um]
    normf_H2O = 4.5*1e3; % Water fugacity [MPa]
    norm_params = [normT, normP, norms, normd, normf_H2O];

% Strain vs. stress
figure(1); hold off;
    for i=1:nout
        %model strain rate
        e_pred = A(i)*xs.^n(i).*normd^-m(i)*normf_H2O^r(i)*exp(-(Q(i)+normP*V(i))/(normT*R));
        loglog(xs, e_pred,'r:');
        if i==1
            hold on; box on;
            ylabel('Strain rate [s^{-1}]')
            xlabel('Stress [MPa]')
            title('Strain rate vs. Stress')
            xlim([10,1500])
            ylim([1e-9,1e-1])
        end
    end
    for i=1:file(end,end)
        % original data
        loglog(s{i},e{i},'go');
        [ex_0{i},ey_0{i}] = calc_err(s{i},ds{i},e{i},(de{i}./e{i}).*e{i});
        loglog(ex_0{i},ey_0{i},'g-');

        % normalized data
        norm_factor{i} = (s{i}.^mean(n).*normd^-mean(m)*normf_H2O^mean(r)*exp(-(mean(Q)+normP*mean(V))/(normT*R)))./...
        (s{i}.^mean(n).*d{i}.^-mean(m).*f_H2O{i}.^mean(r).*exp(-(mean(Q)+P{i}.*mean(V))./(T{i}.*R)));
        loglog(s{i},e{i}.*norm_factor{i},'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(s{i},ds{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        loglog(ex{i},ey{i},'b-');
    end
%}
%{
% Strain vs. temp
figure(2); hold off;
    for i=1:nout
        %model strain rate
        e_pred = A(i)*norms.^n(i).*normd^-m(i)*normf_H2O^r(i).*exp(-(Q(i)+normP*V(i))./(xt.*R));
        semilogy(1e4./xt, e_pred,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('Strain rate [s^{-1}]')
            xlabel('10^4/T [K^{-1}]')
            title('Strain rate vs. Temperature')
            ylim([1e-7,1e-3])
        end
    end

    for i=1:13
        % original data
        semilogy(1e4./T{i}+(0.01*i),e{i},'go');
        [ex_0{i},ey_0{i}] = calc_err(T{i},dT{i},e{i},(de{i}./e{i}).*e{i});
        semilogy(1e4./ex_0{i}+(0.01*i),ey_0{i},'g-');

        % normalized data
        norm_factor{i} = (norms^mean(n).*normd^-mean(m).*normf_H2O^mean(r).*exp(-(mean(Q)+normP*mean(V))./(T{i}.*R)))./(s{i}.^mean(n).*d{i}.^-mean(m).*f_H2O{i}.^mean(r).*exp(-(mean(Q)+P{i}.*mean(V))./(T{i}.*R)));
        semilogy(1e4./T{i}+(0.01*i),e{i}.*norm_factor{i},'bo');
        % error bars
        [ex{i}, ey{i}] = calc_err(T{i}, dT{i}, e{i}.*norm_factor{i}, de{i}.*norm_factor{i});
        semilogy(1e4./ex{i}+(0.01*i),ey{i},'b-');
    end

% Strain vs. grain size
figure(3); hold off;
    for i=1:nout
        %model strain rate
        e_pred = A(i)*norms.^n(i).*xd.^-m(i)*normf_H2O^r(i).*exp(-(Q(i)+normP*V(i))./(normT.*R));
        loglog(xd,e_pred,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('Strain rate [s^{-1}]')
            xlabel('Grain size [um]')
            title('Strain rate vs. Grain size')
            ylim([1e-7,1e-2])
        end
    end
    
    for i=1:file(end,end)
        % original data
        loglog(d{i},e{i},'go');
        % error bars
        [ex_0{i},ey_0{i}] = calc_err(d{i},dd{i},e{i},(de{i}./e{i}).*e{i});
        loglog(ex_0{i},ey_0{i},'g-');

        % normalized data
        norm_factor{i} = (norms^mean(n)*normf_H2O^mean(r).*exp(-(mean(Q)+normP*mean(V))./(normT*R)))./...
        (s{i}.^mean(n).*f_H2O{i}.^mean(r).*exp(-(mean(Q)+P{i}.*mean(V))./(T{i}.*R)));
        loglog(d{i},e{i}.*norm_factor{i},'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(d{i},dd{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        loglog(ex{i},ey{i},'b-');
    end

% Strain vs. pressure
figure(4); hold off;
    for i=1:nout
        %model strain rate
        e_pred = A(i)*norms.^n(i).*normd.^-m(i)*normf_H2O^r(i).*exp(-(Q(i)+xP.*V(i))./(normT.*R));
        semilogy(xP./1e9,e_pred,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('Strain rate [s^{-1}]')
            xlabel('Pressure [GPa]')
            title('Strain rate. vs. Pressure')
            ylim([1e-7,1e-2])
        end
    end
    
    for i=1:file(end,end)
        % original data
        semilogy(P{i}./1e9+(0.005*i),e{i},'go');
        % error bars
        [ex_0{i},ey_0{i}] = calc_err(P{i},dP{i},e{i},(de{i}./e{i}).*e{i});
        semilogy(ex_0{i}./1e9+(0.005*i),ey_0{i},'g-');

        % normalized data
        norm_factor{i} = (norms^mean(n)*normd^-mean(m)*normf_H2O^mean(r).*exp(-(mean(Q)+P{i}.*mean(V))./(normT*R)))./...
        (s{i}.^mean(n).*d{i}.^-mean(m).*f_H2O{i}.^mean(r).*exp(-(mean(Q)+P{i}.*mean(V))./(T{i}.*R)));
        semilogy(P{i}./1e9+(0.005*i),e{i}.*norm_factor{i},'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(P{i},dP{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        semilogy(ex{i}./1e9+(0.005*i),ey{i},'b-');
    end

% Strain vs. water fugacity
figure(5); hold off;
    for i=1:nout
        %model strain rate
        e_pred = A(i)*norms.^n(i).*normd.^-m(i).*xf.^r(i).*exp(-(Q(i)+normP.*V(i))./(normT.*R));
        loglog(xf,e_pred,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('Strain rate [s^{-1}]')
            xlabel('Water fugacity [MPa]')
            title('Strain rate. vs. Water fugacity')
            ylim([1e-7,1e-2])
        end
    end
    
    for i=1:file(end,end)
        % original data
        loglog(f_H2O{i}+(10*i),e{i},'go');
        % error bars
        [ex_0{i},ey_0{i}] = calc_err(f_H2O{i},df_H2O{i},e{i},(de{i}./e{i}).*e{i});
        loglog(ex_0{i}+(10*i),ey_0{i},'g-');

        % normalized data
        norm_factor{i} = (norms^mean(n)*normd^-mean(m).*exp(-(mean(Q)+normP.*mean(V))./(normT*R)))./...
        (s{i}.^mean(n).*d{i}.^-mean(m).*exp(-(mean(Q)+P{i}.*mean(V))./(T{i}.*R)));
        loglog(f_H2O{i}+(10*i),e{i}.*norm_factor{i},'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(f_H2O{i},df_H2O{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        loglog(ex{i}+(10*i),ey{i},'b-');
    end

% Histograms
figure(6); hold off;
 subplot(4,2,1);
 hist(n,20); xlabel('n');
 subplot(4,2,2);
 hist(r,20); xlabel('r');
 subplot(4,2,3);
 hist(m,20); xlabel('m');
 subplot(4,2,4);
 hist(Q/1e3,20); xlabel('Q [kJ/mol]');
 subplot(4,2,5);
 hist(log10(A),20); xlabel('log_{10}(A)');
 subplot(4,2,6);
 hist(V*1e6,20); xlabel('V [m^3/mol]');
 subplot(4,2,7);
 hist(chi2,20); xlabel('\chi^2');

 %1:1 plot
 figure(7); hold off;
 xe = logspace(log10(min(e_all))-1,(log10(max(e_all)))+1,nx);
 for (i = 1:file(end,end))
    e_predict{i} = 10^mean(log10(A)).*f_H2O{i}.^mean(r).*d{i}.^-mean(m).*s{i}.^mean(n).*exp(-(mean(Q)+P{i}.*mean(V))./(T{i}.*R));
    loglog(e_predict{i},e{i},'*');
    if i ==1
        hold on; axis tight;
        ylabel('Observed strain rate [s^{-1}]')
        xlabel('Predicted strain rate [s^{-1}]')
    end
    errorbar(e_predict{i},e{i},de{i},'k*');
end
loglog(xe,xe,'k-');

%}