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
    output = load(['./' run_id '_V.out']);
    % skip first 1000 runs, then select results from every 10 runs
    out = output(100:5:end,1:end);
    %out = output;
    nout = length(out);
    chi2 = out(:,3);
    m = out(:,4); r_dif = out(:,5); Q_dif = out(:,6); V_dif = out(:,7);
    r_dis = out(:,8); n = out(:,9);Q_dis = out(:,10); V_dis = out(:,11);
    A_dif = out(:,end-1); A_dis = out(:,end);
    % 12-24 = Inter-run bias relative to sample 1
    for i=12:(size(out,2)-2)
        X{i-11} = out(:,i);
    end

% Print summary data
    disp(['id=' run_id]);
    disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(1*std(chi2))]);

    disp(['n = ' num2str(mean(n)) ' +/- ' num2str(1*std(n))]);
    disp(['r_dis = ' num2str(mean(r_dis)) ' +/- ' num2str(1*std(r_dis))]);
    disp(['Q_dis = ' num2str(mean(Q_dis)/1e3) ' +/- ' num2str(1*std(Q_dis)/1e3) ' kJ/mol']);
    disp(['V_dis = ' num2str(mean(V_dis)*1e6) ' +/- ' num2str(1*std(V_dis)*1e6) ' m^3/mol']);
    disp(['log10(A_dis) = ' num2str(mean(log10(A_dis))) ' +/- ' num2str(1*std(log10(A_dis)))]);

    disp(['m = ' num2str(mean(m)) ' +/- ' num2str(1*std(m))]);
    disp(['r_dif = ' num2str(mean(r_dif)) ' +/- ' num2str(1*std(r_dif))]);
    disp(['Q_dif = ' num2str(mean(Q_dif)/1e3) ' +/- ' num2str(1*std(Q_dif)/1e3) ' kJ/mol']);
    disp(['log10(A_dif) = ' num2str(mean(log10(A_dif))) ' +/- ' num2str(1*std(log10(A_dif)))]);
    disp(['V_dif = ' num2str(mean(V_dif)*1e6) ' +/- ' num2str(1*std(V_dif)*1e6) ' m^3/mol']);
    for i=2:length(X)
        disp(['X_' num2str(i) ' = ' num2str(mean(X{i})) ' +/- ' num2str(1*std(X{i}))]);
    end

% set constants
    R = 8.3145;
    nx = 20;
    %stress data
    xs = logspace(log10(min(s_all))-1,(log10(max(s_all))+2),nx);
    %temperature data
    xt = linspace(min(T_all)-100,max(T_all)+200,nx);
    %grain size data
    xd = logspace(log10(min(d_all))-1,log10(max(d_all))+2,nx);
    %pressure data
    xP = linspace(min(P_all)-0.25*1e9,max(P_all)+0.5*1e9,nx);
    %water fugacity data
    xf = logspace(log10(min(f_H2O_all))-2,log10(max(f_H2O_all))+3,nx);
    %normalize data
    normT = 1073; % Temperature [K]
    norms = 200; % Stress [MPa]
    normP = 1.5*1e9; % Pressure [Pa]
    normd = 7; % Grain size [um]
    normf_H2O = 4.5*1e3; % Water fugacity [MPa]
    norm_params = [normT, normP, norms, normd, normf_H2O];

% Strain vs. stress
figure(1); hold off;
    for i=1:nout
        %model strain rate
        e_dis = A_dis(i)*xs.^n(i).*normf_H2O^r_dis(i)*exp(-(Q_dis(i)+normP*V_dis(i))/(normT*R));
        e_dif = A_dif(i)*xs.*normd^-m(i).*normf_H2O^r_dif(i)*exp(-(Q_dif(i)+normP*V_dif(i))/(normT*R));
        loglog(xs, e_dis + e_dif,'r:');
        if i==1
            hold on; box on;
            ylabel('Strain rate [s^{-1}]')
            xlabel('Stress [MPa]')
            title('Strain rate vs. Stress')
            xlim([20,2000])
            ylim([1e-7,1e-1])
        end
    end

    for i=1:size(X,2)
        % normalize
        norm_factor = calc_norm("wet","s",out,file,norm_params, "V");
        loglog(s{i},e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(s{i},ds{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        loglog(ex{i},ey{i}.*exp(-mean(X{i})),'b-');
        % original data
        loglog(s{i},e{i}.*exp(-mean(X{i})),'go');
        [ex_0{i},ey_0{i}] = calc_err(s{i},ds{i},e{i},(de{i}./e{i}).*e{i});
        loglog(ex_0{i},ey_0{i}.*exp(-mean(X{i})),'g-');
    end
    % show mean dislocation/diffusion lines
    mean_dis = 10^mean(log10(A_dis))*xs.^mean(n).*normf_H2O^mean(r_dis)*exp(-(mean(Q_dis)+normP*mean(V_dis))/(normT*R));
    mean_dif = 10^mean(log10(A_dif))*xs.*normd^-mean(m).*normf_H2O^mean(r_dif)*exp(-(mean(Q_dif)+normP*mean(V_dif))/(normT*R));
    loglog(xs,mean_dif,'k-');
    loglog(xs,mean_dis,'k--');
    loglog(xs, 10^-2.97*xs.^(1.7).*normf_H2O.*normd^-0.51.*exp(-183000/(normT*R)),'c', 'Linewidth', 1);

%}

% Strain vs. temp
figure(2); hold off;
    for i=1:nout
        %model strain rate
        e_dis = A_dis(i)*norms^n(i).*normf_H2O^r_dis(i)*exp(-(Q_dis(i)+normP*V_dis(i))./(xt.*R));
        e_dif = A_dif(i)*norms*normd^-m(i).*normf_H2O^r_dif(i)*exp(-(Q_dif(i)+normP*V_dif(i))./(xt.*R));
        semilogy(1e4./xt, e_dis + e_dif,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('Strain rate [s^{-1}]')
            xlabel('10^4/T [K^{-1}]')
            title('Strain rate vs. Temperature')
            % ylim([1e-7,1e-1])
        end
    end
    
    for i=1:size(X,2)
        % original data
        semilogy(1e4./T{i}+(0.01*i),e{i}.*exp(-mean(X{i})),'go');
        % error bars
        [ex_0{i},ey_0{i}] = calc_err(T{i},dT{i},e{i},(de{i}./e{i}).*e{i});
        semilogy(1e4./ex_0{i}+(0.01*i),ey_0{i}.*exp(-mean(X{i})),'g-');

        % normalized data
        norm_factor = calc_norm("wet","T",out,file,norm_params, "V");
        semilogy(1e4./T{i}+(0.01*i),e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(T{i},dT{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        semilogy(1e4./ex{i}+(0.01*i),ey{i}.*exp(-mean(X{i})),'b-');
    end
    % show mean dislocation/diffusion lines
    mean_dis = 10^mean(log10(A_dis))*norms.^mean(n).*normf_H2O^mean(r_dis).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(xt.*R));
    mean_dif = 10^mean(log10(A_dif))*norms.*normd^-mean(m).*normf_H2O^mean(r_dif).*exp(-(mean(Q_dif)+normP*mean(V_dif))./(xt.*R));
    semilogy(1e4./xt,mean_dif,'k-');
    semilogy(1e4./xt,mean_dis,'k--');
    semilogy(1e4./xt, 10^-2.97*norms.^(1.7).*normf_H2O.*normd^-0.51.*exp(-183000./(xt.*R)),'c', 'Linewidth', 1);


% Strain vs. grain size

figure(3); hold off;
    for i=1:nout
        %model strain rate
        e_dis = A_dis(i)*norms^n(i).*normf_H2O^r_dis(i)*exp(-(Q_dis(i)+normP*V_dis(i))./(normT.*R));
        e_dif = A_dif(i)*norms.*xd.^-m(i).*normf_H2O^r_dif(i)*exp(-(Q_dif(i)+normP*V_dif(i))./(normT.*R));
        loglog(xd,e_dis + e_dif,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('Strain rate [s^{-1}]')
            xlabel('Grain size [um]')
            title('Strain rate vs. Grain size')
            ylim([1e-7,1e-2])
            xlim([1e0,1e2])
        end
    end
    
    for i=1:size(X,2)
        % original data
        loglog(d{i},e{i}.*exp(-mean(X{i})),'go');
        % error bars
        [ex_0{i},ey_0{i}] = calc_err(d{i},dd{i},e{i},(de{i}./e{i}).*e{i});
        loglog(ex_0{i},ey_0{i}.*exp(-mean(X{i})),'g-');

        % normalized data
        norm_factor = calc_norm("wet","d",out,file,norm_params, "V");
        loglog(d{i},e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(d{i},dd{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        loglog(ex{i},ey{i}.*exp(-mean(X{i})),'b-');
    end
    % get mean dislocation/diffusion lines
    mean_dis = 10^mean(log10(A_dis))*norms.^mean(n).*normf_H2O^mean(r_dis).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(normT*R));
    mean_dif = 10^mean(log10(A_dif))*norms.*xd.^-mean(m).*normf_H2O^mean(r_dif).*exp(-(mean(Q_dif)+normP*mean(V_dif))./(normT*R));
    loglog(xd,mean_dif,'k-');
    loglog(xd,ones(size(xd)).*mean_dis,'k--');
    loglog(xd, ones(size(xd)).*10^-2.97*norms^(1.7).*normf_H2O.*xd.^-0.51.*exp(-183000/(normT*R)),'c', 'Linewidth', 1);

    %}

% Strain vs. pressure
figure(4); hold off;
    for i=1:nout
        %model strain rate
        e_dis = A_dis(i)*norms^n(i).*normf_H2O^r_dis(i).*exp(-(Q_dis(i)+xP.*V_dis(i))./(normT.*R));
        e_dif = A_dif(i)*norms*normd^-m(i).*normf_H2O^r_dif(i).*exp(-(Q_dif(i)+xP.*V_dif(i))./(normT.*R));
        semilogy(xP./1e9,e_dis + e_dif,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('Strain rate [s^{-1}]')
            xlabel('Pressure [GPa]')
            title('Strain rate. vs. Pressure')
            ylim([1e-7,1e-2])
            xlim([0.8,1.8])
        end
    end
    
    for i=1:size(X,2)
        % original data
        semilogy(P{i}./1e9+(0.005*i),e{i}.*exp(-mean(X{i})),'go');
        % error bars
        [ex_0{i},ey_0{i}] = calc_err(P{i},dP{i},e{i},(de{i}./e{i}).*e{i});
        semilogy(ex_0{i}./1e9+(0.005*i),ey_0{i}.*exp(-mean(X{i})),'g-');

        % normalized data
        norm_factor = calc_norm("wet","P",out,file,norm_params, "V");
        semilogy(P{i}./1e9+(0.005*i),e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(P{i},dP{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        semilogy(ex{i}./1e9+(0.005*i),ey{i}.*exp(-mean(X{i})),'b-');
    end
    % get mean dislocation/diffusion lines
    mean_dis = 10^mean(log10(A_dis))*norms.^mean(n).*normf_H2O^mean(r_dis).*exp(-(mean(Q_dis)+xP.*mean(V_dis))./(normT*R));
    mean_dif = 10^mean(log10(A_dif))*norms.*normd^-mean(m).*normf_H2O^mean(r_dif).*exp(-(mean(Q_dif)+xP.*mean(V_dif))./(normT*R));
    semilogy(xP./1e9,mean_dif,'k-');
    semilogy(xP./1e9,mean_dis,'k--');
%}

% Strain vs. water fugacity
figure(5); hold off;
    for i=1:nout
        %model strain rate
        e_dis = A_dis(i)*norms^n(i).*xf.^r_dis(i).*exp(-(Q_dis(i)+normP*V_dis(i))./(normT.*R));
        e_dif = A_dif(i)*norms*normd^-m(i).*xf.^r_dif(i).*exp(-(Q_dif(i)+normP*V_dif(i))./(normT.*R));
        loglog(xf,e_dis + e_dif,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('Strain rate [s^{-1}]')
            xlabel('Water fugacity [MPa]')
            title('Strain rate. vs. Water fugacity')
            ylim([1e-7,1e-2])
            xlim([1e3,1e4])
        end
    end
    
    for i=1:size(X,2)
        % original data
        loglog(f_H2O{i}+(10*i),e{i}.*exp(-mean(X{i})),'go');
        % error bars
        [ex_0{i},ey_0{i}] = calc_err(f_H2O{i},df_H2O{i},e{i},(de{i}./e{i}).*e{i});
        loglog(ex_0{i}+(10*i),ey_0{i}.*exp(-mean(X{i})),'g-');

        % normalized data
        norm_factor = calc_norm("wet","f",out,file,norm_params, "V");
        loglog(f_H2O{i}+(10*i),e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(f_H2O{i},df_H2O{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        loglog(ex{i}+(10*i),ey{i}.*exp(-mean(X{i})),'b-');
    end
    % get mean dislocation/diffusion lines
    mean_dis = 10^mean(log10(A_dis))*norms.^mean(n).*xf.^mean(r_dis).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(normT*R));
    mean_dif = 10^mean(log10(A_dif))*norms.*normd^-mean(m).*xf.^mean(r_dif).*exp(-(mean(Q_dif)+normP*mean(V_dif))./(normT*R));
    loglog(xf,mean_dif,'k-');
    loglog(xf,mean_dis,'k--');
    loglog(xf, 10^-2.97*norms.^(1.7).*xf.*normd^-0.51.*exp(-183000/(normT*R)),'c', 'Linewidth', 1);

%}
%{
% Inter-run bias
figure(6); hold off;
    for i=1:size(X,2)
        plot(i,mean(X{i}),'bo');
        errorbar(i,mean(X{i}),std((X{i})),'b');
        if i==1
            ylabel('X_m')
            xlabel('m')
            xlim([1,size(X,2)+1])
            hold on;
        end
    end

%}

% Histograms
figure(7);
 subplot(4,3,1);
 hist(n,20); xlabel('n');
 subplot(4,3,2);
 hist(r_dis,20); xlabel('r_{dis}');
 subplot(4,3,3);
 hist(m,20); xlabel('m');
 subplot(4,3,4);
 hist(r_dif,20); xlabel('r_{dif}');
 subplot(4,3,5);
 hist(Q_dis/1e3,20); xlabel('Q_{dis} [kJ/mol]');
 subplot(4,3,6);
 hist(Q_dif/1e3,20); xlabel('Q_{dif} [kJ/mol]');
 subplot(4,3,7);
 hist(log10(A_dis),20); xlabel('log_{10}(A_{dis})');
 subplot(4,3,8);
 hist(log10(A_dif),20); xlabel('log_{10}(A_{dif})');
 subplot(4,3,9);
 hist(V_dis*1e6,20); xlabel('V_{dis} [m^3/mol]');
 subplot(4,3,10);
 hist(V_dif*1e6,20); xlabel('V_{dif} [m^3/mol]');
 subplot(4,3,11);
 hist(chi2,20); xlabel('\chi^2');

 %{
figure(8); hold off;
    subplot(1,2,1);
    for i=1:size(X,2)
        norm_factor = calc_norm("wet","s",out,file,norm_params);
        loglog(s{i},e{i}.*norm_factor{i},'-o');
        hold on;
    end
    ylabel('Strain rate [s^{-1}]')
    xlabel('Stress [MPa]')
    ylim([1e-6,1e-2])
    xlim([1e1,5e3])
    title("Without inter-run bias")

    subplot(1,2,2);
    for i=1:size(X,2)
        norm_factor = calc_norm("wet","s",out,file,norm_params);
        loglog(s{i},e{i}.*norm_factor{i}*exp(-mean(X{i})),'-o');
        hold on;
    end
    ylabel('Strain rate [s^{-1}]')
    xlabel('Stress [MPa]')
    ylim([1e-6,1e-2])
    xlim([1e1,5e3])
    title("With inter-run bias")

%}
xe = logspace(log10(min(e_all))-1,(log10(max(e_all)))+1,nx);
figure(9); hold off;
for (i = 1:length(X))
    e_check{i} = (10^mean(log10(A_dis)).*f_H2O{i}.^mean(r_dis).*s{i}.^mean(n).*exp(-(mean(Q_dis)+P{i}.*mean(V_dis))./(T{i}.*R)) + ...
    10^mean(log10(A_dif))*s{i}.*d{i}.^-mean(m).*f_H2O{i}.^mean(r_dif).*exp(-(mean(Q_dif)+P{i}.*mean(V_dif))./(T{i}.*R))...
    ).*exp(-mean(X{i}));
    loglog(e_check{i},e{i},'k*');
    if i ==1
        hold on; axis tight;
        ylabel('Observed strain rate [s^{-1}]')
        xlabel('Predicted strain rate [s^{-1}]')
    end
    errorbar(e_check{i},e{i}, de{i},'k*');
end
loglog(xe,xe,'k-');
%}
%}