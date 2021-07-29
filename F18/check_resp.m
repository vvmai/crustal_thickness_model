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
    output = load(['./' run_id '_d_fixed.out']);
    % skip first 1000 runs, then select results from every 10 runs
    %out = output(100:1:end,1:end);
    out = output;
    nout = length(out);
    chi2 = out(:,3);
    m = out(:,4); r_dif = out(:,5); Q_dif = out(:,6);
    r_dis = out(:,7); n = out(:,8); Q_dis = out(:,9);
    A_dif = out(:,end-1); A_dis = out(:,end);
    % 10-22 = Inter-run bias relative to sample 1
    for i=10:(size(out,2)-2)
        X{i-9} = out(:,i);
    end


A = [log10(A_dis) log10(A_dif) n m r_dis r_dif Q_dis Q_dif];
[Corr] = corrcoef(A)

% Print summary data
    disp(['id=' run_id]);
    disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(1*std(chi2))]);

    disp(['n = ' num2str(mean(n)) ' +/- ' num2str(1*std(n))]);
    disp(['r_dis = ' num2str(mean(r_dis)) ' +/- ' num2str(1*std(r_dis))]);
    disp(['Q_dis = ' num2str(mean(Q_dis)/1e3) ' +/- ' num2str(1*std(Q_dis)/1e3) ' kJ/mol']);
    disp(['log10(A_dis) = ' num2str(mean(log10(A_dis))) ' +/- ' num2str(1*std(log10(A_dis)))]);

    disp(['m = ' num2str(mean(m)) ' +/- ' num2str(1*std(m))]);
    disp(['r_dif = ' num2str(mean(r_dif)) ' +/- ' num2str(1*std(r_dif))]);
    disp(['Q_dif = ' num2str(mean(Q_dif)/1e3) ' +/- ' num2str(1*std(Q_dif)/1e3) ' kJ/mol']);
    disp(['log10(A_dif) = ' num2str(mean(log10(A_dif))) ' +/- ' num2str(1*std(log10(A_dif)))]);
    for i=2:length(X)
        disp(['X_' num2str(i) ' = ' num2str(mean(X{i})) ' +/- ' num2str(1*std(X{i}))]);
    end

% set constants
    R = 8.3145;
    nx = 30;
    %stress data
    xs = logspace(log10(min(s_all))-1,log10(max(s_all))+2,nx);
    %temperature data
    xt = linspace(min(T_all)-200,max(T_all)+500,nx);
    %grain size data
    xd = logspace(log10(min(d_all))-1,log10(max(d_all))+2,nx);
    %pressure data
    xP = linspace(min(P_all)-0.25*1e9,max(P_all)+0.5*1e9,nx);
    %water fugacity data
    xf = logspace(log10(min(f_H2O_all))-0.5,log10(max(f_H2O_all))+0.5,nx);
    %normalize data
    normT = 1073; % Temperature [K]
    norms = 250; % Stress [MPa]
    normP = 1.5*1e9; % Pressure [Pa]
    normd = 10; % Grain size [um]
    normf_H2O = 4.5*1e3; % Water fugacity [MPa]
    norm_params = [normT, normP, norms, normd, normf_H2O];

mat = [log10(A_dis) log10(A_dif) n m r_dis r_dif Q_dis Q_dif];
[cc] = corrcoef(mat)

% Strain vs. stress
figure(1); hold off;
    for i=1:nout
        %model strain rate
        e_dis = A_dis(i)*xs.^n(i).*normf_H2O^r_dis(i).*exp(-(Q_dis(i))./(normT*R));
        e_dif = A_dif(i)*xs.*normd^-m(i).*normf_H2O^r_dif(i).*exp(-(Q_dif(i))./(normT*R));
        loglog(xs, e_dis + e_dif,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('strain rate [s^{-1}]')
            xlabel('stress [MPa]')
            title('F18 (composite)')
            xlim([10^(1),10^(3.5)])
            ylim([1e-7,1e0])
            text(12, 10^(-1.7),['T = ' num2str(normT) ' K' newline 'd = ' num2str(normd) '\mum' newline 'f_{H_2O} = ' num2str(normf_H2O) ' MPa'], 'FontSize', 20)
        end
    end

    for i=1:size(X,2)
        % normalize
        norm_factor = calc_norm("wet","s",out,file,norm_params);
        loglog(s{i},e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(s{i},ds{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        loglog(ex{i},ey{i}.*exp(-mean(X{i})),'b-');
        % original data
        %loglog(s{i},e{i}.*exp(-mean(X{i})),'go');
        %[ex_0{i},ey_0{i}] = calc_err(s{i},ds{i},e{i},(de{i}./e{i}).*e{i});
        %loglog(ex_0{i},ey_0{i}.*exp(-mean(X{i})),'g-');
    end
    % show mean dislocation/diffusion lines
    mean_dis = 10^mean(log10(A_dis))*xs.^mean(n).*normf_H2O^mean(r_dis)*exp(-mean(Q_dis)/(normT*R));
    mean_dif = 10^mean(log10(A_dif))*xs.*normd^-mean(m).*normf_H2O^mean(r_dif)*exp(-mean(Q_dif)/(normT*R));
    loglog(xs, mean_dif,'k-', 'Linewidth', 1.5);
    loglog(xs, mean_dis,'k--', 'Linewidth', 1.5);
loglog(xs, 10^-2.97*xs.^(1.7).*normf_H2O.*normd^-0.51.*exp(-183000/(normT*R)),'c', 'Linewidth', 1.5);
set(gca, 'LineWidth', 1)
set(gca,'FontSize',20)
%}

fprintf("for e = 1e-15, T = %d, K\nCalculated: %d\nFukuda: %d\n", normT, (1e-15./(10^mean(log10(A_dis))*exp(-mean(Q_dis)./(R.*normT)))).^(1/mean(n)), (1e-15./(10^-2.97.*normd^-0.51*normf_H2O.*exp(-183000./(R.*normT)))).^(1/1.7))

% Strain vs. temp
figure(2); hold off;
    for i=1:nout
        %model strain rate
        e_dis = A_dis(i).*norms^n(i).*normf_H2O^r_dis(i).*exp(-(Q_dis(i))./(xt.*R));
        e_dif = A_dif(i).*norms*normd^-m(i).*normf_H2O^r_dif(i).*exp(-(Q_dif(i))./(xt.*R));
        semilogy(1e4./xt, e_dis + e_dif,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('strain rate [s^{-1}]')
            xlabel('10^4/T [K^{-1}]')
            %xlim([7 13])
            %ylim([1e-7,1e0])
            text(9.77, 10^(-0.75),['\sigma = ' num2str(norms) ' MPa' newline 'd = ' num2str(normd) '\mum' newline 'f_{H_2O} = ' num2str(normf_H2O) ' MPa'], 'FontSize', 20)
        end
    end
    
    for i=1:size(X,2)
        % original data
        %semilogy(1e4./T{i}+(0.01*i),e{i}.*exp(-mean(X{i})),'go');
        % error bars
        %[ex_0{i},ey_0{i}] = calc_err(T{i},dT{i},e{i},(de{i}./e{i}).*e{i});
        %semilogy(1e4./ex_0{i}+(0.01*i),ey_0{i}.*exp(-mean(X{i})),'g-');

        % normalized data
        norm_factor = calc_norm("wet","T",out,file,norm_params);
        semilogy(1e4./T{i}+(0.01*i),e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(T{i},dT{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        semilogy(1e4./ex{i}+(0.01*i),ey{i}.*exp(-mean(X{i})),'b-');
    end
    % show mean dislocation/diffusion lines
    mean_dis = 10^mean(log10(A_dis))*norms.^mean(n).*normf_H2O^mean(r_dis).*exp(-(mean(Q_dis))./(xt.*R));
    mean_dif = 10^mean(log10(A_dif))*norms.*normd^-mean(m).*normf_H2O^mean(r_dif).*exp(-(mean(Q_dif))./(xt.*R));
    semilogy(1e4./xt,mean_dif,'k-', 'Linewidth', 1.5);
    semilogy(1e4./xt,mean_dis,'k--', 'Linewidth', 1.5);
    semilogy(1e4./xt, 10^-2.97*norms.^(1.7).*normf_H2O.*normd^-0.51.*exp(-183000./(xt.*R)),'c', 'Linewidth', 1.5);
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',20)

% Strain vs. grain size
figure(3); hold off;
    for i=1:nout
        %model strain rate
        e_dis = A_dis(i).*norms^n(i).*normf_H2O^r_dis(i).*exp(-(Q_dis(i))./(normT.*R));
        e_dif = A_dif(i).*norms.*xd.^-m(i).*normf_H2O^r_dif(i).*exp(-(Q_dif(i))./(normT.*R));
        loglog(xd, e_dis + e_dif,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('strain rate [s^{-1}]')
            xlabel('grain size [\mum]')
            ylim([1e-7,10^(-0.5)])
            xlim([1e0,1e3])
            text(25, 10^(-2),['T = ' num2str(normT) ' K' newline '\sigma = ' num2str(norms) ' MPa' newline 'f_{H_2O} = ' num2str(normf_H2O) ' MPa'], 'FontSize', 20)
        end
    end
    
    for i=1:size(X,2)
        % original data
        %loglog(d{i},e{i}.*exp(-mean(X{i})),'go');
        % error bars
        %[ex_0{i},ey_0{i}] = calc_err(d{i},dd{i},e{i},(de{i}./e{i}).*e{i});
        %loglog(ex_0{i},ey_0{i}.*exp(-mean(X{i})),'g-');

        % normalized data
        norm_factor = calc_norm("wet","d",out,file,norm_params);
        loglog(d{i},e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(d{i},dd{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        loglog(ex{i},ey{i}.*exp(-mean(X{i})),'b-');
    end
    % get mean dislocation/diffusion lines
    mean_dis = 10^mean(log10(A_dis))*norms.^mean(n).*normf_H2O^mean(r_dis).*exp(-(mean(Q_dis))./(normT*R));
    mean_dif = 10^mean(log10(A_dif))*norms.*xd.^-mean(m).*normf_H2O^mean(r_dif).*exp(-(mean(Q_dif))./(normT*R));
    loglog(xd,mean_dif,'k-', 'Linewidth', 1.5);
    loglog(xd,ones(size(xd)).*mean_dis,'k--', 'Linewidth', 1.5);
    loglog(xd, ones(size(xd)).*10^-2.97*norms^(1.7).*normf_H2O.*xd.^-0.51.*exp(-183000/(normT*R)),'c', 'Linewidth', 1.5);
    set(gca, 'LineWidth', 1)
set(gca,'FontSize',20)
%}

% Strain vs. water fugacity
figure(4); hold off;
    for i=1:nout
        %model strain rate
        e_dis = A_dis(i)*norms^n(i).*xf.^r_dis(i).*exp(-(Q_dis(i))./(normT.*R));
        e_dif = A_dif(i)*norms*normd^-m(i).*xf.^r_dif(i).*exp(-(Q_dif(i))./(normT.*R));
        loglog(xf,e_dis + e_dif,'r:');
        if i==1
            hold on; box on; axis tight;
            ylabel('strain rate [s^{-1}]')
            xlabel('water fugacity [MPa]')
            set(gca, 'LineWidth', 1)
            set(gca,'FontSize',20)
            ylim([1e-7 1e-1])
            text(10^(2.65), 10^(-2.25), ['\sigma = ' num2str(norms) ' MPa' newline 'T = ' num2str(normT) ' K' newline 'd = ' num2str(normd) ' \mu m'], 'FontSize', 20)
        end
    end
    
    for i=1:size(X,2)
        % original data
        %loglog(f_H2O{i}+(10*i),e{i}.*exp(-mean(X{i})),'go');
        % error bars
        %[ex_0{i},ey_0{i}] = calc_err(f_H2O{i},df_H2O{i},e{i},(de{i}./e{i}).*e{i});
        %loglog(ex_0{i}+(10*i),ey_0{i}.*exp(-mean(X{i})),'g-');

        % normalized data
        norm_factor = calc_norm("wet","f",out,file,norm_params);
        loglog(f_H2O{i}+(10*i),e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
        % error bars
        [ex{i},ey{i}] = calc_err(f_H2O{i},df_H2O{i},e{i}.*norm_factor{i},(de{i}./e{i}).*e{i}.*norm_factor{i});
        loglog(ex{i}+(10*i),ey{i}.*exp(-mean(X{i})),'b-');
    end
    % get mean dislocation/diffusion lines
    mean_dis = 10^mean(log10(A_dis))*norms.^mean(n).*xf.^mean(r_dis).*exp(-(mean(Q_dis))./(normT*R));
    mean_dif = 10^mean(log10(A_dif))*norms.*normd^-mean(m).*xf.^mean(r_dif).*exp(-(mean(Q_dif))./(normT*R));
    loglog(xf,mean_dif,'k-', 'Linewidth', 1.5);
    loglog(xf,mean_dis,'k--', 'Linewidth', 1.5);
    loglog(xf, 10^-2.97*norms.^(1.7).*xf.*normd^-0.51.*exp(-183000/(normT*R)),'c', 'Linewidth', 1.5);
%}


% Histograms
figure(5);
 subplot(3,3,1);
 hist(n,20); xlabel('n');
 subplot(3,3,2);
 hist(r_dis,20); xlabel('r_{dis}');
 subplot(3,3,3);
 hist(m,20); xlabel('m');
 subplot(3,3,4);
 hist(r_dif,20); xlabel('r_{dif}');
 subplot(3,3,5);
 hist(Q_dis/1e3,20); xlabel('Q_{dis} [kJ/mol]');
 subplot(3,3,6);
 hist(Q_dif/1e3,20); xlabel('Q_{dif} [kJ/mol]');
 subplot(3,3,7);
 hist(log10(A_dis),20); xlabel('log_{10}(A_{dis})');
 subplot(3,3,8);
 hist(log10(A_dif),20); xlabel('log_{10}(A_{dif})');
 subplot(3,3,9);
 hist(chi2,20); xlabel('\chi^2');
 %}
 %}