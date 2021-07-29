function norm = calc_norm(c,out,dat,norm_params)
    % dat = data file
    R = 8.3145;
    for i=1:dat(end,end)
        data = dat(dat(:,end) == i,:);
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
    %normalize data
    normT = norm_params(1);
    normP = norm_params(2);
    norms = norm_params(3);
    normd = norm_params(4);
    chi2 = out(:,3); m = out(:,4); Q_dif = out(:,5); V_dif = out(:,6); n = out(:,7); Q_dis = out(:,8); V_dis = out(:,9);
    A_dif = out(:,end-1); A_dis = out(:,end);
    for i=1:dat(end,end)
        e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*exp(-(mean(Q_dis)+P{i}.*mean(V_dis))./(T{i}.*R)) + ...
        mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*exp(-(mean(Q_dif)+P{i}.*mean(V_dif))./(T{i}.*R));
        if c == "s"
            e_pred_T_d_P{i} = mean(log10(A_dis))*s{i}.^mean(n).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(normT.*R)) + ...
            mean(log10(A_dif))*s{i}.*normd^-mean(m).*exp(-(mean(Q_dif)+normP*mean(V_dif))/(normT*R));
            norm{i} = e_pred_T_d_P{i}./e_pred0{i};

        elseif c == "T"
            e_pred_s_d_P{i} = mean(log10(A_dis))*norms.^mean(n).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(T{i}.*R)) + ...
            mean(log10(A_dif))*norms.*normd^-mean(m).*exp(-(mean(Q_dif)+normP*mean(V_dif))./(T{i}.*R));
            norm{i} = e_pred_s_d_P{i}./e_pred0{i};

        elseif c == "d"
            e_pred_s_T_P{i} = mean(log10(A_dis))*norms.^mean(n).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(normT.*R)) + ...
            mean(log10(A_dif))*norms.*d{i}.^-mean(m).*exp(-(mean(Q_dif)+normP*mean(V_dif))./(normT.*R));
            norm{i} = e_pred_s_T_P{i}./e_pred0{i};

        elseif c == "P"
            e_pred_s_T_d{i} = mean(log10(A_dis))*norms.^mean(n).*exp(-(mean(Q_dis)+P{i}.*mean(V_dis))./(normT.*R)) + ...
            mean(log10(A_dif))*norms.*normd.^-mean(m).*exp(-(mean(Q_dif)+P{i}.*mean(V_dif))./(normT.*R));
            norm{i} = e_pred_s_T_d{i}./e_pred0{i};
        end
    end
    