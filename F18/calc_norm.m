% b: "wet" or "dry"
% c: parameter to be varied (all else normalized)
    %s = stress; T = temperature; d = grain size; P = pressure; f = water fugacity
% out: output file
% dat: data file
% norm_params: vector of normalized parameters
function norm = calc_norm(b,c,out,dat,norm_params, V)
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
        if (b == "wet")
            f_H2O{i} = data(:,11); % water fugacity [MPa]
            df_H2O{i} = data(:,12);
        end
    end
    %normalize data
    normT = norm_params(1);
    normP = norm_params(2);
    norms = norm_params(3);
    normd = norm_params(4);
    if exist('V','var')
        switch b
            case "wet"
                normf_H2O = norm_params(5);
                r_dif = out(:,4); m = out(:,5); Q_dif = out(:,6); V_dif = out(:,7);
                r_dis = out(:,8); n = out(:,9);Q_dis = out(:,10); V_dis = out(:,11);
                A_dif = out(:,end-1); A_dis = out(:,end);
                for i=1:dat(end,end)
                    e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*f_H2O{i}.^mean(r_dis).*exp(-(mean(Q_dis)+P{i}.*mean(V_dis))./(T{i}.*R)) + ...
                    mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*f_H2O{i}.^mean(r_dif).*exp(-(mean(Q_dif)+P{i}.*mean(V_dif))./(T{i}.*R));
                    switch c
                        case "s"
                            e_pred_T_d_P_f{i} = mean(log10(A_dis))*s{i}.^mean(n).*normf_H2O.^mean(r_dis).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(normT.*R)) + ...
                            mean(log10(A_dif))*s{i}.*normd^-mean(m).*normf_H2O.^mean(r_dif).*exp(-(mean(Q_dif)+normP*mean(V_dif))/(normT*R));
                            norm{i} = e_pred_T_d_P_f{i}./e_pred0{i};
                        case "T"
                            e_pred_s_d_P_f{i} = mean(log10(A_dis))*norms.^mean(n).*normf_H2O.^mean(r_dis).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(T{i}.*R)) + ...
                            mean(log10(A_dif))*norms.*normd^-mean(m).*normf_H2O.^mean(r_dif).*exp(-(mean(Q_dif)+normP*mean(V_dif))./(T{i}.*R));
                            norm{i} = e_pred_s_d_P_f{i}./e_pred0{i};
                        case "d"
                            e_pred_s_T_P_f{i} = mean(log10(A_dis))*norms.^mean(n).*normf_H2O.^mean(r_dis).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(normT.*R)) + ...
                            mean(log10(A_dif))*norms.*d{i}.^-mean(m).*normf_H2O.^mean(r_dif).*exp(-(mean(Q_dif)+normP*mean(V_dif))./(normT.*R));
                            norm{i} = e_pred_s_T_P_f{i}./e_pred0{i};
                        case "P"
                            e_pred_s_T_d_f{i} = mean(log10(A_dis))*norms.^mean(n).*normf_H2O.^mean(r_dis).*exp(-(mean(Q_dis)+P{i}.*mean(V_dis))./(normT.*R)) + ...
                            mean(log10(A_dif))*norms.*normd.^-mean(m).*normf_H2O.^mean(r_dif).*exp(-(mean(Q_dif)+P{i}.*mean(V_dif))./(normT.*R));
                            norm{i} = e_pred_s_T_d_f{i}./e_pred0{i};
                        case "f"
                            e_pred_s_T_d_P{i} = mean(log10(A_dis))*norms.^mean(n).*f_H2O{i}.^mean(r_dis).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(normT.*R)) + ...
                            mean(log10(A_dif))*norms.*normd.^-mean(m).*f_H2O{i}.^mean(r_dif).*exp(-(mean(Q_dif)+normP*mean(V_dif))./(normT.*R));
                            norm{i} = e_pred_s_T_d_P{i}./e_pred0{i};
                    end
                end

            case "dry"
                m = out(:,4); Q_dif = out(:,5); V_dif = out(:,6);
                n = out(:,7); Q_dis = out(:,8); V_dis = out(:,9);
                A_dif = out(:,end-1); A_dis = out(:,end);
                for i=1:dat(end,end)
                    e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*exp(-(mean(Q_dis)+P{i}.*mean(V_dis))./(T{i}.*R)) + ...
                    mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*exp(-(mean(Q_dif)+P{i}.*mean(V_dif))./(T{i}.*R));
                    switch c
                        case "s"
                            e_pred_T_d_P{i} = mean(log10(A_dis))*s{i}.^mean(n).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(normT.*R)) + ...
                            mean(log10(A_dif))*s{i}.*normd^-mean(m).*exp(-(mean(Q_dif)+normP*mean(V_dif))/(normT*R));
                            norm{i} = e_pred_T_d_P{i}./e_pred0{i};
                        case "T"
                            e_pred_s_d_P{i} = mean(log10(A_dis))*norms.^mean(n).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(T{i}.*R)) + ...
                            mean(log10(A_dif))*norms.*normd^-mean(m).*exp(-(mean(Q_dif)+normP*mean(V_dif))./(T{i}.*R));
                            norm{i} = e_pred_s_d_P{i}./e_pred0{i};
                        case "d"
                            e_pred_s_T_P{i} = mean(log10(A_dis))*norms.^mean(n).*exp(-(mean(Q_dis)+normP*mean(V_dis))./(normT.*R)) + ...
                            mean(log10(A_dif))*norms.*d{i}.^-mean(m).*exp(-(mean(Q_dif)+normP*mean(V_dif))./(normT.*R));
                            norm{i} = e_pred_s_T_P{i}./e_pred0{i};
                        case "P"
                            e_pred_s_T_d{i} = mean(log10(A_dis))*norms.^mean(n).*exp(-(mean(Q_dis)+P{i}.*mean(V_dis))./(normT.*R)) + ...
                            mean(log10(A_dif))*norms.*normd.^-mean(m).*exp(-(mean(Q_dif)+P{i}.*mean(V_dif))./(normT.*R));
                            norm{i} = e_pred_s_T_d{i}./e_pred0{i};
                    end
                end
    end
    else
        switch b
        case "wet"
            normf_H2O = norm_params(5);
            r_dif = out(:,4); m = out(:,5); Q_dif = out(:,6);
            r_dis = out(:,7); n = out(:,8);Q_dis = out(:,9);
            A_dif = out(:,end-1); A_dis = out(:,end);
            for i=1:dat(end,end)
                e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*f_H2O{i}.^mean(r_dis).*exp(-mean(Q_dis)./(T{i}.*R)) + ...
                mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*f_H2O{i}.^mean(r_dif).*exp(-mean(Q_dif)./(T{i}.*R));
                switch c
                    case "s"
                        e_pred_T_d_P_f{i} = mean(log10(A_dis))*s{i}.^mean(n).*normf_H2O.^mean(r_dis).*exp(-mean(Q_dis)./(normT.*R)) + ...
                        mean(log10(A_dif))*s{i}.*normd^-mean(m).*normf_H2O.^mean(r_dif).*exp(-mean(Q_dif)/(normT*R));
                        norm{i} = e_pred_T_d_P_f{i}./e_pred0{i};
                    case "T"
                        e_pred_s_d_P_f{i} = mean(log10(A_dis))*norms.^mean(n).*normf_H2O.^mean(r_dis).*exp(-(mean(Q_dis))./(T{i}.*R)) + ...
                        mean(log10(A_dif))*norms.*normd^-mean(m).*normf_H2O.^mean(r_dif).*exp(-mean(Q_dif)./(T{i}.*R));
                        norm{i} = e_pred_s_d_P_f{i}./e_pred0{i};
                    case "d"
                        e_pred_s_T_P_f{i} = mean(log10(A_dis))*norms.^mean(n).*normf_H2O.^mean(r_dis).*exp(-(mean(Q_dis))./(normT.*R)) + ...
                        mean(log10(A_dif))*norms.*d{i}.^-mean(m).*normf_H2O.^mean(r_dif).*exp(-(mean(Q_dif))./(normT.*R));
                        norm{i} = e_pred_s_T_P_f{i}./e_pred0{i};
                    case "P"
                        e_pred_s_T_d_f{i} = mean(log10(A_dis))*norms.^mean(n).*normf_H2O.^mean(r_dis).*exp(-(mean(Q_dis))./(normT.*R)) + ...
                        mean(log10(A_dif))*norms.*normd.^-mean(m).*normf_H2O.^mean(r_dif).*exp(-(mean(Q_dif))./(normT.*R));
                        norm{i} = e_pred_s_T_d_f{i}./e_pred0{i};
                    case "f"
                        e_pred_s_T_d_P{i} = mean(log10(A_dis))*norms.^mean(n).*f_H2O{i}.^mean(r_dis).*exp(-(mean(Q_dis))./(normT.*R)) + ...
                        mean(log10(A_dif))*norms.*normd.^-mean(m).*f_H2O{i}.^mean(r_dif).*exp(-(mean(Q_dif))./(normT.*R));
                        norm{i} = e_pred_s_T_d_P{i}./e_pred0{i};
                end
            end

        case "dry"
            m = out(:,4); Q_dif = out(:,5);
            n = out(:,6); Q_dis = out(:,7);
            A_dif = out(:,end-1); A_dis = out(:,end);
            for i=1:dat(end,end)
                e_pred0{i} = mean(log10(A_dis)).*s{i}.^mean(n).*exp(-(mean(Q_dis))./(T{i}.*R)) + ...
                mean(log10(A_dif)).*s{i}.*d{i}.^-mean(m).*exp(-(mean(Q_dif))./(T{i}.*R));
                switch c
                    case "s"
                        e_pred_T_d_P{i} = mean(log10(A_dis))*s{i}.^mean(n).*exp(-(mean(Q_dis))./(normT.*R)) + ...
                        mean(log10(A_dif))*s{i}.*normd^-mean(m).*exp(-(mean(Q_dif))/(normT*R));
                        norm{i} = e_pred_T_d_P{i}./e_pred0{i};
                    case "T"
                        e_pred_s_d_P{i} = mean(log10(A_dis))*norms.^mean(n).*exp(-(mean(Q_dis))./(T{i}.*R)) + ...
                        mean(log10(A_dif))*norms.*normd^-mean(m).*exp(-(mean(Q_dif))./(T{i}.*R));
                        norm{i} = e_pred_s_d_P{i}./e_pred0{i};
                    case "d"
                        e_pred_s_T_P{i} = mean(log10(A_dis))*norms.^mean(n).*exp(-(mean(Q_dis))./(normT.*R)) + ...
                        mean(log10(A_dif))*norms.*d{i}.^-mean(m).*exp(-(mean(Q_dif))./(normT.*R));
                        norm{i} = e_pred_s_T_P{i}./e_pred0{i};
                    case "P"
                        e_pred_s_T_d{i} = mean(log10(A_dis))*norms.^mean(n).*exp(-(mean(Q_dis))./(normT.*R)) + ...
                        mean(log10(A_dif))*norms.*normd.^-mean(m).*exp(-(mean(Q_dif))./(normT.*R));
                        norm{i} = e_pred_s_T_d{i}./e_pred0{i};
                end
            end
    end

end 