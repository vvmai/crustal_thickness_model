function IQR = calc_IQR(file, flow_law, water, refwater, T, V)
    % file: .out file
    % flow_law: flow law (GBS, composite)
    % T: temperature
    % water: "wet" or "dry" inversion
    % refwater: COH/fH2O value to use as reference
    % V: include action volume? (optional)
    refE = 1e-15;
    R = 8.3145;
    if exist('V','var')
        switch flow_law
        case "GBS"
        case "composite"
            switch water
            case "wet"
            case "dry"
                refP = 1.5e9;
                out = file;
                m = out(:,4); Q_dif = out(:,5); V_dif = out(:,6); n = out(:,7); Q_dis = out(:,8);
                V_dis = out(:,9); A_dif = out(:,end-1); A_dis = out(:,end);
                for i=1:length(T)
                    for j=1:length(out)
                        s{i}(j) = (refE./(A_dis(j)*exp(-(Q_dis(j)+V_dis(j)*refP)./(R.*T(i))))).^(1/n(j));
                    end
                    IQR_difference = iqr(log10(s{i}));
                    IQR(i,1) = median(log10(s{i}))-IQR_difference/2;
                    IQR(i,2) = median(log10(s{i}));
                    IQR(i,3) = median(log10(s{i}))+IQR_difference/2;
                    %{
                    s_pd{i} = makedist('Normal','mu', mean(log10(s{i})), 'sigma', std(log10(s{i})));
                    xs{i} = linspace(min(log10(s{i})), max(log10(s{i})), 100);
                    s_cdf{i} = cdf(s_pd{i}, xs{i});
                    IQR(i,1) = interp1(s_cdf{i}, xs{i}, 0.25);
                    IQR(i,2) = interp1(s_cdf{i}, xs{i}, 0.50);
                    IQR(i,3) = interp1(s_cdf{i}, xs{i}, 0.75);
                    %}
                end
            end
        case "dislocation"
        case "diffusion"
    end
    else
        switch flow_law
            case "GBS"
                out = file;
                refF = refwater;
                refd = water;
                m = out(:,4); n = out(:,5); r = out(:,6);
                Q = out(:,7); A = out(:,end);
                for i=1:length(T)
                    for j=1:length(out)
                        s{i}(j) = (refE./(A(j).*refd^m(j).*refF.^(r(j)).*exp(-Q(j)./(R.*T(i))))).^(1/n(j));
                    end
                    IQR_difference = iqr(log10(s{i}));
                    IQR(i,1) = median(log10(s{i}))-IQR_difference/2;
                    IQR(i,2) = median(log10(s{i}));
                    IQR(i,3) = median(log10(s{i}))+IQR_difference/2;
                end
            case "composite"
                switch water
                case "wet"
                    refCOH = refwater;
                    out = file;
                    r_dif = out(:,4); m = out(:,5); Q_dif = out(:,6);
                    r_dis = out(:,7); n = out(:,8); Q_dis = out(:,9);
                    A_dif = out(:,end-1); A_dis = out(:,end);
                    for i=1:length(T)
                        for j=1:length(out)
                            s{i}(j) = (refE./(A_dis(j).*refCOH.^(r_dis(j)).*exp(-Q_dis(j)./(R.*T(i))))).^(1/n(j));
                        end
                        IQR_difference = iqr(log10(s{i}));
                        IQR(i,1) = median(log10(s{i}))-IQR_difference/2;
                        IQR(i,2) = median(log10(s{i}));
                        IQR(i,3) = median(log10(s{i}))+IQR_difference/2;
                    end
                case "dry"
                    out = file;
                    m = out(:,4); Q_dif = out(:,5); n = out(:,6); Q_dis = out(:,7);
                    A_dif = out(:,end-1); A_dis = out(:,end);
                    for i=1:length(T)
                        for j=1:length(out)
                            s{i}(j) = (refE./(A_dis(j)*exp(-Q_dis(j)./(R.*T(i))))).^(1/n(j));
                
                        end
                        IQR_difference = iqr(log10(s{i}));
                        IQR(i,1) = median(log10(s{i}))-IQR_difference/2;
                        IQR(i,2) = median(log10(s{i}));
                        IQR(i,3) = median(log10(s{i}))+IQR_difference/2;
                    end
                end
            case "dislocation"
            case "diffusion"
                switch water
                case "wet"
                case "dry"
                    out = file;
                    refd = refwater;
                    m = out(:,4); Q_dif = out(:,5); n = out(:,6); Q_dis = out(:,7);
                    A_dif = out(:,10); A_dis = out(:,11);
                    for i=1:length(T)
                        for j=1:length(out)
                            s{i}(j) = (refE./(A_dif(j)*refd^-m(j)*exp(-Q_dif(j)./(R.*T(i)))));
                        end
                        IQR_difference = iqr(log10(s{i}));
                        IQR(i,1) = median(log10(s{i}))-IQR_difference/2;
                        IQR(i,2) = median(log10(s{i}));
                        IQR(i,3) = median(log10(s{i}))+IQR_difference/2;
                    end
                end
            end
        end
    end
    