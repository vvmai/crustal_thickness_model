clear
% Geotherm
    % Depth from 0 to 40 km, in increments of 1 km
    z = linspace(0,120000,100);
    % Heat flow equation
    k = 2.7; % Thermal conductivity [W/mK] (assume constant)
    D = 10000; % Depth of upper crust producting radioactive heat [m]
    T_0 = 0; % surface temperature [C]
    Q_0 = 0.046; % Surface heat flow [W/m2]
    Q_r = 0.025; % Nonradiogenic heat flow [W/m2]
    A_0 = (Q_0 - Q_r)/D; % Surface heat production [W/m3]
    T = zeros(1, length(z));
    Q = zeros(1, length(z));
    Q(1) = Q_r;
    T(1) = T_0;
    for i=2:length(z)
        % Exponentially decreasing rate of radiogenic heat production
        A(i) = A_0*exp(-z(i)/D);
        % Heat flow as a function of depth
        Q(i) = Q(1) - A(i)*D;
        % Temperature as a function of depth
        T(i) = T(1) + (Q(1)/k)*z(i) + (A_0*D^2)/(k) - (A(i).*D^2)/(k); 
    end
    %{
    % Plot geotherm
    figure(1);
    set(gca,'xaxislocation','top'); hold on;
    plot(T,-z/1000);
    plot(400*ones(length(z)),-z/1000,'k--');
    plot(200*ones(length(z)),-z/1000,'k--');
    xlabel('Temperature (C)')
    ylabel('Depth (km)')
    %}
%}

%constants
    R = 8.3145;
    refd = 1000; % Reference grain size [um]
    refE = 1e-15; % Reference strain rate [s-1]
    refF = 300; % Reference water fugacity [MPa]
    refP = 1.5e9; % Reference pressure [Pa]
    refCOH = 4000;
    T = T + 273.15; % Switch to K

%Flow laws from papers
    % RB04a:
    %Q = 242000 kJ/mol; A = 10^-4.93; n = 2.97;
    s_dis_RB = (refE./(exp(-11.35)*refF*exp(-242000./(R.*T)))).^(1/2.97);

    % GT95:
    %n = 4.0 +/- 0.9;
    %Q = 223 +/- 56 kJ/mol (without melt); A = 1.1e(-4 +/- 2) (without melt)
    s_GT = (refE./(exp(-7.58).*exp(-223000./(R.*T)))).^(1/4);

    % LP92:
    %n = 4.0 +/- 0.8; Q = 152 +/- 71 kJ/mol; A = 1.2e-08
    s_LP = (refE./(1.2e-8.*exp(-152000./(R.*T)))).^(1/4);

    % Fukuda18:
    % include water fugacity/grain size??
    %A = 10^-2.97 +/- 0.23; n = 1.7 +/- 0.2; m = -0.51 +/- 0.13; r = 1.0 +/- 0.2; Q = 183 +/- 25 kJ/mol
    s_fukuda = (refE./(10^-2.97.*refd^-0.51*refF.*exp(-183000./(R.*T)))).^(1/1.7);

    % Renner02:
    % Q = 200 kJ/mol; m = 0.5; A = 10; K = 115; Sigma = 7.8
    T_m = 1600; % calcite melting T [K]
    s = logspace(-5,6,100);
    A = 10;
    K = 115;
    for i=1:100
        e{i} = A*s.^2.*exp(s./((7.8 + K*refd^-0.5)*(T_m - T(i)))).*exp(-200000./(R.*T(i)));
        s_renner(i) = exp(interp1(log(e{i}),log(s),log(refE)));
    end

% Inversion flow laws
%RB04ab
% composite inversion, use only dislocation part
output_RB = load('./rutter04ab_X.out');
out_RB = output_RB(100:2:end,1:end);
m_RB = out_RB(:,4); Q_dif_RB = out_RB(:,5); n_RB = out_RB(:,6); Q_dis_RB = out_RB(:,7);
A_dif_RB = log10(out_RB(:,10)); A_dis_RB = log10(out_RB(:,11));
s_RB_MCMC = (refE./(10^mean(A_dis_RB)*exp(-mean(Q_dis_RB)./(R.*T)))).^(1/mean(n_RB));
s_RB_MCMC_dif = (refE./(10^mean(A_dif_RB).*refd^mean(-m_RB).*exp(-mean(Q_dif_RB)./(R.*T))));
IQR_RB = calc_IQR(out_RB, "composite", "dry", 0, T);
IQR_RB_dif = calc_IQR(out_RB, "diffusion", "dry", refd, T);

%GT95
% composite inversion, using dislocation part, no avtivation volume
% (dislocation-only inversion did not use corrected stress)
output_GT = load('./GT95_d_fixed.out');
out_GT = output_GT(100:2:end,1:end);
m_GT = out_GT(:,4); Q_dif_GT = out_GT(:,5); 
%V_dif_GT = out_GT(:,6);
n_GT = out_GT(:,6); Q_dis_GT = out_GT(:,7);
%V_dis_GT = out_GT(:,9);
A_dif_GT = out_GT(:,end-1); A_dis_GT = out_GT(:,end);
%s_GT_MCMC = (refE./(10^mean(log10(A_dis_GT))*exp(-(mean(Q_dis_GT)+mean(V_dis_GT)*refP)./(T.*R)))).^(1/mean(n_GT));
s_GT_MCMC = (refE./(10^mean(log10(A_dis_GT))*exp(-(mean(Q_dis_GT))./(T.*R)))).^(1/mean(n_GT));
IQR_GT = calc_IQR(out_GT, "composite", "dry", 0, T);
%LP92
% composite inversion, no activation volume, use only dislocation part
output_LP = load('./LP92_pXC.out');
%output_LP = load('./LP92_COH.out');
out_LP = output_LP(100:2:end,1:end);
r_dif_LP = out_LP(:,4); m_LP = out_LP(:,5); Q_dif_LP = out_LP(:,6);
r_dis_LP = out_LP(:,7); n_LP = out_LP(:,8); Q_dis_LP = out_LP(:,9);
A_dif_LP = out_LP(:,24); A_dis_LP = out_LP(:,25);

s_LP_MCMC = (refE./(10^mean(log10(A_dis_LP)).*refCOH.^mean(r_dis_LP).*exp(-mean(Q_dis_LP)./(R.*T)))).^(1/mean(n_LP));
IQR_LP = calc_IQR(out_LP, "composite", "wet", 4000, T);


%Fukuda18
% composite inversion, no activation volume, use only dislocation part
% COH: assumes samples are saturated
output_fukuda = load('./fukuda18_d_fixed.out');
out_fukuda = output_fukuda(100:2:end,1:end);
r_dis_fukuda = out_fukuda(:,7);
n_fukuda = out_fukuda(:,8); Q_dis_fukuda = out_fukuda(:,9); A_dis_fukuda = out_fukuda(:,end);
s_fukuda_MCMC = (refE./(10^mean(log10(A_dis_fukuda))*refF^mean(r_dis_fukuda)*exp(-mean(Q_dis_fukuda)./(R.*T)))).^(1/mean(n_fukuda));
IQR_fukuda = calc_IQR(out_fukuda, "composite", "wet", refF, T);

%Fukuda18_f
% GBS inversion, no activation volume
output_fukuda_f = load('./fukuda18_f.out');
out_fukuda_f = output_fukuda_f(100:2:end,1:end);
m_fukuda_f = out_fukuda_f(:,4); n_fukuda_f = out_fukuda_f(:,5); r_fukuda_f = out_fukuda_f(:,6);
Q_fukuda_f = out_fukuda_f(:,7); A_fukuda_f = out_fukuda_f(:,end);
s_fukuda_MCMC_f = (refE./(10^mean(log10(A_fukuda_f)).*refd^mean(m_fukuda_f).*refF^mean(r_fukuda_f).*exp(-mean(Q_fukuda_f)./(R.*T)))).^(1/mean(n_fukuda_f));
IQR_fukuda_f = calc_IQR(out_fukuda_f, "GBS", refd, refF, T);
%}

figure(3);
%T = T - 273.15; % Switch to C
h=[];

h(1) = semilogx(s_dis_RB,T,'b');
set(gca,'yDir','reverse');
hold on; axis tight;
h(2) = semilogx(s_RB_MCMC,T,'b--');
h(3) = semilogx(s_RB_MCMC_dif,T,'b-.');
X=[10.^IQR_RB(:,1).',fliplr(10.^IQR_RB(:,3).')];  
Y=[T,fliplr(T)];          
fill(X,Y,'b', 'EdgeColor', 'None', 'FaceAlpha', 0.15);
X=[10.^IQR_RB_dif(:,1).',fliplr(10.^IQR_RB_dif(:,3).')];  
Y=[T,fliplr(T)];          
fill(X,Y,'b', 'EdgeColor', 'None', 'FaceAlpha', 0.15);

h(4) = semilogx(s_GT,T, 'g');
h(5) = semilogx(s_GT_MCMC,T,'g--');
X=[10.^IQR_GT(:,1).',fliplr(10.^IQR_GT(:,3).')];  
Y=[T,fliplr(T)];          
fill(X,Y,'g', 'EdgeColor', 'None', 'FaceAlpha', 0.15);  

h(6) = semilogx(s_LP,T,'m');
h(7) = semilogx(s_LP_MCMC,T,'m--');
X=[10.^IQR_LP(:,1).',fliplr(10.^IQR_LP(:,3).')];  
Y=[T,fliplr(T)];          
fill(X,Y,'m', 'EdgeColor', 'None', 'FaceAlpha', 0.15);  

h(8) = semilogx(s_fukuda,T,'r');
h(9) = semilogy(s_fukuda_MCMC,T,'r--');
h(10) = semilogy(s_fukuda_MCMC_f,T,'r-.');
X=[10.^IQR_fukuda(:,1).',fliplr(10.^IQR_fukuda(:,3).')];  
Y=[T,fliplr(T)];          
fill(X,Y,'r', 'EdgeColor', 'None', 'FaceAlpha', 0.15);
X=[10.^IQR_fukuda_f(:,1).',fliplr(10.^IQR_fukuda_f(:,3).')];  
Y=[T,fliplr(T)];          
fill(X,Y,'r', 'EdgeColor', 'None', 'FaceAlpha', 0.15);

h(11) = semilogx(s_renner,T,'k','Linewidth', 1.25);
xlim([1e-4 1e4])
ylim([273 1373])
xlabel('Stress [MPa]')
ylabel('Temperature [K]')
title(['Stress vs. Temperature [K]'])

%text(10^(2.4),400,'LP92','Color','m')
%text(10^(-1.9),600,sprintf('LP92\n(this study)'),'Color','m', 'HorizontalAlignment', 'center')

%text(1e2,530,'GT95','Color','g')
%text(10^(-1),1000,sprintf('GT95\n(this study)'),'Color','g', 'HorizontalAlignment', 'center')

%text(10^(-2.5),1000,'F18','Color','r')
%text(10^(0.3),400,sprintf('F18\n(this study)'),'Color','r', 'HorizontalAlignment', 'center')
%text(10^(-3),400,sprintf('F18\n(this study, GBS)'),'Color','r', 'HorizontalAlignment', 'center')

%text(10^(3.2),580,'RB04','Color','b')
%text(10^(1.8),580,sprintf('RB04\n(this study)'),'Color','b', 'HorizontalAlignment', 'center')
%text(10^(1.5),950,sprintf('RB04\n(this study, diffusion)'),'Color','b', 'HorizontalAlignment', 'center')

%text(10^(-2.8),900,'R02')

legend(h,{'RB04','RB04 (this study)','RB04 (this study, diffusion)','GT95','GT95 (this study)','LP92','LP92 (this study)','F18','F18 (this study)','F18 (this study, GBS)','R02'},'Location','southoutside','NumColumns',2)

%}
