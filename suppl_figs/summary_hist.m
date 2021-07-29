% Inversion flow laws
%RB04ab
% composite inversion, use only dislocation part
output_RB = load('./rutter04ab_X.out');
out_RB = output_RB(100:10:end,1:end);
m_RB = out_RB(:,4); Q_dif_RB = out_RB(:,5); n_RB = out_RB(:,6); Q_dis_RB = out_RB(:,7);
A_dif_RB = log10(out_RB(:,10)); A_dis_RB = log10(out_RB(:,11)); chi2_RB = out_RB(:,3);

%GT95
% composite inversion, using dislocation part, no avtivation volume
% (dislocation-only inversion did not use corrected stress)
output_GT = load('./GT95_d_fixed.out');
out_GT = output_GT(100:10:end,1:end);
m_GT = out_GT(:,4); Q_dif_GT = out_GT(:,5); 
%V_dif_GT = out_GT(:,6);
n_GT = out_GT(:,6); Q_dis_GT = out_GT(:,7);
%V_dis_GT = out_GT(:,9);
A_dif_GT = out_GT(:,end-1); A_dis_GT = out_GT(:,end);
chi2_GT = out_GT(:,3);

%LP92
% composite inversion, no activation volume, use only dislocation part
%output_LP = load('./LP92_pXC.out');
output_LP = load('./LP92_COH.out');
out_LP = output_LP(100:10:end,1:end);
% 4 = Water exponent, 5 = Stress exponent, 6 = Activation energy (J/mol), 21 = A
chi2_LP = out_LP(:,3); r_LP = out_LP(:,4); n_LP = out_LP(:,5); Q_LP = out_LP(:,6); A_LP = out_LP(:,21);

%Fukuda18
% composite inversion, no activation volume, use only dislocation part
% COH: assumes samples are saturated
output_fukuda = load('./fukuda18_d_fixed.out');
out_fukuda = output_fukuda(100:2:end,1:end);
r_dif_fukuda = out_fukuda(:,4); m_fukuda = out_fukuda(:,5); Q_dif_fukuda = out_fukuda(:,6);
r_dis_fukuda = out_fukuda(:,7); n_fukuda = out_fukuda(:,8); Q_dis_fukuda = out_fukuda(:,9);
A_dis_fukuda = out_fukuda(:,end); A_dif_fukuda = out_fukuda(:,end-1);
chi2_fukuda = out_fukuda(:,3);

%Fukuda18_f
% GBS inversion, no activation volume
output_fukuda_f = load('./fukuda18_f_d_fixed.out');
out_fukuda_f = output_fukuda_f(100:10:end,1:end);
m_fukuda_f = out_fukuda_f(:,4); n_fukuda_f = out_fukuda_f(:,5); r_fukuda_f = out_fukuda_f(:,6);
Q_fukuda_f = out_fukuda_f(:,7); A_fukuda_f = out_fukuda_f(:,end);
chi2_fukuda_f = out_fukuda_f(:,3);
%}

% Histograms
figure(1);
% n
    subplot(9,5,1); 
    hist(n_fukuda,20); title('F18'); ylabel('n', 'fontweight', 'bold');
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,2);  
    hist(n_GT,20); title('GT95');
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,3); 
    hist(n_RB,20); title('RB04ab');
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,4); 
    hist(n_LP,20); title('LP92');
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,5);
    hist(n_fukuda_f,20); title('F18 (GBS)');
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
% m 
    subplot(9,5,6); 
    hist(m_fukuda,20); 
    ylabel('m', 'fontweight', 'bold');
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,7);
    hist(m_GT,20); 
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,8);
    hist(m_RB,20); 
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,10);
    hist(m_fukuda_f,20); 
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
% r_dis
    subplot(9,5,11); 
    hist(r_dis_fukuda,20); 
    ylabel('r_{dis}', 'fontweight', 'bold'); 
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,14);
    hist(r_LP,20);
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,15);
    hist(r_fukuda_f,20);
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
% r_dif
    subplot(9,5,16);
    hist(r_dif_fukuda,20);
    ylabel('r_{dif}', 'fontweight', 'bold');
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
% Q_dis
    subplot(9,5,21); 
    hist(Q_dis_fukuda/1e3,20); 
    ylabel('Q_{dis}', 'fontweight', 'bold');
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,22);
    hist(Q_dis_GT/1e3,20);
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,23);
    hist(Q_dis_RB/1e3,20); 
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,24);
    hist(Q_LP/1e3,20);
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,25);
    hist(Q_fukuda_f/1e3,20); 
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
% Q_dif
    subplot(9,5,26);
    hist(Q_dif_fukuda/1e3,20);
    ylabel('Q_{dif}', 'fontweight', 'bold');
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,27);
    hist(Q_dif_GT/1e3,20); 
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,28);
    hist(Q_dif_RB/1e3,20); 
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
% A_dis
    subplot(9,5,31);
    hist(log10(A_dis_fukuda),20); 
    ylabel('logA_{dis}', 'fontweight', 'bold');
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,32);
    hist(log10(A_dis_GT),20); 
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,33);
    hist(A_dis_RB,20); 
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,34);
    hist(log10(A_LP),20);
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
    subplot(9,5,35);
    hist(log10(A_fukuda_f),20); 
    set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
    set(gca,'FontSize',10)
% A_dif
   subplot(9,5,36);
   fukuda_A = A_dif_fukuda > 10^-10;
   %hist(log10(A_dif_fukuda(fukuda_A)),20);
    hist(log10(A_dif_fukuda),20);
    ylabel('logA_{dif}', 'fontweight', 'bold');
   set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
   set(gca,'FontSize',10)
   subplot(9,5,37);
   GT_A = A_dif_GT > 10^-20;
    %hist(log10(A_dif_GT(GT_A)),20);
   hist(log10(A_dif_GT),20);
   set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1) 
   set(gca,'FontSize',10)
   subplot(9,5,38);
   hist(A_dif_RB,20); 
   set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
   set(gca,'FontSize',10)
% chi_2
   subplot(9,5,41);
   hist(chi2_fukuda,20);
   ylabel('\chi^2/N', 'fontweight', 'bold');
   set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
   set(gca,'FontSize',10)
   subplot(9,5,42);
   hist(chi2_GT,20);
   set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
   set(gca,'FontSize',10)
   subplot(9,5,43);
   hist(chi2_RB,20);
   set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
   set(gca,'FontSize',10)
   subplot(9,5,44);
   hist(chi2_LP,20);
   set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
   set(gca,'FontSize',10)
   subplot(9,5,45);
   hist(chi2_fukuda_f,20);
   set(gca,'TickLength',[0.1 0.1])
    set(gca, 'LineWidth', 1)
   set(gca,'FontSize',10)