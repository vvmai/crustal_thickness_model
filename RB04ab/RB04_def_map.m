clear
% Rutter & Brodie's parameters
% Dislocation:
    % Q = 242000;
    % A = 10^-4.93;
    % n = 2.97;
% Diffusion
    % Q = 220000; 
    % A = 0.4;
    % m = 2;
output = load(['./rutter04ab_X.out']);
% skip first 1000 runs, then select results from every 100 runs
out = output(100:10:end,1:end);
%1 = id, 2 = simplified chi2, 3 = original chi2, 4 = grain size exponent, 5 = Q diffusion 
%6 = stress exp., 7 = Q dislocation, 8 = A diffusion, 9 = A dislocation
m = mean(out(:,4)); Q_dif = mean(out(:,5)); n = mean(out(:,6)); Q_dis = mean(out(:,7));
A_dif = mean(log10(out(:,10))); A_dis = mean(log10(out(:,11)));

R = 8.3145;
%reference temperature
refT = 800;

figure(1);
%RB04 deformation map
subplot(1,2,1);
%stress data
xs = logspace(-2,4,50);
%grain size data
xd = logspace(-1,6,50);
[XS,XD] = meshgrid(xs,xd);
%diffusion flow law
E_dif = 0.4*XS.*XD.^-2.*exp(-220000/(R*refT));
%dislocation flow law
E_dis = 10^-4.93*XS.^2.97.*exp(-242000/(R*refT));
%magnitude difference between diffusion and dislocation
crit_line = log10(E_dif)-log10(E_dis);
%delete diffusion data where dislocation dominates
E_dif(crit_line <= 0) = NaN;
%delete dislocation data where diffusion dominates
E_dis(crit_line >= 0) = NaN;
%grain size vs. stress contour plot for strain rates
contour(log10(XD),log10(XS),log10(E_dif),'k--','ShowText','on'); hold on;
contour(log10(XD),log10(XS),log10(E_dis),'k--','ShowText','on');
contour(log10(XD),log10(XS),(crit_line),[0 0],'k-','LineWidth',2);
title(['RB04 Data [T = ' num2str(refT) ' K]'])
ylabel('log_{10}(stress) [MPa]')
xlabel('log_{10}(d) [\mum]')

%Inversion deformation map
subplot(1,2,2);
%diffusion flow law
E_dif = 10^A_dif*XS.*XD.^-m.*exp(-Q_dif/(R*refT));
%dislocation flow law
E_dis = 10^A_dis*XS.^n.*exp(-Q_dis/(R*refT));
%magnitude difference between diffusion and dislocation
crit_line = log10(E_dif)-log10(E_dis);
%delete diffusion data where dislocation dominates
E_dif(crit_line <= 0) = NaN;
%delete dislocation data where diffusion dominates
E_dis(crit_line >= 0) = NaN;
%grain size vs. stress contour plot for strain rates
contour(log10(XD),log10(XS),log10(E_dif),'k--','ShowText','on'); hold on;
contour(log10(XD),log10(XS),log10(E_dis),'k--','ShowText','on');
contour(log10(XD),log10(XS),(crit_line),[0 0],'k-','LineWidth',2);
title(['Inversion Data [T = ' num2str(refT) ' K]'])
ylabel('log_{10}(stress) [MPa]')
xlabel('log_{10}(d) [\mum]')