% check_resp.m
clear
run_id = 'rutter04ab';

data = load(['./' run_id '.dat']);
T = data(:,1); dT = data(:,2);
e = data(:,5); de = data(:,6); % strain rate
s = data(:,7); ds = data(:,8); % stress
d = data(:,9); dd = data(:,10); % grain size
id = data(:,11);
%dislocation data (from Rutter04a)
data_dis = data(data(:,11) == 1,:);
T_dis = data_dis(:,1);
dT_dis = data_dis(:,2);
e_dis = data_dis(:,5);
de_dis = data_dis(:,6);
s_dis = data_dis(:,7);
ds_dis = data_dis(:,8);
d_dis = data_dis(:,9);
dd_dis = data_dis(:,10);
%diffusion data (from Rutter04b)
data_dif = data(data(:,11) == 2,:);
T_dif = data_dif(:,1);
dT_dif = data_dif(:,2);
e_dif = data_dif(:,5);
de_dif = data_dif(:,6);
s_dif = data_dif(:,7);
ds_dif = data_dif(:,8);
d_dif = data_dif(:,9);
dd_dif = data_dif(:,10);

output = load(['./rutter04ab_X.out']);
% skip first 1000 runs, then select results from every 100 runs
out = output(100:20:end,1:end);

%1 = id, 2 = simplified chi2, 3 = original chi2, 4 = grain size exponent, 5 = Q diffusion 
%6 = stress exp., 7 = Q dislocation, 8 = A diffusion, 9 = A dislocation
nout = length(out);
chi2 = out(:,3); m = out(:,4); Q_dif = out(:,5); n = out(:,6); Q_dis = out(:,7);
X_dif = out(:,9); X_dis = out(:,8); A_dif = out(:,10); A_dis = out(:,11);
%A_dif = out(:,8); A_dis = out(:,9);

%print summary data
disp(['id=' run_id]);
disp(['m = ' num2str(mean(m)) ' +/- ' num2str(2*std(m))]);
disp(['n = ' num2str(mean(n)) ' +/- ' num2str(2*std(n))]);
disp(['Q_dif = ' num2str(mean(Q_dif)/1e3) ' +/- ' num2str(2*std(Q_dif)/1e3) ' kJ/mol']);
disp(['Q_dis = ' num2str(mean(Q_dis)/1e3) ' +/- ' num2str(2*std(Q_dis)/1e3) ' kJ/mol']);
disp(['log10(A_dif) = ' num2str(mean(log10(A_dif))) ' +/- ' num2str(2*std(log10(A_dif)))]);
disp(['log10(A_dis) = ' num2str(mean(log10(A_dis))) ' +/- ' num2str(2*std(log10(A_dis)))]);
disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(2*std(chi2))]);

R = 8.3145;
nx = 50;
%stress data
xs = logspace(log10(min(s))-2,log10(max(s))+2,nx);
%temperature data
xt = linspace(min(T)-200,max(T)+200,nx);
%grain size data
xd = logspace(log10(min(d))-1,log10(max(d))+2,nx);

mat = [log10(A_dis) log10(A_dif) n m Q_dis Q_dif];
[cc] = corrcoef(mat)


figure(1);
set(gca,'FontSize',20)
%normalize data
normT = 1373;
norms = 250;
normd = 1.0;
normf = 300; % MPa

%strain vs. stress
for i=1:nout
    %model strain rate
    e_pred = A_dif(i)*xs.*normd^-m(i).*exp(-Q_dif(i)/(normT*R)) + A_dis(i)*xs.^n(i).*exp(-Q_dis(i)/(normT*R));
    loglog(xs,e_pred,'r:');
    if i==1
        hold on; box on; axis tight;
    end
end
%mean dislocation and diffusion predictions from inversion
log_avg_dif = mean(log10(A_dif));
log_avg_dis = mean(log10(A_dis));
ye_dis_T = 10^log_avg_dis.*xs.^mean(n).*exp(-mean(Q_dis)/(normT*R));
ye_dif_T_d = 10^log_avg_dif.*xs.*normd^-mean(m).*exp(-mean(Q_dif)/(normT*R));
loglog((xs),(ye_dif_T_d),'k-', 'LineWidth', 1.5);
loglog((xs),(ye_dis_T),'k--', 'LineWidth', 1.5);
%flow law predictions from RB04ab
loglog((xs), (10^(-4.93).*normf.*xs.^(2.97).*exp(-242000/(R*normT))), 'c--', 'LineWidth', 1.5);
loglog((xs), (0.4.*normd^(-2).*xs.*exp(-220000/(R*normT))), 'c-', 'LineWidth', 1.5);

%using Rutter's diffusion data
e_pred_T_d_dif = mean(A_dif)*s_dif.*normd^-mean(m).*exp(-mean(Q_dif)/(normT*R)) + mean(A_dis)*s_dif.^mean(n).*exp(-mean(Q_dis)/(normT*R));
e_pred0_dif = mean(A_dif)*s_dif.*d_dif.^-mean(m).*exp(-mean(Q_dif)./(T_dif*R)) + mean(A_dis)*s_dif.^mean(n).*exp(-mean(Q_dis)./(T_dif*R));

%using Rutter's dislocation data
e_pred_T_d_dis = mean(A_dif)*s_dis.*normd^-mean(m).*exp(-mean(Q_dif)/(normT*R)) + mean(A_dis)*s_dis.^mean(n).*exp(-mean(Q_dis)/(normT*R));
e_pred0_dis = mean(A_dif)*s_dis.*d_dis.^-mean(m).*exp(-mean(Q_dif)./(T_dis*R)) + mean(A_dis)*s_dis.^mean(n).*exp(-mean(Q_dis)./(T_dis*R));

%plot(0.985.*log10(s_dif),log10(e_dif),'g*');
%plot(0.99.*log10(s_dis),log10(e_dis),'go');
%normalized diffusion data
loglog(0.995.*(s_dif),(e_dif.*(e_pred_T_d_dif./e_pred0_dif).*exp(-mean(X_dif))),'bo');
%normalized dislocation data
loglog((s_dis),(e_dis.*(e_pred_T_d_dis./e_pred0_dis)),'bo');

%errorbars
[ex_dis, ey_dis] = calc_err(s_dis, ds_dis, e_dis.*(e_pred_T_d_dis./e_pred0_dis), de_dis.*(e_pred_T_d_dis./e_pred0_dis));
loglog((ex_dis),(ey_dis),'b-');
[ex_dif, ey_dif] = calc_err(s_dif, ds_dif, e_dif.*(e_pred_T_d_dif./e_pred0_dif).*exp(-mean(X_dif)), de_dif.*(e_pred_T_d_dif./e_pred0_dif).*exp(-mean(X_dif)));
loglog(0.995*(ex_dif),(ey_dif),'b-');

title('RB04ab')
ylabel('strain rate [s^{-1}]')
xlabel('stress [MPa]')
text(12, 0.1,['T = ' num2str(normT) ' K' newline 'd = ' num2str(normd) '\mum'], 'FontSize', 20)
xlim([10^(1) 10^(3.5)])
ylim([10^-8 10^(0)])
set(gca, 'LineWidth', 1)
set(gca,'FontSize',20)

e_pred_s_d_dif = mean(A_dif)*norms*normd^-mean(m).*exp(-mean(Q_dif)./(T_dif.*R)) + mean(A_dis)*norms^mean(n).*exp(-mean(Q_dis)./(T_dif.*R));
e_pred_s_d_dis = mean(A_dif)*norms*normd^-mean(m).*exp(-mean(Q_dif)./(T_dis*R)) + mean(A_dis)*norms^mean(n).*exp(-mean(Q_dis)./(T_dis.*R));
figure(2); hold off;

%strain vs. temp
for i=1:nout
    %model strain rate
    e_pred = A_dif(i)*norms*normd^-m(i).*exp(-Q_dif(i)./(xt.*R))+ A_dis(i)*norms.^n(i).*exp(-Q_dis(i)./(xt.*R));
    semilogy(10^4./xt,(e_pred),'r:');
    hold on; box on; axis tight;
end

ye_dis_s = 10^log_avg_dis.*norms^mean(n).*exp(-mean(Q_dis)./(xt.*R));
ye_dif_s_d = 10^log_avg_dif.*norms*normd^-mean(m).*exp(-mean(Q_dif)./(xt.*R));
semilogy(10^4./xt,(ye_dif_s_d),'k-', 'LineWidth', 1.5);
semilogy(10^4./xt,(ye_dis_s),'k--', 'LineWidth', 1.5);
semilogy(10^4./xt, (10^-4.93*normf.*norms.^(2.97).*exp(-242000./(R.*xt))), 'c--', 'LineWidth', 1.5);
semilogy(10^4./xt, (0.4.*normd^(-2).*norms.*exp(-220000./(R.*xt))), 'c-', 'LineWidth', 1.5);

%plot(0.985.*10^4./T_dif,log10(e_dif),'g*');
%plot(0.99.*10^4./T_dis,log10(e_dis),'go');
%normalized diffusion data
semilogy(0.995.*10^4./T_dif,(e_dif.*(e_pred_s_d_dif./e_pred0_dif).*exp(-mean(X_dif))),'bo');
%normalized dislocation data
semilogy(10^4./T_dis,(e_dis.*(e_pred_s_d_dis./e_pred0_dis)),'bo');
%errorbars
[ex_dif, ey_dif] = calc_err(T_dif, dT_dif, e_dif.*(e_pred_s_d_dif./e_pred0_dif).*exp(-mean(X_dif)), de_dif.*(e_pred_s_d_dif./e_pred0_dif).*exp(-mean(X_dif)));
semilogy(0.995.*10^4./ex_dif,(ey_dif),'b-');
[ex_dis, ey_dis] = calc_err(T_dis, dT_dis, e_dis.*(e_pred_s_d_dis./e_pred0_dis), de_dis.*(e_pred_s_d_dis./e_pred0_dis));
semilogy(10^4./ex_dis,(ey_dis),'b-');

ylabel('strain rate [s^{-1}]')
xlabel('10^4/T [K^{-1}]')
text(8.35, 0.09,['\sigma = ' num2str(norms) ' MPa' newline 'd = ' num2str(normd) '\mum'], 'FontSize', 20)
xlim([6 10])
ylim([10^-10 10^(0)])
set(gca, 'LineWidth', 1)
set(gca,'FontSize',20)

e_pred_s_T_dif = mean(A_dif)*norms*d_dif.^-mean(m).*exp(-mean(Q_dif)/(normT*R)) + mean(A_dis)*norms^mean(n)*exp(-mean(Q_dis)/(normT*R));
e_pred_s_T_dis = mean(A_dif)*norms*d_dis.^-mean(m).*exp(-mean(Q_dif)/(normT*R)) + mean(A_dis)*norms^mean(n)*exp(-mean(Q_dis)/(normT*R));

%strain vs. grain size
figure(3); hold off;
for i=1:nout
    %model strain rate
    e_pred = A_dif(i)*norms.*xd.^-m(i).*exp(-Q_dif(i)/(normT*R))+ A_dis(i)*norms^n(i)*exp(-Q_dis(i)/(normT*R));
    loglog((xd),(e_pred),'r:');
    hold on; box on; axis tight;
end

ye_dis_s_T = 10^log_avg_dis.*norms^mean(n)*exp(-mean(Q_dis)/(normT*R));
ye_dif_s_T = 10^log_avg_dif.*norms.*xd.^-mean(m).*exp(-mean(Q_dif)/(normT*R));
loglog((xd),(ye_dif_s_T),'k-', 'LineWidth', 1.5);
loglog((xd), ones(size(xd)).*(ye_dis_s_T),'k--', 'LineWidth', 1.5);
loglog((xd), ones(size(xd)).*(10^-4.93*normf.*norms.^(2.97).*exp(-242000./(R.*normT))), 'c--', 'LineWidth', 1.5);
loglog((xd), ones(size(xd)).*(0.4.*xd.^(-2).*norms.*exp(-220000./(R.*normT))), 'c-', 'LineWidth', 1.5);

%plot(0.99.*log10(d_dif),log10(e_dif),'g*');
%plot(0.985.*log10(d_dis),log10(e_dis),'go');
%normalized diffusion data
loglog(0.995.*(d_dif),(e_dif.*(e_pred_s_T_dif./e_pred0_dif).*exp(-mean(X_dif))),'bo');
%normalized dislocation data
loglog((d_dis),(e_dis.*(e_pred_s_T_dis./e_pred0_dis)),'bo');
%errorbars
[ex_dif, ey_dif] = calc_err(d_dif, dd_dif, e_dif.*(e_pred_s_T_dif./e_pred0_dif).*exp(-mean(X_dif)), de_dif.*(e_pred_s_T_dif./e_pred0_dif).*exp(-mean(X_dif)));
loglog(0.995.*(ex_dif),(ey_dif),'b-');
[ex_dis, ey_dis] = calc_err(d_dis, dd_dis, e_dis.*(e_pred_s_T_dis./e_pred0_dis), de_dis.*(e_pred_s_T_dis./e_pred0_dis));
loglog((ex_dis),(ey_dis),'b-');
ylabel('strain rate [s^{-1}]')
xlabel('grain size [\mum]')
text(6, 0.1,['T = ' num2str(normT) ' K' newline '\sigma = ' num2str(norms) ' MPa'], 'FontSize', 20)
xlim([10^(-1) 10^(2)])
ylim([10^-8 10^(0)])
set(gca, 'LineWidth', 1)
set(gca,'FontSize',20)

figure(4); 
subplot(2,1,1)
%normalized diffusion data
loglog(0.995.*(s_dif),(e_dif.*(e_pred_T_d_dif./e_pred0_dif)),'b*','LineStyle','none'); hold on;
%normalized dislocation data
loglog((s_dis),(e_dis.*(e_pred_T_d_dis./e_pred0_dis)),'r*','LineStyle','none');
title('Without inter-run bias')
ylabel('strain rate [s^{-1}]')
xlabel('stress [MPa]')
ylim([1e-8 1e-0])
set(gca, 'LineWidth', 1)
set(gca,'FontSize',20)

subplot(2,1,2)
%normalized diffusion data
loglog(0.995.*(s_dif),(e_dif.*(e_pred_T_d_dif./e_pred0_dif)).*exp(-mean(X_dif)),'b*','LineStyle','none'); hold on;
%normalized dislocation data
loglog((s_dis),(e_dis.*(e_pred_T_d_dis./e_pred0_dis)),'r*','LineStyle','none');
title('With inter-run bias')
ylabel('strain rate [s^{-1}]')
xlabel('stress [MPa]')
ylim([1e-8 1e-0])
set(gca, 'LineWidth', 1)
set(gca,'FontSize',20)

%{
%histograms
figure(4);
subplot(4,2,1);
hist(m,20); xlabel('m');
subplot(4,2,2);
hist(n,20); xlabel('n');
subplot(4,2,3);
hist(Q_dif,20); xlabel('Q_{dif} [kJ/mol]');
subplot(4,2,4);
hist(Q_dis,20); xlabel('Q_{dis} [kJ/mol]');
subplot(4,2,5);
hist(log10(A_dif),20); xlabel('log_{10}(A_{dif})');
subplot(4,2,6);
hist(log10(A_dis),20); xlabel('log_{10}(A_{dis})');
subplot(4,2,7);
hist(chi2,20); xlabel('\chi^2');
subplot(4,2,8);
hist(X_dif,20); xlabel('X_{dif}');
%}
%}
