% check_res_n.m

run_id = 'rutter04a';

data = load(['./' run_id '.dat']);
T = data(:,1); dT = data(:,2);
p = data(:,3); dp = data(:,4);
e = data(:,5); de = data(:,6); % strain rate
s = data(:,7); ds = data(:,8); % stress
d = data(:,9); dd = data(:,10); % grain size
id = data(:,11);

out = load(['./' run_id '.out']);
%1 = id, 2 = simplified chi2, 3 = original chi2, 4 = n, 5 = A
nout = length(out(:,1));
chi2 = out(:,3); n = out(:,4); A = out(:,6);
%6 = activation energy (J)
Q = out(:,5);

%print summary data
disp(['id=' run_id]);
disp(['n = ' num2str(mean(n)) ' +/- ' num2str(2*std(n))]);
disp(['Q = ' num2str(mean(Q)/1e3) ' +/- ' num2str(2*std(Q)/1e3) ' kJ/mol']);
disp(['log10(A) = ' num2str(mean(log10(A))) ' +/- ' num2str(2*std(log10(A)))]);
disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(2*std(chi2))]);

figure(1); hold off; 
%flow law based on model predictions
nx = 20;
%stress data
xs = linspace(min(s)-10,max(s)+50,nx);
%temperature data
t = linspace(min(T)-10,max(T)+25,nx);
R = 8.3145;
%normalize data
normT = 1473;
normT_factor = (exp(-mean(Q)/(normT*R))./exp(-mean(Q)./(T.*R)));
for i=1:length(chi2)
    %model strain rate
    xe = A(i)*xs.^n(i).*exp(-Q(i)./(normT*R));
    loglog(xs,xe,'r:');
    if (i == 1)
        hold on; axis tight;
    end
end

%loglog(s*0.98,e,'bo');
loglog(s, e.*normT_factor,'bo');
loglog(s*0.98, e, 'ro');
ylabel('log(strain rate)')
xlabel('log(stress)')
loglog(xs, 10^-4.93*300.*xs.^(2.97).*exp(-242000/(R*normT)), 'g');

%error
[ex, ey_normT] = calc_err(s,ds,e.*normT_factor,de.*normT_factor);
[ex, ey] = calc_err(s,ds,e,de);
loglog(ex,ey_normT,'b-');
loglog(ex*0.98,ey,'r-');
%[ex_normT,ey_normT] = calc_err(s,ds,e_normT,(de./e).*e_normT);
%loglog(ex_normT,ey_normT,'g-');
% Rutter04a flow law

%temp vs. strain
figure(2); hold off;

norms = 180;
norms_factor = (norms^mean(n))./(s.^mean(n));
e_norms = e.*norms_factor;
for i=1:length(chi2)
    %model strain rate
    xe = A(i)*norms^n(i).*exp(-Q(i)*(t.*R).^-1);
    semilogy(10^4./t, xe,'r:');
    if (i == 1)
        hold on; axis tight;
    end
end
semilogy(10^4./T.*0.998, e,'bo');
semilogy(10^4./T, e_norms,'ro');
xlabel('10^4/T');
ylabel('log_{10}(strain rate)');

%error
[ex,ey] = calc_err(T,dT,e,de);
plot(10^4./ex*0.998,log10(ey),'b-');
[ex_norms,ey_norms] = calc_err(T,dT,e_norms,(de./e).*e_norms);
plot(10^4./ex_norms,log10(ey_norms),'g-');

%histograms
figure(3);
subplot(2,2,1);
hist(n,20); xlabel('n');
subplot(2,2,2);
hist(log10(A),20); xlabel('log_{10}(A)');
subplot(2,2,3);
hist(chi2,20); xlabel('\chi^2');
subplot(2,2,4);
hist(Q/1e3,20); xlabel('Q [kJmol^{-1}]');