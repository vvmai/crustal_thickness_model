% check_res_n.m

run_id = 'rutter04b';

data = load(['./' run_id '.dat']);
T = data(:,1); dT = data(:,2);
p = data(:,3); dp = data(:,4);
e = data(:,5); de = data(:,6); % strain rate
s = data(:,7); ds = data(:,8); % stress
d = data(:,9); dd = data(:,10); % grain size
id = data(:,11);

output = load(['./' run_id '.out']);
out = output(100:20:end,1:end);
%1 = id, 2 = simplified chi2, 3 = original chi2, 4 = grain size exponent, 5 = activation energy (J/mol), 6 = A
nout = length(out);
chi2 = out(:,3); m = out(:,4); A = out(:,6); Q = out(:,5);

%print summary data
disp(['id=' run_id]);
disp(['m = ' num2str(mean(m)) ' +/- ' num2str(2*std(m))]);
disp(['Q = ' num2str(mean(Q)/1e3) ' +/- ' num2str(2*std(Q)/1e3) ' kJ/mol']);
disp(['log10(A) = ' num2str(mean(log10(A))) ' +/- ' num2str(2*std(log10(A)))]);
disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(2*std(chi2))]);

R = 8.3145;
nx = 40;
%stress data
xs = linspace(min(s)-20,max(s)+200,nx);
%temperature data
xt = linspace(min(T)-100,max(T)+150,nx);
%grain size data
xd = linspace(min(d)-0.1,max(d)+0.5,nx);
normT = 1373;
normT_factor =(exp(-mean(Q)/(normT*R))./exp(-mean(Q)./(T.*R)));
norms = 220;
norms_factor = norms./s; %assuming n=1
normd = 1.3;
normd_factor = (d./normd).^mean(m);

figure(1);
for i=1:nout
    %model strain rate
    ye = A(i)*xs*normd^-m(i).*exp(-Q(i)./(normT.*R));
    plot(log10(xs),log10(ye),'r:');
    hold on; axis tight; box on;
end
%strain vs. stress
plot(log10(s*0.97),log10(e),'gd');
plot(log10(s),log10(e.*normT_factor.*normd_factor),'bo');
plot(log10(xs), log10(0.4.*normd^(-2).*xs.*exp(-220000/(R*normT))), 'c-', 'LineWidth', 1.5);

ylabel('log(strain rate)')
xlabel('log(stress)')
ylim([-9,-2]);
xlim([1.3,2.9]);


hold off; figure(2);
for i=1:nout
    %model strain rate
    ye = A(i)*norms*normd^-m(i).*exp(-Q(i)./(xt.*R));
    plot(1e4.*xt.^-1,log10(ye),'r:');
    hold on; axis tight; box on;
end
%strain vs. temp
plot(1e4./T*0.997,log10(e),'gd');
plot(1e4./T, log10(e.*norms_factor.*normd_factor),'bo');
plot(1e4./T, log10(0.4.*normd^(-2).*norms.*exp(-220000./(R.*T))), 'c-', 'LineWidth', 1.5);
ylabel('log(strain rate)')
xlabel('10^4/T [K^{-1}]')
ylim([-9,-2.5])
xlim([6.5,9.5])

hold off; figure(3);
for i=1:nout
    %model strain rate
    ye = A(i)*norms.*xd.^-m(i).*exp(-Q(i)./(normT*R));
    plot(log10(xd),log10(ye),'r:');
    hold on; axis tight; box on;
end

%strain vs. d
plot(log10(d*0.97),log10(e),'gd');
plot(log10(d),log10(e.*norms_factor.*normT_factor),'bo'); 
plot(log10(d), log10(0.4.*d.^(-2).*norms.*exp(-220000./(R.*normT))), 'c-', 'LineWidth', 1.5);
ylabel('log(strain rate)')
xlabel('log(d) [\mum]')
ylim([-10,-1]);
xlim([-.4,.7]);
%}
%{
histograms
figure(4);
subplot(2,2,1);
hist(m,20); xlabel('m');
subplot(2,2,2);
hist(log10(A),20); xlabel('log_{10}(A)');
subplot(2,2,3);
hist(chi2,20); xlabel('\chi^2');
subplot(2,2,4);
hist(Q/1e3,20); xlabel('Q [kJmol^{-1}]');
%}