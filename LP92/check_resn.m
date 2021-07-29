% check_resp.m
clear
run_id = 'LP92';

file = load(['./' run_id '.dat']);
%T = data(:,1); dT = data(:,2);
%e = data(:,5); de = data(:,6); % strain rate (1/s)
%s = data(:,7); ds = data(:,8); % stress (MPa)
%d = data(:,9); dd = data(:,10); % grain size (um)
%COH = data(:,11); dCOH = data(:,12); % water content (H/Si)
%id = data(:,13);

% Separate samples into individual arrays
for i=1:14
    data = file(file(:,13) == i,:);
    T{i} = data(:,1);
    dT{i} = data(:,2);
    e{i} = data(:,5);
    de{i} = data(:,6);
    s{i} = data(:,7);
    ds{i} = data(:,8);
    d{i} = data(:,9);
    dd{i} = data(:,10);
    COH{i} = data(:,11);
    dCOH{i} = data(:,12);
end

output = load(['./' run_id '.out']);
% skip first 1000 runs, then select results from every 100 runs
out = output(100:10:end,1:end);
%1 = id, 2 = simplified chi2, 3 = original chi2
%4 = water exponent, 5 = stress exponent, 6 = Q (J/mol), 7 = A
nout = length(out);
chi2 = out(:,3); r = out(:,4); n = out(:,5); Q = out(:,6); A = out(:,7);

%print summary data
disp(['id=' run_id]);
disp(['n = ' num2str(mean(n)) ' +/- ' num2str(2*std(n))]);
disp(['r = ' num2str(mean(r)) ' +/- ' num2str(2*std(r))]);
disp(['Q = ' num2str(mean(Q)/1e3) ' +/- ' num2str(2*std(Q)/1e3) ' kJ/mol']);
disp(['log10(A) = ' num2str(mean(log10(A))) ' +/- ' num2str(2*std(log10(A)))]);
disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(2*std(chi2))]);

%normalize data
normT = 1373;
norms = 220;
normd = 1.3;
normCOH = 4000e-6;
R = 8.3145;

figure(1)

for i=1:14
    norm_factor{i} = (normCOH^mean(r)*exp(-mean(Q)/(normT*R)))./(COH{i}.^mean(r).*exp(-mean(Q)./(T{i}.*R)));
    loglog(s{i},e{i}.*norm_factor{i},'-o');
    %loglog(s{i}*0.99,e{i},'ro');
    if i == 1
        hold on;
    end
end

ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')
xlim([95,800])

%histograms
figure(4);
subplot(3,2,1);
hist(n,20); xlabel('n');
subplot(3,2,2);
hist(r,20); xlabel('r');
subplot(3,2,3);
hist(Q/1e3,20); xlabel('Q [kJ/mol]');
subplot(3,2,4);
hist(log10(A),20); xlabel('log_{10}(A)');
subplot(3,2,5);
hist(chi2,20); xlabel('\chi^2');