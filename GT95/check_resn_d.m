% check_resn.m
clear
run_id = 'GT95_d';

file = load(['./' run_id '.dat']);
T_all = file(:,1);
P_all = file(:,3)*1e9; % Pressure (Pa)
e_all = file(:,5); % Strain rate (s^-1)
s_all = file(:,7); % Stress (MPa)
d_all = file(:,9); % Stress (MPa)

% Separate samples into arrays nested in a cell array
for i=1:5
    data = file(file(:,11) == i,:);
    T{i} = data(:,1);
    dT{i} = data(:,2);
    P{i} = data(:,3)*1e9; % Pressure [Pa]
    dP{i} = data(:,4)*1e9;
    e{i} = data(:,5);
    de{i} = data(:,6);
    s{i} = data(:,7);
    ds{i} = data(:,8);
    d{i} = data(:,9); % Grain size [um]
    dd{i} = data(:,10);
end

output = load(['./' run_id '_dif.out']);
% skip first 1000 runs, then select results from every 100 runs
%out = output(100:10:end,1:end);
out=output;
nout = length(out);
% 4 = Stress exponent, 5 = Activation energy (J/mol), 21 = A
chi2 = out(:,3); m = out(:,4); Q = out(:,5); A = out(:,11);
% 6-10 = Inter-run bias relative to sample 1
for i=6:10
    X{i-5} = out(:,i);
end

% Print summary data
disp(['id=' run_id]);
disp(['m = ' num2str(mean(m)) ' +/- ' num2str(1*std(m))]);
disp(['Q = ' num2str(mean(Q)/1e3) ' +/- ' num2str(1*std(Q)/1e3) ' kJ/mol']);
disp(['log10(A) = ' num2str(mean(log10(A))) ' +/- ' num2str(1*std(log10(A)))]);
disp(['chi2 = ' num2str(mean(chi2)) ' +/- ' num2str(1*std(chi2))]);
for i=2:5
    disp(['X_' num2str(i) ' = ' num2str(mean(X{i})) ' +/- ' num2str(1*std(X{i}))]);
end

R = 8.3145;
nx = 20;
%stress data
xs = linspace(min(s_all)-15,max(s_all)+150,nx);
%temperature data
xt = linspace(min(T_all)-100,max(T_all)+300,nx);
%grain size data
xd = linspace(min(d_all)-1,max(d_all)+30,nx);
%normalize data
normT = 1373;
norms = 220;
normP = 1.5*1e9; % Pressure [Pa]
normd = 1.3; % Grain size [um]

figure(1);
% Strain vs. stress

for i=1:nout
    %model strain rate
    e_pred = A(i)*xs.*normd^-m(i).*exp(-Q(i)./(normT*R));
    loglog(xs,e_pred,'r:');
    if i==1
        hold on; box on; axis tight;
    end
end


for i=1:5
    e_pred_T_d_P{i} = A(i)*s{i}.*normd^-m(i).*exp(-Q(i)./(normT*R));
    e_pred0{i} = A(i)*s{i}.*d{i}.^-m(i).*exp(-Q(i)./(T{i}.*R));
    norm_factor{i} = e_pred_T_d_P{i}./e_pred0{i};
    loglog(s{i},e{i}.*norm_factor{i}.*exp(-mean(X{i})),'bo');
    if i == 1
        hold on; box on; axis tight;
    end
end
ylabel('Strain rate [s^{-1}]')
xlabel('Stress [MPa]')

figure(2);
subplot(2,2,1);
for i=1:nout
    param_1 = log10(A(i));
    plot(chi2(i),param_1,'r*');
    if i==1
        hold on;
    end
end
xlabel('chi2')
ylabel('A')

subplot(2,2,2);
for i=1:nout
    param_2 = Q(i)/1e3;
    plot(chi2(i),param_2,'r*');
    if i==1
        hold on;
    end
end
xlabel('chi2')
ylabel('Q')

subplot(2,2,3); hold off;
for i=1:nout
    param_3 = m(i);
    plot(chi2(i),param_3,'r*');
    if i==1
        hold on;
    end
end
xlabel('chi2')
ylabel('m')

subplot(2,2,4); hold off;
for i=1:nout
    param_4 = A(i)*norms*normd^-m(i).*exp(-Q(i)/(normT*R));
    plot(chi2(i),param_4,'r*');
    hold on;
end
xlabel('chi2')
ylabel('Strain')

 % Histograms
figure(4);
subplot(3,2,1);
hist(m,20); xlabel('m');
subplot(3,2,2);
hist(Q/1e3,20); xlabel('Q [kJ/mol]');
subplot(3,2,3);
hist(log10(A),20); xlabel('log_{10}(A)');
subplot(3,2,4);
hist(chi2,20); xlabel('\chi^2');
%}