clear

% Generate synthetic dataset with 15 datapoints
n_data = 15;
% Fukuda 18 values
A = 10^-2.97; n = 1.7;
m = 0.51; r = 1.0;
Q = 183000;
% Set constants
normT = 1473; % Temperature [K]
normd = 10; % Grain size [um]
norms = 350; % Stress [MPa]
norme = 3e-5;
normP = 300;
normf = 4000; % Water fugacity [MPa]
R = 8.3145;

% Strain stepping data
% Set stress ranges
smin = 20;
smax = 600;
% Set standard dev
sstd = 0.5;
% Vary stress, keep T, d, and f_H2O constant, solve for e
xs = logspace(log10(smin),log10(smax),n_data)+(rand(1,n_data)-0.5).*sstd;
xe = A*xs.^n.*normd^-m.*normf^r.*exp(-Q./(R*normT));
sigma = xs.*(1+0.01*(rand(1, n_data)-0.5)*2);
e_dot = xe.*(1+0.03*(rand(1, n_data)-0.5)*2);

T = normT*ones(n_data, 1);
dT = 7*ones(n_data, 1);
P = normP*ones(n_data, 1);
dP = 0.01*P;
d = normd*ones(n_data, 1);
dd = 0.01*d;
f = normf*ones(n_data, 1);
df = 0.02*f;
id = 1*ones(n_data, 1);

data_1 = [T dT P dP e_dot' 0.03*e_dot' sigma' 0.01*sigma' ...
    d dd f df id];
    
% Temperature stepping data
% Set temperature ranges
tmin = 873;
tmax = 1473;
% Set standard dev
sstd = 0.5;
% Vary T, keep e_dot, d, and f_H2O constant, solve for s
xt = linspace(tmin,tmax,n_data)+(rand(1,n_data)-0.5).*sstd;
xs = (norme./(A.*normd^-m.*normf^r.*exp(-Q./(R.*xt)))).^(1/n);
t = xt.*(1+0.01*(rand(1, n_data)-0.5)*2);
s = xs.*(1+0.01*(rand(1, n_data)-0.5)*2);

e = norme*ones(n_data, 1);
de = 0.01*e;
P = normP*ones(n_data, 1);
dP = 0.01*P;
d = normd*ones(n_data, 1);
dd = 0.01*d;
f = normf*ones(n_data, 1);
df = 0.02*f;
id_2 = 2*ones(n_data, 1);


data_2 = [t' dT P dP e de s' 0.01.*s' ...
    d dd f df id_2];

% grain size stepping data
% set grain size ranges
dmin = 1;
dmax = 100;
% Set standard dev
sstd = 0.5;
% Vary d, keep T, e_dot, f_H2O constant, solve for e
% make grain size a function of T
xd = logspace(log10(dmin),log10(dmax),n_data)+(rand(1,n_data)-0.5).*sstd;
xe = A*norms.^n.*xd.^-m.*normf^r.*exp(-Q./(R*normT));
d = xd.*(1+0.01*(rand(1, n_data)-0.5)*2);
e_dot = xe.*(1+0.03*(rand(1, n_data)-0.5)*2);

T = normT*ones(n_data, 1);
dT = 7*ones(n_data, 1);
s = norms*ones(n_data, 1);
ds = 0.01*s;
P = normP*ones(n_data, 1);
dP = 0.01*P;
f = normf*ones(n_data, 1);
df = 0.02*f;
id_3 = 3*ones(n_data, 1);

data_3 = [T dT P dP e_dot' 0.03*e_dot' s ds ...
    d' dd f df id_3];
%}

% water fugacity stepping data
% set water fugacity ranges
fmin = 1000;
fmax = 5000;
% Set standard dev
sstd = 0.5;
% Vary f_H2O, keep T, d, f_H2O constant, solve for e
xf = logspace(log10(fmin),log10(fmax),n_data)+(rand(1,n_data)-0.5).*sstd;
xe = A*norms.^n.*normd.^-m.*xf.^r.*exp(-Q./(R*normT));
f = xf.*(1+0.01*(rand(1, n_data)-0.5)*2);
e_dot = xe.*(1+0.03*(rand(1, n_data)-0.5)*2);

T = normT*ones(n_data, 1);
dT = 7*ones(n_data, 1);
d = normd*ones(n_data, 1);
dd = 0.01*d;
s = norms*ones(n_data, 1);
ds = 0.01*s;
P = normP*ones(n_data, 1);
dP = 0.01*P;
id_4 = 4*ones(n_data, 1);

data_4 = [T dT P dP e_dot' 0.03*e_dot' s ds ...
    d dd f' df id_4];
%}

data = [data_1; data_2; data_3; data_4];

for i=1:data(end,end)
    out = data(data(:,end) == i,:);
    T_out{i} = out(:,1); % Temperature [K]
    dT_out{i} = out(:,2); 
    e_out{i} = out(:,5); % Strain rate [s^-1]
    de_out{i} = out(:,6);
    s_out{i} = out(:,7); % Stress [MPa]
    ds_out{i} = out(:,8);
    d_out{i} = out(:,9); % Grain size [um]
    dd_out{i} = out(:,10);
    f_out{i} = out(:,11); % Water fugacity [MPa]
    df_out{i} = out(:,12);
end

%%
% Strain vs. stress
figure(1); hold off;
for i=1:out(end,end)
    % normalized data
    norm_factor{i} = (normd^-m.*normf^r.*exp(-Q./(normT.*R)))./...
    (d_out{i}.^-m.*f_out{i}.^r.*exp(-Q./(T_out{i}.*R)));
    loglog(s_out{i},e_out{i}.*norm_factor{i},'o');
    if i == 1
        hold on;
    end
end
title('Strain rate vs. stress')
 
%}
% Strain vs. temp
figure(2); hold off;
for i=1:out(end,end)
    % normalized data
    norm_factor{i} = (norms^n.*normd^-m.*normf^r.*exp(-Q./(T_out{i}.*R)))./...
    (s_out{i}.^n.*d_out{i}.^-m.*f_out{i}.^r.*exp(-Q./(T_out{i}.*R)));
    semilogy(1e4./T_out{i},e_out{i}.*norm_factor{i},'o'); 
    if i == 1
        hold on;
    end
end
title('Strain rate vs. temp')

% Strain vs. grain size
figure(3); hold off;
for i=1:out(end,end)
    % normalized data
    norm_factor{i} = (norms^n*normf^r.*exp(-Q./(normT*R)))./...
    (s_out{i}.^n.*f_out{i}.^r.*exp(-Q./(T_out{i}.*R)));
    loglog(d_out{i},e_out{i}.*norm_factor{i},'o');
    if i == 1
        hold on;
    end
end     
title('Strain rate vs. grain size')

%}

% Strain vs. water fugacity
figure(4); hold off;
for i=1:out(end,end)
    % normalized data
    norm_factor{i} = (norms^n*normd^-m.*exp(-Q./(normT*R)))./...
    (s_out{i}.^n.*d_out{i}.^-m.*exp(-Q./(T_out{i}.*R)));
    loglog(f_out{i},e_out{i}.*norm_factor{i},'o');
    if i == 1
        hold on;
    end
end     
title('Strain rate vs. f_{H2O}')

save('syn_d.out','data','-ascii');
