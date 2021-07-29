clear
run_id = 'fukuda18';
% load data file
    file = load(['./' run_id '.dat']);
    T = file(:,1); % Temperature (K)
    P = file(:,3); % Pressure (GPa)
    f = file(:,11); % water fugacity [MPa]

for i = 1:length(file)
    calc_f(i) = calc_fH2O(P(i),T(i));
end

disp("Caculated/fukuda18 f_H2O: ")
fprintf('%2.4f\n',transpose(calc_f)./f)