%
function [diff] = calc_PS94_eq2RHS(rho,c,LHS)

rho2=rho^2;

diff = rho + rho2*(c(1) ...
	 - ((c(3)+2*c(4)*rho+(3*c(5)+4*c(6)*rho)*rho2) ...
	 / (c(2)+c(3)*rho+(c(4)+c(5)*rho+c(6)*rho2)*rho2)^2) ...
	 + c(7)*exp(-c(8)*rho) ...
	 + c(9)*exp(-c(10)*rho)) - LHS;
diff = abs(diff);