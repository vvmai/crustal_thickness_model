% check_PS94.m

Ps = []; 
rhos1 = []; rhors1 = []; PrhoRTs1=[];
rhos2 = []; rhors2 = []; PrhoRTs2=[];

for P=10.^(-5:0.2:1.6) % GPa
  [rho,rhor,PrhoRT] = calc_rhoH2O(P,800);
  Ps = [Ps P];
  rhos1 = [rhos1 rho];
  rhors1 = [rhors1 rhor];
  PrhoRTs1 = [PrhoRTs1 PrhoRT];

  [rho,rhor,PrhoRT] = calc_rhoH2O(P,1600);
  rhos2 = [rhos2 rho];
  rhors2 = [rhors2 rhor];
  PrhoRTs2 = [PrhoRTs2 PrhoRT];
end

figure(1);
subplot(2,1,1);
plot(rhors1,PrhoRTs1, 'b-', rhors2,PrhoRTs2,'r--');
axis([0 7 0 25]);
xlabel('reduced density');
ylabel('P/\rho RT');
legend('T=800K', 'T=1600K');

subplot(2,1,2);
rhoc = 322; % kg/m^3;
plot(Ps,rhors1*rhoc,'b-', Ps,rhors2*rhoc,'r--');
xlabel('P [GPa]');
ylabel('rho [kg/m^3]');
axis([0 25 0 3000]);
