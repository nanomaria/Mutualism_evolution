function yout = Func_AMF_Plant_evol_alpha_comp_continuous(X)
global q_hp q_cm q_hm q_cp mup mui d rp ALPHA Ad_alpha dm Aa BETA dalpha
Nm = length(X);
N_plant = length(BETA);
N_AMF = Nm-N_plant;

yout = zeros(Nm,1);
P = X(1:N_plant);
M = X(N_plant+1:end);
alpha = ALPHA';
beta  = BETA';
sM = sum(M*dalpha)';
Gamma_m = Aa*M*dalpha./(sM+(sM<=0));
Pi = P.*(q_hp*rp + q_hp*sum(alpha.*M*dalpha )./(d+P) ...
    -q_cp*beta.*sum(M.*Gamma_m*dalpha./(Gamma_m + sM + ((Gamma_m + sM)<=0) )) ...
    -mup*P) ;
Mj = dm*Ad_alpha*M + M.*( (q_cm*beta.*Gamma_m./(Gamma_m + sM + ((Gamma_m + sM)<=0) )...
    - q_hm*alpha./(d+P)).*P  - mui.*M ); %

yout(1:N_plant) = Pi;
yout(N_plant+1:end) = Mj;

end