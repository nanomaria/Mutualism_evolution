function yout = Func_AMF_Plant_evol_alpha_beta_comp_nodisp(X)
global q_hp q_cm q_hm q_cp mup mui d rp ALPHA Ad_alpha dm Aa BETA Ap Ad_beta dp
Nm = length(X);
N_plant = length(BETA);
N_AMF = Nm-N_plant;

yout = zeros(Nm,1);
P = X(1:N_plant);
M = X(N_plant+1:end);
alpha = ALPHA';
beta  = BETA';
sM = sum(M.*ones(1,N_AMF) - diag(M))';
Gamma_m = Aa*M./(sM+(sM<=0));
sP = sum(P.*ones(1,N_plant) - diag(P))';
Gamma_p = Ap*P./(sP+(sP<=0));
Pi = dp*Ad_beta*P + P.*(q_hp*rp + q_hp*sum(alpha.*M )./(d+P).*Gamma_p./(Gamma_p + sP + ((Gamma_p + sP)<=0)) ...
    -q_cp*beta.*sum(M.*Gamma_m./(Gamma_m + sM + ((Gamma_m + sM)<=0) )) ...
    -mup*P) ;
Mj = dm*Ad_alpha*M + M.*( q_cm*sum(beta.*P).*Gamma_m./(Gamma_m + sM + ((Gamma_m + sM)<=0) )...
    - q_hm*alpha.*sum(P./(d+P).*Gamma_p./(Gamma_p + sP + ((Gamma_p + sP)<=0)))  - mui.*M );

yout(1:N_plant) = Pi;
yout(N_plant+1:end) = Mj;

end