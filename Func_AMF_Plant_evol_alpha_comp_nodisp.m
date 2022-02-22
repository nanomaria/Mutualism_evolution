function yout = Func_AMF_Plant_evol_alpha_comp_nodisp(X)
global q_hp q_cm q_hm q_cp beta mup mui d rp ALPHA Ad_alpha dm Aa 
Nm = length(X);
N_AMF = Nm-1;
yout = zeros(Nm,1);
P = X(1);
M = X(2:end);
alpha = ALPHA';
sM = sum(M.*ones(1,N_AMF) - diag(M))';
Gamma = Aa*M./(sM+(sM<=0));
pnew =  P.*(q_hp*rp + q_hp*sum(alpha.*M)./(d+P) ...
    -q_cp*beta*sum(M.*Gamma./(Gamma + sM + ((Gamma + sM)<=0) )) ...
    -mup*P) ;
Mj = dm*Ad_alpha*M + M.*( (q_cm*beta*Gamma./(Gamma + sM + ((Gamma + sM)<=0) )...
    - q_hm*alpha./(d+P)).*P - mui.*M );

yout(1) = pnew;
yout(2:end) = Mj;

end