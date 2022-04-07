function yout = Func_AMF_Plant_evol_alpha_beta_comp_nodisp_comm(X)
global q_hp q_cm q_hm q_cp mup mui d rp ALPHA1 ALPHA2 Ad_alpha1 Ad_alpha2 dm1 dm2 Aa1 Aa2 BETA1 BETA2 Ap1 Ap2 Ad_beta1 Ad_beta2 dp1 dp2 C1 C2

Nm = length(X);
N_plant1 = length(BETA1);
N_AMF1 = length(ALPHA1);
N_plant2 = length(BETA2);
N_AMF2 = length(ALPHA2);

yout = zeros(Nm,1);



P1= X(1:N_plant1);
M1 = X(N_plant1+1:(N_plant1+N_AMF1));
P2 = X((N_plant1+N_AMF1+1):(N_plant1+N_AMF1+N_plant2));
M2 = X((N_plant1+N_AMF1+N_plant2+1):(N_plant1+N_AMF1+N_plant2+N_AMF2));


alpha1 = ALPHA1';
beta1 = BETA1';
alpha2 = ALPHA2';
beta2  = BETA2';

sM1 = sum(M1.*ones(1,N_AMF1) - diag(M1))';
Gamma_m1 = Aa1*M1./(sM1+(sM1<=0));
sP1 = sum(P1.*ones(1,N_plant1) - diag(P1))';
Gamma_p1 = Ap1*P1./(sP1+(sP1<=0));
sM2 = sum(M2.*ones(1,N_AMF2) - diag(M2))';
Gamma_m2 = Aa2*M2./(sM2+(sM2<=0));
sP2 = sum(P2.*ones(1,N_plant2) - diag(P2))';
Gamma_p2 = Ap2*P2./(sP2+(sP2<=0));

%.*Gamma_m2./(Gamma_m2 + sM2 + ((Gamma_m2 + sM2)<=0))

P1i = dp1*Ad_beta1*P1 + P1.*(q_hp*rp + q_hp*sum(alpha1.*M1 )./(d+P1).*C2./(C2+sM2) ...
    -q_cp*beta1.*M1...
    -mup*P1) ;
M1j = dm1*Ad_alpha1*M1 + M1.*( q_cm*sum(beta1.*P1) ...
    - q_hm*alpha1.*sum(P1./(d+P1)) - mui.*M1);

P2i = dp2*Ad_beta2*P2 + P2.*(q_hp*rp + q_hp*sum(alpha2.*M2 )./(d+P2).*C1./(C1+sM1)...
    -q_cp*beta2.*M2 ...
    -mup*P2) ;
M2j = dm2*Ad_alpha2*M2 + M2.*( q_cm*sum(beta2.*P2) ...
    - q_hm*alpha2.*sum(P2./(d+P2))  - mui.*M2);


yout(1:N_plant1) = P1i;
yout(N_plant1+1:N_plant1+N_AMF1) = M1j;

yout(N_plant1+N_AMF1+1:N_plant1+N_AMF1+N_plant2) = P2i;
yout(N_plant1+N_AMF1+N_plant2+1:N_plant1+N_AMF1+N_plant2+N_AMF2) = M2j;

end