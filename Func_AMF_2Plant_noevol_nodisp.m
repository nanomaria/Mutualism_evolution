function yout = Func_AMF_2Plant_noevol_nodisp(X)
global q_hp q_cm q_hm q_cp mup mui d rp ALPHA  Aa BETA Ap 
Nm = length(X);
N_plant = length(BETA);
N_AMF = length(ALPHA);

yout = zeros(Nm,1);
P = X(1:N_plant);

alpha = ALPHA';
beta  = BETA';

%% Plant
sP = sum(P.*ones(1,N_plant) - diag(P))';
Gamma_p = Ap;% *P./(sP+(sP<=0))+(P<=0);

%% AMF
for i = 1:N_plant
    M = X(N_plant+1+(i-1)*N_AMF:N_plant+i*N_AMF);
    sM = sum(M.*ones(1,N_AMF) - diag(M))';
    Gamma_m = Aa*M./(sM+(sM<=0));
    gamma_p = Gamma_p(i);
    Mj =  M.*( q_cm*beta(i).*P(i).*Gamma_m./(Gamma_m + sM + ((Gamma_m + sM)<=0)) ...
        - q_hm*alpha.*P(i)./(d+P(i))...
        .*(gamma_p./(gamma_p + sP(i) + ((gamma_p + sP(i))<=0)).*(sP(i)>0) + (sP(i)<=0))  ...
    - mui.*M );
    Pi =  P(i).*(q_hp*rp + q_hp*sum(alpha.*M )./(d+P(i))...
        .*(gamma_p./(gamma_p + sP(i) + ((gamma_p + sP(i))<=0)).*(sP(i)>0) + (sP(i)<=0)) ...
    -q_cp*beta(i).*sum(M.*Gamma_m./(Gamma_m + sM + ((Gamma_m + sM)<=0) )) ...
    -mup*P(i)) ;

    yout(i) = Pi;
    yout(N_plant+1+(i-1)*N_AMF:N_plant+i*N_AMF) = Mj;
end


end