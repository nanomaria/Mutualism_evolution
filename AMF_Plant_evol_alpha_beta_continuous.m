clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Tf = 1000;
%% Parameter of the model
global q_hp q_cm q_hm q_cp mup mui d rp Aa a Ad_alpha ALPHA dm BETA dalpha
q_hp  = 3; % q>1
q_cm  = 2;
q_hm  = 1;
q_cp  = 1;
%% Trait beta
beta_min = 0;
beta_max = 10;
dbeta = 0.1;
BETA = beta_min:dbeta:beta_max;
Nbeta = length(BETA);
% alpha = beta;
%% Trait alpha
alphamin = 0;
alphamax = 10;
dalpha = 0.1;
ALPHA  = alphamin:dalpha:alphamax;
Nalpha = length(ALPHA);

mup = 1/100; % 1/100 %0.3
mui = 1/20; % 1/20 % 0.03

d = 1.2;
rp = 0;
%% Competition terms
a = 0.2;
ap = 0.2;

N_AMF = Nalpha;
N_plant = Nbeta;
Aa = a*(ones(N_AMF,N_AMF)-diag(ones(1,N_AMF)));
Ap = ap*(ones(N_plant,N_plant)-diag(ones(1,N_plant)));


% Diffusion matrix alpha 
e = ones(Nalpha,1);
I_alpha  = spdiags(e,0,Nalpha,Nalpha);
Ad_alpha = spdiags([e -2*e e],-1:1,Nalpha,Nalpha);
Ad_alpha(1,1) = -1;
Ad_alpha(end,end) = -1;
Ad_alpha = Ad_alpha/(dalpha^2);

% Diffusion matrix beta 
e = ones(Nbeta,1);
I_beta  = spdiags(e,0,Nbeta,Nbeta);
Ad_beta = spdiags([e -2*e e],-1:1,Nbeta,Nbeta);
Ad_beta(1,1) = -1;
Ad_beta(end,end) = -1;
Ad_beta = Ad_beta/(dbeta^2);

dm = @(m) 1e-3;  % mutation rate AMF
dm_min = 1e-3;
dm_max = 1;
% dm = @(m) (1-exp(-m))*(dm_max-dm_min)+dm_min;
dp = 1e-3;  % mutation rate plant

%% Functionnal response and interactions MARIA
choice = 3;
sM = @(M) sum((ones(N_AMF,N_AMF)-diag(ones(1,N_AMF)))*M)*dalpha;
sP = @(P) sum((ones(N_plant,N_plant)-diag(ones(1,N_plant)))*P)*dbeta;
m = @(alpha) abs(alpha);
if (choice == 1)
    %% Competition between AMF and Plant
    Gamma_m = @(M) Aa*M*dalpha./(sM(M)+(sM(M)<=0));
    Gamma_p = @(P) Ap*P*dbeta./(sP(P)+(sP(P)<=0));
    fp = @(alpha,beta,P,M) P.*(q_hp*rp + q_hp*sum(m(alpha).*M)*dalpha./(d+P).*Gamma_p(P)./(Gamma_p(P) + sP(P) + ((Gamma_p(P) + sP(P))<=0)) ...
        -q_cp*beta.*sum(M.*Gamma_m(M)./(Gamma_m(M) + sM(M) + ((Gamma_m(M) + sM(M))<=0) ))*dalpha ...
        -mup*sP(P)) ;
    fm = @(alpha,beta,P,M) M.*( (q_cm*sum(beta.*P*dbeta).*Gamma_m(M)./(Gamma_m(M) + sM(M) + ((Gamma_m(M) + sM(M))<=0) )...
        - q_hm*m(alpha).*sum(P*dbeta./(d+P).*Gamma_p(P)./(Gamma_p(P) + sP(P) + ((Gamma_p(P) + sP(P))<=0)))) - mui.*sM(M) );
elseif (choice == 2)
    %% No competition
    fp = @(alpha,beta,P,M) P.*(q_hp*rp + q_hp*sum(alpha.*M*dalpha)./(d+P) - q_cp*beta*sum(M)*dalpha ...
        -mup*sP(P)) ;
    fm = @(alpha,beta,P,M) M.*( (q_cm*sum(beta.*P*dbeta) - q_hm*m(alpha).*sum(P./(d+P)*dbeta)) - mui.*sM(M) );
elseif(choice == 3)
        %% Competition between AMF and Plant
    Gamma_m = @(M) Aa*M*dalpha./(sM(M)+(sM(M)<=0));
    Gamma_p = @(P) Ap*P*dbeta./(sP(P)+(sP(P)<=0));
    fp = @(alpha,beta,P,M) P.*(q_hp*rp + q_hp*sum(m(alpha).*M)*dalpha./(d+P).*ap./(ap + sP(P)) ...
        -q_cp*beta.*sum(M.*a./(a + sM(M)))*dalpha ...
        -mup*P) ;
    fm = @(alpha,beta,P,M) M.*( q_cm*sum(beta.*P*dbeta).*a./(a + sM(M) )...
        - q_hm*m(alpha).*sum(P*dbeta./(d+P).*ap./(ap + sP(P))) - mui.*M );
    elseif(choice == 4)
        %% Competition between AMF and Plant + desnity dependence global
    Gamma_m = @(M) Aa*M*dalpha./(sM(M)+(sM(M)<=0));
    Gamma_p = @(P) Ap*P*dbeta./(sP(P)+(sP(P)<=0));
    fp = @(alpha,beta,P,M) P.*(q_hp*rp + q_hp*sum(m(alpha).*M)*dalpha./(d+P).*ap./(ap + sP(P)) ...
        -q_cp*beta.*sum(M.*a./(a + sM(M)))*dalpha ...
        -mup*sum(P)*beta) ;
    fm = @(alpha,beta,P,M) M.*( q_cm*sum(beta.*P*dbeta).*a./(a + sM(M) )...
        - q_hm*m(alpha).*sum(P*dbeta./(d+P).*ap./(ap + sP(P))) - mui.*sum(M)*dalpha );
end


%% Initial data
P0 = 0.1*rand(1,N_plant);
M0 = 0.1*rand(1,N_AMF);
% M0 = 0.1*(abs(ALPHA)<.1);
t = 0; it = 0; 
tt = t;
dt = 0.1;
Pnew = P0';  PP = P0;  
Mnew = M0'; MM = M0;
%% Explicit scheme
MM_b = sum(M0*dalpha);
MM_d_new = MM./MM_b; MM_d_old = 0;
PP_b = sum(P0*dbeta);
PP_d_new = PP./PP_b; PP_d_old = 0;
while (t<Tf)&&(sum(abs(MM_d_new-MM_d_old))>1e-6)
    Pold = Pnew; Mold = Mnew; MM_d_old = MM_d_new;
    Pnew = (I_beta -  dt*dp*Ad_beta)\(Pold + dt*fp(ALPHA',BETA',Pold,Mold));
    Mnew = (I_alpha -  dt*dm(Mold).*Ad_alpha)\(Mold + dt*fm(ALPHA',BETA',Pold,Mold));
    MM_d_new = Mnew'./sum(Mnew*dalpha);
 
    it = it + 1; t = t+dt;
    tt = [tt;t];
    PP = [PP;Pnew'];
    MM = [MM;Mnew'];
    PP_bt = sum(Pnew*dbeta);
    PP_b = [PP_b;PP_bt];
    MM_bt = sum(Mnew*dalpha);
    MM_b = [MM_b;MM_bt];
    
end
%% ode45 scheme competition
% X0 = [P0,M0];
% [t,X] = ode45(@(t,y) Func_AMF_Plant_evol_alpha_comp_continuous(y),[0,Tf],X0);
% tt = 0:dt:Tf;
% PP = interp1(t,X(:,1),tt);
% MM = interp1(t,X(:,2:end),tt);
% MM_b = sum(MM,2)*dalpha;

%% ode45 scheme no competition
% X0 = [P0,M0];
% [t,X_nc] = ode45(@(t,y) Func_AMF_Plant_evol_alpha_nodisp(y),[0,Tf],X0);
% tt_nc = 0:dt:Tf;
% PP_nc = interp1(t,X(:,1),tt);
% MM_nc = interp1(t,X(:,2:end),tt);
% MM_b_nc = sum(MM,2);

%% Plot biommass 
% PP_b = sum(PP,2)*dx;
figure(1)
clf
hold on
plot(tt,PP_b,'-')
plot(tt,MM_b,'-')
% xlim([0,Tf])
drawnow
hold off
legend('plant','AMF')

%% Plot of P and M distribution over time over space trait
MM_d = MM./MM_b;
PP_d = PP./PP_b;
% MM_d_half = MM(:,ALPHA>=0)./sum(MM(:,ALPHA>=0)*dalpha,2);
figure(2)
clf

for i = 1:100:length(tt)
    yyaxis right
    plot(ALPHA,MM_d(end,:))
    % plot(ALPHA(ALPHA>=0),MM_d_half(end,:))
    hold on
    plot(ALPHA,MM_d(i,:))
    drawnow
    ylabel('AMF density distribution')
    hold off
    
    yyaxis left
    plot(BETA,PP_d(end,:))
    % plot(ALPHA(ALPHA>=0),MM_d_half(end,:))
    hold on
    plot(BETA,PP_d(i,:))
    drawnow
        ylabel('Plant density distribution')
    hold off
    % pause
end
% ylim([0,3e-3])

% for It = 1:10:length(tt)
%     figure(3)
%     clf
%     hold on
%     plot(xx,MM_d(It,:),'-')
%     drawnow
%     pause(0.1)
%     hold off
%     
% %     figure(5)
% %     clf
% %     yyaxis left
% %     Fm = fm(PP(It,:),MM(It,:));
% %     plot(xx,Fm)
% % %     axis([xmin,xmax,0,max(Fm,1)])
% %     yyaxis right
% %     plot(xx,MM(It,:))
% %     axis([xmin,xmax,0,1.01*mc_sstar])
% %     drawnow
% end





