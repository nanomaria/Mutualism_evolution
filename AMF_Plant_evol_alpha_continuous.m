clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Tf = 600;
%% Parameter of the model
global q_hp q_cm q_hm q_cp mup mui d rp Aa a Ad_alpha ALPHA dm BETA dalpha
q_hp  = 3; % q>1
q_cm  = 2;
q_hm  = 1;
q_cp  = 1;
%% Trait beta
BETA = 0.4;
beta = BETA;
Nbeta = length(BETA);
% alpha = beta;
%% Trait alpha
alphamin = 0;
alphamax = 10;
dalpha = 0.01;
ALPHA  = alphamin:dalpha:alphamax;
Nalpha = length(ALPHA);

mup = 1/100; % 1/100 %0.3
mui = 1/20; % 1/20 % 0.03

d = 1.2;
rp = 0;
%% Competition terms
a = 0.2;

N_AMF = Nalpha;
N_plant = Nbeta;
Aa = a*(ones(N_AMF,N_AMF)-diag(ones(1,N_AMF)));

% Diffusion matrix alpha 
e = ones(Nalpha,1);
I_alpha  = spdiags(e,0,Nalpha,Nalpha);
Ad_alpha = spdiags([e -2*e e],-1:1,Nalpha,Nalpha);
Ad_alpha(1,1) = -1;
% Ad_alpha(end,end) = -1;
Ad_alpha = Ad_alpha/(dalpha^2);

dm = 1e-1;  % mutation rate AMF

%% Functionnal response and interactions MARIA
choice = 1;
sM = @(M) sum(M)*dalpha;
m = @(alpha) abs(alpha);
if (choice == 1)
    %% Competition between AMF
    Gamma = @(M) Aa*M*dalpha./(sM(M)+(sM(M)<=0));
    fp = @(alpha,P,M) P.*(q_hp*rp + q_hp*sum(m(alpha).*M)*dalpha./(d+P) ...
        -q_cp*beta.*sum(M.*Gamma(M)./(Gamma(M) + sM(M) + ((Gamma(M) + sM(M))<=0) ))*dalpha ...
        -mup*P) ;
    fm = @(alpha,P,M) M.*( (q_cm*beta*Gamma(M)./(Gamma(M) + sM(M) + ((Gamma(M) + sM(M))<=0) )...
        - q_hm*m(alpha)./(d+P)).*P - mui.*sM(M) );
elseif (choice == 2)
    %% No competition
    fp = @(alpha,P,M) P.*(q_hp*rp + q_hp*sum(alpha.*M*dalpha)./(d+P) - q_cp*beta*sum(M)*dalpha ...
        -mup*P) ;
    fm = @(alpha,P,M) M.*( (q_cm*beta - q_hm*m(alpha)./(d+P)).*P - mui.*sM(M) );
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
MM_b = sum(M0);
MM_d_new = MM./MM_b; MM_d_old = 0;
while (t<Tf)&&(sum(abs(MM_d_new-MM_d_old))>1e-6)
    Pold = Pnew; Mold = Mnew; MM_d_old = MM_d_new;
    Pnew = Pold + dt*fp(ALPHA',Pold,Mold);
    Mnew = (I_alpha -  dt*dm*Ad_alpha)\(Mold + dt*fm(ALPHA',Pold,Mold));
    MM_d_new = Mnew'./sum(Mnew*dalpha);
 
    it = it + 1; t = t+dt;
    tt = [tt;t];
    PP = [PP;Pnew];
    MM = [MM;Mnew'];
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
plot(tt,PP,'-')
plot(tt,MM_b,'-')
% xlim([0,Tf])
drawnow
hold off
legend('plant','AMF')

%% Plot of M distribution over time over space trait
MM_d = MM./MM_b;
% MM_d_half = MM(:,ALPHA>=0)./sum(MM(:,ALPHA>=0)*dalpha,2);
figure(2)
clf

for i = 1:10:length(tt)
    plot(ALPHA,MM_d(end,:))
    % plot(ALPHA(ALPHA>=0),MM_d_half(end,:))
    hold on
    plot(ALPHA,MM_d(i,:))
    drawnow
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





