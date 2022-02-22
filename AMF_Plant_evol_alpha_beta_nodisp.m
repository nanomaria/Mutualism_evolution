clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Tf = 1000;
%% Parameter of the model
global q_hp q_cm q_hm q_cp mup mui d rp Aa Ad_alpha ALPHA dm BETA Ap Ad_beta dp dalpha dbeta
q_hp  = 3; % q>1
q_cm  = 2;
q_hm  = 1;
q_cp  = 1;
%% Trait beta
% beta  = 0.4;
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

mup = 0.3; % 1/100 %0.3
mui = 0.03; % 1/20 % 0.03

d = 1.2;

rp = 0;
%% Competition terms
a_min = 0.02;
a_max = 0.2;
a = @(x) (a_max-a_min).*(1-exp(-abs(x)))+a_min;
ap_min = 0.02;
ap_max = 0.2;
ap = @(x) (ap_max-ap_min).*(1-exp(-abs(x)))+ap_min;

[alpha_x,alpha_y]=meshgrid(ALPHA,ALPHA);
AA = a(min(alpha_x,alpha_y));
Aa = AA-diag(diag(AA));
[beta_x,beta_y]=meshgrid(BETA,BETA);
AP = a(min(beta_x,beta_y));
Ap = AP-diag(diag(AP));

N_AMF = Nalpha;
N_plant = Nbeta;
% Aa = a*(ones(N_AMF,N_AMF)-diag(ones(1,N_AMF)));
% Ap = ap*(ones(N_plant,N_plant)-diag(ones(1,N_plant)));


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


dm = 0.01;  % mutation rate AMF
dp = 0.001;  % mutation rate plant
% %% Space x
% xmin = -10;
% xmax = 10;
% dx = 0.1;
% xx = xmin:dx:xmax;
% Nx = length(xx);
% 
% % Diffusion matrix
% e = ones(Nx,1);
% I  = spdiags(e,0,Nx,Nx);
% Ad_x = spdiags([e -2*e e],-1:1,Nx,Nx);
% Ad_x(1,1) = -1;
% Ad_x(end,end) = -1;
% Ad_x = Ad_x/(dx^2);

% D_p = 0.1;
% D_m = 0.1;  % diffusion rate
%% Functionnal response and interactions MARIA
% sM = @(M) sum(M.*ones(1,N_AMF) - diag(M))'*dalpha;
% Gamma = @(M) Aa*M./(sM(M)+(sM(M)<=0));
% fp = @(alpha,P,M) P.*(q_hp*rp + q_hp*sum(alpha.*M)*dalpha./(d+P) ...
%     -q_cp*beta*sum(M.*Gamma(M)./(Gamma(M) + sM(M) + ((Gamma(M) + sM(M))<=0) ))*dalpha ...
%     -mup*P) ;
% sM = @(M) sum(M.*ones(1,N_AMF) - diag(M))';
% Gamma = @(M) Aa*M./(sM(M)+(sM(M)<=0));
% fp = @(alpha,P,M) P.*(q_hp*rp + q_hp*sum(alpha.*M)./(d+P) ...
%     -q_cp*beta*sum(M.*Gamma(M)./(Gamma(M) + sM(M) + ((Gamma(M) + sM(M))<=0) )) ...
%     -mup*P) ;
% fm = @(alpha,P,M) M.*( (q_cm*beta*Gamma(M)./(Gamma(M) + sM(M) + ((Gamma(M) + sM(M))<=0) )...
%     - q_hm*alpha./(d+P)).*P - mui.*M );

%% Initial data
%% Random
% P0 = 0.1*rand(1,N_plant);
% M0 = 0.1*rand(1,N_AMF);
P0 = .1*(BETA<0.1);
M0 = .1*(ALPHA<0.5);
t = 0; it = 0; 
tt = t;
dt = 0.1;
Pnew = P0';  PP = P0;  
Mnew = M0'; MM = M0;
%% Explicit sceme
% MM_b = sum(M0);
% MM_d_new = MM./MM_b; MM_d_old = 0;
% while (t<Tf)&&(sum(abs(MM_d_new-MM_d_old))>1e-6)
%     Pold = Pnew; Mold = Mnew; MM_d_old = MM_d_new;
%     Pnew = Pold + dt*fp(ALPHA',Pold,Mold);
%     Mnew = (I -  dt*dm*Ad_alpha)\(Mold + dt*fm(ALPHA',Pold,Mold));
%     MM_d_new = Mnew'./sum(Mnew);
%  
%     it = it + 1; t = t+dt;
%     tt = [tt;t];
%     PP = [PP;Pnew];
%     MM = [MM;Mnew'];
%     MM_bt = sum(Mnew);
%     MM_b = [MM_b;MM_bt];
%     
% end
%% ode45 scheme competition
X0 = [P0,M0];
[t,X] = ode45(@(t,y) Func_AMF_Plant_evol_alpha_beta_comp_nodisp(y),[0,Tf],X0);
tt_disc = t;
PP_disc = X(:,1:N_plant);
PP_b_disc = sum(PP_disc,2);
MM_disc = X(:,N_plant+1:end);
MM_b_disc = sum(MM_disc,2);

%% ode45 scheme competition continuous
X0 = [P0,M0];
% [t,X] = ode45(@(t,y) Func_AMF_Plant_evol_alpha_beta_comp_nodisp_continuous(y),[0,Tf],X0);
% tt = 0:dt:Tf;
% PP = interp1(t,X(:,1:N_plant),tt);
% MM = interp1(t,X(:,N_plant+1:end),tt);
% tt_cont = t;
% PP_cont = X(:,1:N_plant);
% PP_b_cont = sum(PP_cont,2)*dbeta;
% MM_cont = X(:,N_plant+1:end);
% MM_b_cont = sum(MM_cont,2)*dalpha;

%% ode45 scheme no competition
% X0 = [P0,M0];
% [t,X_nc] = ode45(@(t,y) Func_AMF_Plant_evol_alpha_nodisp(y),[0,Tf],X0);
% tt_nc = 0:dt:Tf;
% PP_nc = interp1(t,X(:,1),tt);
% MM_nc = interp1(t,X(:,2:end),tt);
% MM_b_nc = sum(MM,2);

%% Figures
Color = get(gca,'colororder');
Marker = ['o','*','d','^'];
%% Plot biommass 
% PP_b = sum(PP,2)*dx;

figure(1)
clf
hold on
plot(tt_disc,PP_b_disc,'--','color',Color(1,:))
plot(tt_disc,MM_b_disc,'-o','color',Color(1,:))
% plot(tt_cont,PP_b_cont,'--','color',Color(2,:))
% plot(tt_cont,MM_b_cont,'-o','color',Color(2,:))
% xlim([0,Tf])
xlabel('time $t$','interpreter','latex')
ylabel('biomass of plant and AMF','interpreter','latex')

drawnow
hold off

%% Plot of M distribution over time over space trait
PP_d_disc = PP_disc./PP_b_disc;
MM_d_disc = MM_disc./MM_b_disc;
% PP_d_cont = PP_cont./PP_b_cont;
% MM_d_cont = MM_cont./MM_b_cont;
figure(2)
clf
hold on
    plot(BETA,PP_d_disc(end,:),'--','linewidth',2)
        plot(ALPHA,MM_d_disc(end,:),'o-','linewidth',2,'MarkerIndices',1:10:length(ALPHA))
legend('plant','AMF')
% for i = 1:10:length(tt)
% for i = length(tt)
%     subplot(1,2,1)
%     hold on
%     plot(ALPHA,MM_d_disc(end,:),'o-','linewidth',2,'MarkerIndices',1:10:length(ALPHA))
% %     plot(ALPHA,MM_d_cont(1,:),'-')
%     
% %     plot(ALPHA,MM_d_cont(end,:),'o-','linewidth',2,'MarkerIndices',1:10:length(ALPHA))
% %     plot(ALPHA,MM_d_cont(1,:),'-')
%     drawnow
% %     ylim([0,max(MM_d(end,:))])
%     ylabel('AMF distribution','interpreter','latex')
%     xlabel('phosphorous supply $\alpha$','interpreter','latex')
% 
%     hold off
%     
%     subplot(1,2,2)
%     plot(BETA,PP_d_disc(end,:),'--','linewidth',2)
%     hold on
% %     plot(BETA,PP_d(1,:),'--','linewidth',2)
% %     plot(BETA,PP_d_cont(end,:),'--','linewidth',2)
%     
%     drawnow
% %      ylim([0,max(PP_d(end,:))])
%     ylabel('Plant distribution','interpreter','latex')
%     xlabel('carbon supply $\beta$','interpreter','latex')
%     hold off
    % pause
% end

%% Mean trait over time
alpha_bar_disc = sum(ALPHA.*MM_d_disc,2);
beta_bar_disc = sum(BETA.*PP_d_disc,2);
% alpha_bar_cont = sum(ALPHA.*MM_d_cont*dalpha,2);
% beta_bar_cont = sum(BETA.*PP_d_cont*dbeta,2);
figure(3)
clf
hold on
plot(tt_disc,beta_bar_disc,'--','color',Color(1,:))
plot(tt_disc,alpha_bar_disc,'o','color',Color(1,:))

% plot(tt_cont,beta_bar_cont,'--','color',Color(2,:))
% plot(tt_cont,alpha_bar_cont,'o','color',Color(2,:))
ylabel('mean trait in the AMF and plant','interpreter','latex')
xlabel('time $t$','interpreter','latex')

%% Mean fitness of plant and AMF over time
% P = PP_cont';
% M = MM_cont';
% alpha = ALPHA';
% beta  = BETA';
% sM = sum(M*dalpha);
% Gamma_m = Aa*M./(sM+(sM<=0));
% sP = sum(P*dbeta);
% Gamma_p = Ap*P./(sP+(sP<=0));
% fp = (q_hp*rp + q_hp*sum(alpha.*M*dalpha )./(d+P).*Gamma_p./(Gamma_p + sP + ((Gamma_p + sP)<=0)) ...
%     -q_cp*beta.*sum(M.*Gamma_m*dalpha./(Gamma_m + sM + ((Gamma_m + sM)<=0) )) ); % ...
% %     -mup*P);
% fm = ( q_cm*sum(beta.*P*dbeta).*Gamma_m./(Gamma_m + sM + ((Gamma_m + sM)<=0) )...
%     - q_hm*alpha.*sum(P*dbeta./(d+P).*Gamma_p./(Gamma_p + sP + ((Gamma_p + sP)<=0))));%  - mui.*M );
% f_bar_P = sum(fp.*P)./sum(P);
% f_bar_M = sum(fm.*M)./sum(M);
% 
% figure(4)
% clf
% hold on
% plot(tt_cont,f_bar_P,'--','color',Color(2,:))
% plot(tt_cont,f_bar_M,'-o','color',Color(2,:))

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





