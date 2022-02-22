clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Tf = 500;
%% Parameter of the model
global q_hp q_cm q_hm q_cp beta mup mui d rp Aa a
q_hp  = 3; % q>1
q_cm  = 2;
q_hm  = 1;
q_cp  = 1;

beta  = 0.8;
% alpha = beta;

mup = 0.3; % 1/100
mui = 0.03; % 1/20

d = 1.2;

rp = 0.02;
%% Choice case 3
choice = 1;
% if (choice==3)       % Case 3(iii) Invasion
%     a_ww = 2.2;
%     a_cw = 0.1; % 0.001
%     a_wc = 2.2; % 2.2, 0.2
% elseif (choice == 2) % Case 3(ii) Exclusion
%     a_ww = 2.2;
%     a_cw = 0.3;
%     a_wc = 0.3;
% elseif (choice == 1) % Case 3(i) Coexistence
    a = 0.2;
    a_ww = a;
    a_cw = a;
    a_wc = a;
    a_cc = a;
% end

N_AMF = 2;
if (choice==1)
    Aa = a*ones(N_AMF+1,N_AMF+1);
else
    Aa = [0 , a_wc*ones(1,N_AMF);...
        a_cw*ones(N_AMF,1),a_ww*(ones(N_AMF,N_AMF)-diag(ones(1,N_AMF)))];
end

%% Trait alpha
xmin = 0;
xmax = 1;
dx = 1e-2;
xx = xmin:dx:xmax;
Nx = length(xx);

% Diffusion matrix
e = ones(Nx,1);
I  = spdiags(e,0,Nx,Nx);
Ad = spdiags([e -2*e e],-1:1,Nx,Nx);
Ad(1,1) = -1;
Ad(end,end) = -1;
Ad = Ad/(dx^2);

dm = 1e-2;  % mutation rate

%% Functionnal response and interactions MARIA
% fp = @(alpha,P,M)  (q_hp*rp + (q_hp*sum(alpha.*M*dx)./(d+P) -q_cp*beta*a*sum(M*dx)./(a+sum(M*dx)) )-mup*P).*P ;
% fm = @(alpha,P,M)  M.*((q_cm*beta*a./(a+sum(M*dx)) - q_hm*alpha./(d+P)).*P - mui);

%% Competition + maintenance cost
fp = @(alpha,P,M)  (q_hp*rp + (q_hp*sum(alpha.*M*dx)./(d+P) -q_cp*beta*a*sum(M*dx)./(a+sum(M*dx)) )-mup*P).*P ;
fm = @(alpha,P,M)  M.*((q_cm*beta*a./(a+sum(M*dx)) - q_hm*alpha./(d+P)).*P - mui.*sum(M)*dx);

%% No Competition + maintenance cost
% fp = @(alpha,P,M)  (q_hp*rp + (q_hp*sum(alpha.*M*dx)./(d+P) -q_cp*beta*sum(M*dx) )-mup*P).*P ;
% fm = @(alpha,P,M)  M.*((q_cm*beta - q_hm*alpha./(d+P)).*P - mui.*sum(M)*dx);

%% Initial data
P0 = 0.15;
% M0 = (xx<=0.55).*(xx>=0.45);
M0 = (xx<0.01);
it = 0; tt = it;
dt = 0.05;
Pnew = P0'; PP = P0;  
Mnew = M0'; MM = M0;
while (it<Tf)
    Pold = Pnew;Mold = Mnew;
    Pnew   = (Pold + dt*fp(xx',Pold,Mold));
    Mnew   = (I-dt*dm*Ad)\(Mold + dt*fm(xx',Pold,Mold));
    
    it = it + dt; tt = [tt;it];
    PP = [PP;Pnew'];
    MM = [MM; Mnew'];
end

%% Plot biommass 
MM_b = sum(MM,2)*dx;
mean_alpha = sum(xx.*MM,2)./sum(MM,2);
figure(1)
clf
hold on
plot(tt,PP,'--')
plot(tt,MM_b,'-o')
xlim([0,Tf])
drawnow
hold off

%% Plot of M distribution over time over space trait
figure(1)
clf
plot(tt,mean_alpha)
ylabel('Mean trait of AMF $\displaystyle mean(\alpha) =\int_{\overline\alpha}^{\underline\alpha} {\alpha\,m(t,\alpha)\,d\alpha}$','Interpreter', 'latex','FontSize',16)
xlabel('time $t$','Interpreter','latex','FontSize',16)
axis([0,Tf,0,1])


MM_d = MM./MM_b;
figure(2)
clf
plot(xx,MM_d(end,:))
ylabel('Trait distribution of AMF $m(t,\alpha)$','Interpreter', 'latex','FontSize',16)
xlabel('trait $\alpha$','Interpreter','latex','FontSize',16)
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





