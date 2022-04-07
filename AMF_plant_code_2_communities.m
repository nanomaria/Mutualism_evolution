clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Michorizer model with frequency dependence and evolution of alpha no space   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Time parameter
Tf = 1000;
%% Parameter of the model
global q_hp q_cm q_hm q_cp mup mui d rp Aa1 Aa2 Ad_alpha1 Ad_alpha2  ALPHA1 ALPHA2 dm1 dm2 BETA1 BETA2 Ap1 Ap2 Ad_beta1 Ad_beta2 dp1 dp2 dalpha1 dalpha2 dbeta1 dbeta2 C1 C2
q_hp  = 3; % q>1
q_cm  = 2;
q_hm  = 1;
q_cp  = 1;

%% Trait beta
beta_min1 = 0;
beta_max1 = 1;
dbeta1 = 0.1;
BETA1 = beta_min1:dbeta1:beta_max1;
Nbeta1 = length(BETA1);
beta_min2 = 0;
beta_max2 = 1;
dbeta2 = 0.1;
BETA2 = beta_min2:dbeta2:beta_max2;
Nbeta2 = length(BETA2);

%% Trait alpha
alphamin1 = 0;
alphamax1 = 1;
dalpha1 = 0.1;
ALPHA1  = alphamin1:dalpha1:alphamax1;
Nalpha1 = length(ALPHA1);

alphamin2 = 0;
alphamax2 = 1;
dalpha2 = 0.1;
ALPHA2  = alphamin2:dalpha2:alphamax2;
Nalpha2 = length(ALPHA2);


mup = 0.3; % 1/100 %0.3
mui = 0.3; % 1/20 % 0.03

d = 1.2;

rp = 0;
%% Competition terms
a1_min = 0.2;
a1_max = 0.3;
a1 = @(x) (a1_max-a1_min).*(1-exp(-abs(x)))+a1_min;
ap1_min = 0.2;
ap1_max = 0.3;
ap1 = @(x) (ap1_max-ap1_min).*(1-exp(-abs(x)))+ap1_min;

a2_min = 0.2;
a2_max = 0.3;
a2 = @(x) (a2_max-a2_min).*(1-exp(-abs(x)))+a2_min;
ap2_min = 0.2;
ap2_max = 0.3;
ap2 = @(x) (ap2_max-ap2_min).*(1-exp(-abs(x)))+ap2_min;

C1 = 0.1;
C2 = 0.1;

[alpha1_x,alpha1_y]=meshgrid(ALPHA1,ALPHA1);
AA1 = a1(min(alpha1_x,alpha1_y));
Aa1 = AA1-diag(diag(AA1));
[beta1_x,beta1_y]=meshgrid(BETA1,BETA1);
AP1 = ap1(min(beta1_x,beta1_y));
Ap1 = AP1-diag(diag(AP1));

[alpha2_x,alpha2_y]=meshgrid(ALPHA2,ALPHA2);
AA2 = a2(min(alpha2_x,alpha2_y));
Aa2 = AA2-diag(diag(AA2));
[beta2_x,beta2_y]=meshgrid(BETA2,BETA2);
AP2 = ap2(min(beta2_x,beta2_y));
Ap2 = AP2-diag(diag(AP2));




N_AMF1 = Nalpha1;
N_plant1 = Nbeta1;
N_AMF2 = Nalpha2;
N_plant2 = Nbeta2;


% Diffusion matrix alpha 
e1 = ones(Nalpha1,1);
I_alpha1  = spdiags(e1,0,Nalpha1,Nalpha1);
Ad_alpha1 = spdiags([e1 -2*e1 e1],-1:1,Nalpha1,Nalpha1);
Ad_alpha1(1,1) = -1;
Ad_alpha1(end,end) = -1;
Ad_alpha1 = Ad_alpha1/(dalpha1^2);

e2 = ones(Nalpha1,1);
I_alpha2  = spdiags(e2,0,Nalpha2,Nalpha2);
Ad_alpha2 = spdiags([e2 -2*e2 e2],-1:1,Nalpha2,Nalpha2);
Ad_alpha2(1,1) = -1;
Ad_alpha2(end,end) = -1;
Ad_alpha2 = Ad_alpha2/(dalpha2^2);

% Diffusion matrix beta 
e1 = ones(Nbeta1,1);
I_beta1  = spdiags(e1,0,Nbeta1,Nbeta1);
Ad_beta1 = spdiags([e1 -2*e1 e1],-1:1,Nbeta1,Nbeta1);
Ad_beta1(1,1) = -1;
Ad_beta1(end,end) = -1;
Ad_beta1 = Ad_beta1/(dbeta1^2);

e2 = ones(Nbeta2,1);
I_beta2  = spdiags(e2,0,Nbeta2,Nbeta2);
Ad_beta2 = spdiags([e2 -2*e2 e1],-1:1,Nbeta2,Nbeta2);
Ad_beta2(1,1) = -1;
Ad_beta2(end,end) = -1;
Ad_beta2 = Ad_beta2/(dbeta2^2);


dm1 = 0.1;  % mutation rate AMF
dp1 = 0.1;  % mutation rate plant
dm2 = 0.1;  % mutation rate AMF
dp2 = 0.1;  % mutation rate plant

P01 = 0.1*rand(1,N_plant1);
M01 = 0.1*rand(1,N_AMF1);
P02 = 0.1*rand(1,N_plant2);
M02 = 0.1*rand(1,N_AMF2);

t = 0; it = 0; 
tt = t;
dt = 0.1;
Pnew1 = P01';  PP1 = P01;  
Mnew1 = M01'; MM1 = M01;
Pnew2 = P02';  PP2 = P02;  
Mnew2 = M02'; MM2 = M02;


%% ode45 scheme competition
X0 = [P01,M01,P01,M02];
[t,X] = ode45(@(t,y) Func_AMF_Plant_evol_alpha_beta_comp_nodisp_comm(y),[0,Tf],X0);
tt_disc = t;

PP_disc1 = X(:,1:N_plant1);
PP_b_disc1 = sum(PP_disc1,2);
MM_disc1 = X(:,N_plant1+1:(N_plant1+N_AMF1));
MM_b_disc1 = sum(MM_disc1,2);

PP_disc2 = X(:,(N_plant1+N_AMF1+1):(N_plant1+N_AMF1+N_plant2));
PP_b_disc2 = sum(PP_disc2,2);
MM_disc2 = X(:,(N_plant1+N_AMF1+N_plant2+1):(N_plant1+N_AMF1+N_plant2+N_AMF2));
MM_b_disc2 = sum(MM_disc2,2);


%% ode45 scheme competition continuous
X0 = [P01,M01,P02,M02];

%% Figures
Color = get(gca,'colororder');
Marker = ['o','0','d','^',];
%% Plot biommass 
% PP_b = sum(PP,2)*dx;

figure(1)
clf
hold on
plot(tt_disc,PP_b_disc1,'--','color',Color(1,:))
plot(tt_disc,MM_b_disc1,'-o','color',Color(1,:))
plot(tt_disc,PP_b_disc2,'--','color',Color(2,:))
plot(tt_disc,MM_b_disc2,'-o','color',Color(2,:))

% plot(tt_cont,MM_b_cont,'-o','color',Color(2,:))
% xlim([0,Tf])
xlabel('time $t$','interpreter','latex')
ylabel('biomass of plant and AMF','interpreter','latex')
legend('plant 1', 'AMF 1', 'plant 2', ' AMF 2')
drawnow
hold off

%% Plot of M distribution over time over space trait
PP_d_disc1 = PP_disc1./(PP_b_disc1+PP_b_disc2);
MM_d_disc1 = MM_disc1./(MM_b_disc1+MM_b_disc2);
PP_d_disc2 = PP_disc2./(PP_b_disc1+PP_b_disc2);
MM_d_disc2 = MM_disc2./(MM_b_disc1+MM_b_disc2);


%% Mean trait over time
alpha_bar_disc = sum(ALPHA1.*MM_d_disc1+ALPHA2.*MM_d_disc2,2);
beta_bar_disc = sum(BETA1.*PP_d_disc1+BETA2.*PP_d_disc2,2);
%alpha_bar_disc2 = sum(ALPHA2.*MM_d_disc2,2);
%beta_bar_disc2 = sum(BETA2.*PP_d_disc2,2);

% alpha_bar_cont = sum(ALPHA.*MM_d_cont*dalpha,2);
% beta_bar_cont = sum(BETA.*PP_d_cont*dbeta,2);
figure(3)
clf
hold on
plot(tt_disc,beta_bar_disc,'--','color',Color(1,:))
plot(tt_disc,alpha_bar_disc,'o','color',Color(2,:))

% plot(tt_cont,beta_bar_cont,'--','color',Color(2,:))
% plot(tt_cont,alpha_bar_cont,'o','color',Color(2,:))
ylabel('mean trait in the AMF and plant','interpreter','latex')
xlabel('time $t$','interpreter','latex')
legend('mean trait in plant','mean trait in AMF')

figure(2)
clf
hold on
        plot(BETA1,PP_d_disc1(1,:),'--','linewidth',1,'Color','k')
        plot(BETA1,PP_d_disc1(end,:),'-o','linewidth',2, 'Color','b')
        plot(ALPHA1,MM_d_disc1(end,:),'o-','linewidth',2,'MarkerIndices',1:10:length(ALPHA1),'Color','r')
        plot(ALPHA1,MM_d_disc1(1,:),'--','linewidth',1,'MarkerIndices',1:10:length(ALPHA1),'Color','r')
        plot(BETA2,PP_d_disc2(end,:),'-*','linewidth',2, 'Color','g')
        plot(ALPHA2,MM_d_disc2(end,:),'*-','linewidth',2,'MarkerIndices',1:10:length(ALPHA2),'Color','c')
        plot(BETA2,PP_d_disc2(1,:),'--','linewidth',1,'Color','g')
        plot(ALPHA2,MM_d_disc2(1,:),'--','linewidth',1,'MarkerIndices',1:10:length(ALPHA2),'Color','c')

        xlabel('trait value')
        ylabel('proportion of plants/AMF with that trait')
legend('initial mean trait (all the same)', 'plant 1 end','AMF 1 end', 'plant 2 end', 'AMF 2 end')

