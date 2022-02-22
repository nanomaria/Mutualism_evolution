clear all
close all
% format LongE;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       2 plants with same community for each plant                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Time parameter
Tf = 1000;
%% Parameter of the model
global q_hp q_cm q_hm q_cp mup mui d rp Aa ALPHA BETA Ap  dbeta
q_hp  = 3; % q>1
q_cm  = 2;
q_hm  = 1;
q_cp  = 1;
%% Trait beta
% beta  = 0.4;
beta_min = 0.2;
beta_max = 0.4;
dbeta = 0.1;
N_plant = 2;
BETA = linspace(beta_min,beta_max,N_plant); %[beta_min,(beta_min+beta_max)/2,beta_max];  % beta_min:dbeta:beta_max;
Nbeta = length(BETA);
N_plant = Nbeta;

%% Trait alpha
alphamin = beta_min;
alphamax = 10;
Nalpha = 5;
ALPHA  = linspace(alphamin,alphamax,Nalpha); % alphamin:dalpha:alphamax;
N_AMF = Nalpha;

mup = 0.03; % 1/100 %0.3
mui = 0.05; % 1/20 % 0.03

d = 1.2;

rp = 0;

%% Competition terms
a = 0.02;
Aa = a*(ones(N_AMF,N_AMF)-diag(ones(1,N_AMF)));

ap = 0.02;
Ap = ap; %*(ones(N_plant,N_plant)-diag(ones(1,N_plant)));
%% Initial data
%% Random
% P0 = zeros(1,N_plant); 
% Ip = 1;
% P0(Ip) = 1;
% M0 = zeros(1,N_plant*N_AMF);
% M0(1+(Ip-1)*N_AMF:Ip*N_AMF) = 0.5*rand(1,N_AMF);

P0 =  10*rand(1,N_plant); %[1.7,0.3];
M0 = rand(1,N_AMF*N_plant);

% P0 = .1*(BETA<0.1);
% M0 = [.1*(ALPHA<0.5),.1*(ALPHA<0.5)];
% t = 0; it = 0; 
% tt = t;
% dt = 0.1;
% Pnew = P0';  PP = P0;  
% Mnew = M0'; MM = M0;

%% ode45 scheme competition
X0 = [P0,M0];
[t,X] = ode45(@(t,y) Func_AMF_2Plant_noevol_nodisp(y),[0,Tf],X0);
tt = t;
PP = X(:,1:N_plant);
PP_b = sum(PP,2);
MM = zeros(length(t),N_AMF,N_plant);
for i = 1:N_plant
    MM(:,:,i) = X(:,N_plant+1+(i-1)*N_AMF:N_plant+i*N_AMF);
end
MM_b = permute(sum(MM,2),[1,3,2]);

%% Plot biommass 
% PP_b = sum(PP,2)*dx;
Color = get(gca,'colororder');
Marker = ['o','*','d','^'];
figure(1)
clf
hold on
for i = 1:N_plant
    plot(tt,PP(:,i),'--','color',Color(i,:))
    plot(tt,MM_b(:,i),'-o','color',Color(i,:))
    % xlim([0,Tf])
    drawnow
end
hold off













































