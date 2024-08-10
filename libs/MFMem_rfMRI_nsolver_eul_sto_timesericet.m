function [BOLD_TR,neuro_act,yT] = MFMem_rfMRI_nsolver_eul_sto_timesericet(parameter,SC,Nstate,Tepochlong,TR,Si,Sti_time)
% function [BOLD_TR,neuro_act,yT] = MFMem_rfMRI_nsolver_eul_sto_timeserice(parameter,SC,Nstate,Tepochlong,TR,I)

%-----------------------------------------------------------------------------
% [FC_sim, RSN_FC_sim,metastable_sim,synchrony_sim] = MFMem_rfMRI_nsolver_eul_sto(parameter,SC,Nstate,Tepochlong,TR)
%
% Function to 
%  (a)solve diffitial equation of dynamic mean field and hemodynamic model using stochastic Euler method 
%    
%  (b)caculate simulated functional connectivity (FC) from simulated BOLD
%  (c)caculate simulated RSN_FC
%
% Input:
%     - SC:        structural connectivity matrix
%     - Nstate:    noise randon seed
%     - Tepochlong:simulation long in [min], exclusive 1min
%     pre-simulation(casted)  仿真长度，不包括1分钟的初始状态
% Output:
%     - BOLD_d:  simulated bold
%   
%
%----------------------------------------------------------------------------

%caculate first [BOLD,yT,fT,qT,vT,zT,Time], then caculate FC_stim


%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
%(a) solve diffitial equation of dynamic mean field and hemodynamic model using stochastic Euler method 
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if size(parameter, 2) > 1
    error('Input argument ''parameter'' should be a column vector');
end


%-----------------------------------------------------------------------
%initial system
%-----------------------------------------------------------------------
%simulation time
kstart = 0;  %s
Tpre = 60; %s
kend = Tpre+60*Tepochlong; %s

dt_l = 0.001;    %s  integration time step

dt = dt_l;  %s  time step for neuro
dtt = 0.001; %s, time step for BOLD

%sampling ratio
k_P = kstart:dt:kend;   
k_PP = kstart:dtt:kend;

%initial
Nnodes = size(SC,1);
Nsamples = length(k_P);
Bsamples = length(k_PP);

%for neural activity y0 = 0
yT = zeros(Nnodes,1);

%for hemodynamic activity z0 = 0, f0 = v0 = q0 =1
zT = zeros(Nnodes,Bsamples);
fT = zeros(Nnodes,Bsamples);
fT(:,1) = 1;
vT = zeros(Nnodes,Bsamples);
vT(:,1) = 1;
qT = zeros(Nnodes,Bsamples);
qT(:,1) = 1;

F = [zT(:,1) fT(:,1) vT(:,1) qT(:,1)];
yT(:,1) = 0.001;


%wiener process
w_coef = parameter(end)/sqrt(0.001); 
w_dt = dt; %s
w_L = length(k_P);
rng(Nstate);% rng(seed) 使用非负整数 seed 为随机数生成器提供种子
dW = sqrt(w_dt)*randn(Nnodes,w_L+1000);  %plus 1000 warm-up

j = 0;

%--------------------------------------------------------------------------
%solver: Euler
%-------------------------------------------------------------------------- 
tic

%% warm-up


for i = 1:1000
        %%此处属于预热阶段
        dy = MFMem_rfMRI_mfm_ode1(yT,parameter,SC);
        yT = yT + dy*dt + w_coef*dW(:,i);

end


%% main body: caculation 

for i = 1:length(k_P)
        
        
        dy = MFMem_rfMRI_mfm_ode1(yT,parameter,SC); 
        if( 60 > k_P(i))
              yT = yT + dy*dt + w_coef*dW(:,i+1000);
        elseif  k_P(i) <= Sti_time + 60
            yT = yT + dy*dt + w_coef*dW(:,i+1000)+Si;
        else
             yT = yT + dy*dt + w_coef*dW(:,i+1000);
        end

        if mod(i,dtt/dt) == 0
            j = j+1;
            y_neuro(:,j) = yT;
        end
        
end

for i = 2:length(k_PP)

            dF = MFMem_rfMRI_rfMRI_BW_ode1(y_neuro(:,i-1),F,Nnodes);
            F = F + dF*dtt;
            zT(:,i) = F(:,1);
            fT(:,i) = F(:,2);
            vT(:,i) = F(:,3);
            qT(:,i) = F(:,4);

end


%Parameter for Balloom-Windkessel model, we updated  this model and its
%parameter according to Stephan et al 2007, NeuroImage 38:387-401 and
%Heinzle et al. 2016 NeuroImage 125:556-570, Parameter are set for the 3T
%and TE=0.0331s
p = 0.34; 
v0 = 0.02;
k1 = 4.3*28.265*3*0.0331*p;
k2 = 0.47*110*0.0331*p;
k3 = 1-0.47;
y_BOLD = 100/p*v0*( k1*(1-qT) + k2*(1-qT./vT) + k3*(1-vT) );     

Time = k_PP;
toc
 


%--------------------------------------------------------------------------
%(b)&(c) compute simulated FC and correlation of 2 FCs
%--------------------------------------------------------------------------
%get the static part
cut_indx = find(Time == Tpre);% after xx s
BOLD_act = y_BOLD(:,cut_indx:end);
neuro_act = y_neuro(:,cut_indx:end);
Time_act = Time(:,cut_indx:end);
[BOLD_TR]=rfMRI_simBOLD_downsampling(BOLD_act,TR/dt); %down sample 
end


