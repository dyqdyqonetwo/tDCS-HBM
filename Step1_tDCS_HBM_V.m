clear all
close all


%% setup directory for lib, data, cluster, save
main_dir = pwd;
data_dir = fullfile(main_dir,'data'); 
save_dir_input = fullfile(main_dir,'result_tDCS_F3a_Fp2c');
save_dir = save_dir_input;

lib_dir = fullfile(pwd, 'libs');
addpath(lib_dir);


load([data_dir,'\','MFM_para\FCSC_Desikan68_ave_2_para4_3.mat'],'SC','Para_E','FC_emp');
SC = SC./max(max(SC)).*0.2;
Tmax=14.4;
TR=0.72;
%load F3a-Fp2c  E_normal
Sti_time = 60
V = load([data_dir,'\FEM_E\F3a_Fp2c.mat']).V
Si =0.13 * V 
mkdir(save_dir)
for trials=1:30
    Nstate = rng;
    [BOLD_TR] = MFMem_rfMRI_nsolver_eul_sto_timesericet(Para_E,SC,Nstate,Tmax,TR,Si,Sti_time);
    [FC_sim,metastable_sim,synchrony_sim,FC_cor] = estimation_corr_emp_sim_noRSN(FC_emp,BOLD_TR);
    save( [save_dir ,'\',num2str(trials),'.mat'] ,'FC_sim','metastable_sim','synchrony_sim','FC_cor','BOLD_TR');  
end

rmpath(lib_dir);


