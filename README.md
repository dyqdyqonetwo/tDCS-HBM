
We present a novel tDCS hybrid brain model (tDCS-HBM) that incorporates tDCS-induced gray matter electric fields into a large-scale brain network model, considering their relationship with membrane potential to effectively predict spatiotemporal dynamics. Using this model, we simulated brain activity in response to tDCS over the left (anode F3, cathode Fp2) and right DLPFC (anode F4, cathode Fp1).

Step1_tDCS_HBM_V.m is used to simulate brain activity in response to tDCS.

Step1_Without_stimulation.m is used to simulate brain activity without receiving tDCS stimulation.

data/FEM_E/F3a-Fp2c.mat refers to the normal electric field under the F3-Fp2 montage modelled by tDCS finite elements.

data/FEM_E/F4a-Fp1c.mat refers to the normal electric field under the F4-Fp1 montage modelled by tDCS finite elements.

data/MFM_para/FCSC_Desikan68_ave_2_para4_3.mat contains the structural connectivity matrix (SC), the functional connectivity matrix (FC_emp), and the mean-field model (para_E) trained parameters.
