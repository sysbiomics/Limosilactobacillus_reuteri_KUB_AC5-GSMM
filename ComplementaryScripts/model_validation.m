%% Validation of iTN656 
% Purpose: to validate iTN656 by comparing in silico 
% and in vivo growth of L. reuteri KUB-AC5 on different carbon sources 
% i.e. glucose, sucrose, maltose and lactose
% 
% Written by Nachon Raethong, 5-DEC-2021
%% WORKSPACE
cd '/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM';
%initCobraToolbox;
model_final=importModel('model/xml/iTN656.xml',false);
new7maxsOUT=model_final;

%% Constraints
%  The simulation was constrainted according to the experiment:
%  (i)   Anaerobic condition 
%  (ii)  Complex media composed of a complement of amino acids, vitamins,
%        lipids and ions
%  (iii) Used a single carbon substrate i.e. glucose, sucrose, maltose and lactose

amino_acids = {'EXC_BOTH_arg__L_e', 'EXC_BOTH_cys__L_e', 'EXC_BOTH_ile__L_e', 'EXC_BOTH_leu__L_e', 'EXC_BOTH_lys__L_e', 'EXC_BOTH_met__L_e', 'EXC_BOTH_thr__L_e', 'EXC_BOTH_tyr__L_e', 'EXC_BOTH_val__L_e', 'EXC_BOTH_ala__L_e', 'EXC_BOTH_asn__L_e', 'EXC_BOTH_phe__L_e', 'EXC_BOTH_trp__L_e', 'EXC_BOTH_pro__L_e', 'EXC_BOTH_gln__L_e', 'EXC_BOTH_asp__L_e', 'EXC_BOTH_gly_e', 'EXC_BOTH_ser__L_e', 'EXC_BOTH_glu__L_e', 'EXC_BOTH_his__L_e'};
for i = 1:numel(amino_acids)
    new7maxsOUT=setParam(new7maxsOUT,'lb',amino_acids(i),[-5]);
end
lipid= {'EXC_BOTH_hdcea_e', 'EXC_BOTH_ocdcya_e', 'EXC_BOTH_ocdctr_e'};
for i = 1:numel(lipid)
    new7maxsOUT=setParam(new7maxsOUT,'lb',lipid(i),[-5]);
end

vitamin = {'EXC_BOTH_4abz_e', 'EXC_BOTH_pydam_e', 'EXC_BOTH_thm_e', 'EXC_BOTH_pydam_e', 'EXC_BOTH_pnto__R_e', 'EXC_BOTH_nac_e', 'EXC_BOTH_btn_e'};
for i = 1:numel(vitamin)
    new7maxsOUT=setParam(new7maxsOUT,'lb',vitamin(i),[-0.0001]);
end

carbon_sources = {'EXC_BOTH_melib_e', 'EXC_BOTH_raffin_e',...
    'EXC_BOTH_glc__D_e','EXC_BOTH_13ppd_e',...
    'EXC_BOTH_glcn__D_e', 'EXC_BOTH_glyc_e', 'EXC_BOTH_lcts_e',...
    'EXC_BOTH_sucr_e', 'EXC_BOTH_tre_e', ...
    'EXC_BOTH_rib__D_e',...
    'EXC_BOTH_gal_e', 'EXC_BOTH_malt_e', 'EXC_BOTH_arab__L_e',...
    'EXC_BOTH_tre_e', 'EXC_BOTH_sucr_e',...
    'EXC_BOTH_raffin_e', 'EXC_BOTH_melib_e',...
    'EXC_BOTH_malt_e', 'EXC_BOTH_lcts_e',...
    'EXC_BOTH_fuc__L_e', 'EXC_BOTH_arab__L_e', 'EXC_BOTH_glyc_e',...
    'EXC_BOTH_rib__D_e', 'EXC_BOTH_mnl_e',...
    'EXC_BOTH_glc__D_e', 'EXC_BOTH_glcn__D_e', 'EXC_BOTH_gal_e',...
    'EXC_BOTH_fru_e'};

for i = 1:numel(carbon_sources)
    new7maxsOUT=setParam(new7maxsOUT,'lb',carbon_sources(i),[0]);
end


products = {'EXC_BOTH_hxan_e' 'EXC_BOTH_xan_e',...
    'EXC_BOTH_gcald_e', 'EXC_BOTH_btd__RR_e',...
    'EXC_BOTH_ac_e', 'EXC_BOTH_lac__D_e',...
    'EXC_BOTH_orot_e','EXC_BOTH_etoh_e', ...
    'EXC_BOTH_mal__L_e', 'EXC_BOTH_lac__L_e'};
for i = 1:numel(products)
    new7maxsOUT=setParam(new7maxsOUT,'lb',products(i),[0]);
    new7maxsOUT=setParam(new7maxsOUT,'ub',products(i),[1000]);
end

new7maxsOUT2 = setParam(new7maxsOUT,'lb',{'EXC_BOTH_o2_e'},[-0.000000000001]); 

%% The simulation
%  The growth simulation was carried out by constrainted the model into 4 conditions
%  accroding to the measured uptake rate of the each carbon substrate 
%  i.e. glucose, sucrose, maltose and lactose

new7maxsOUT3 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_glc__D_e'},[-6.634]); 
new7maxsOUT4 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_sucr_e'},[-5.033]); 
new7maxsOUT5 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_lcts_e'},[-0.940]); 
new7maxsOUT6 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_malt_e'},[-3.310]); 

sol3 = solveLP(new7maxsOUT3);
sol4 = solveLP(new7maxsOUT4);
sol5 = solveLP(new7maxsOUT5);
sol6 = solveLP(new7maxsOUT6);
fprintf(['umax glucose = ' num2str(sol3.f*-1) ' per hour' '\n']);
fprintf(['umax sucrose = ' num2str(sol4.f*-1) ' per hour' '\n']);
fprintf(['umax lactose = ' num2str(sol5.f*-1) ' per hour' '\n']);
fprintf(['umax maltose = ' num2str(sol6.f*-1) ' per hour' '\n']);


% in vivo
%EXC_BOTH_glc__D_e 6.634±0.684 > umax 0.151±0.004 
%EXC_BOTH_lcts_e 0.940±0.322 > umax 0.078±0.005 
%EXC_BOTH_sucr_e 5.033±0.310 > umax 0.247±0.003 c%EXC_BOTH_malt_e 3.310±0.764 > umax 0.199±0.009 

%in silico
%umax glucose = 0.15237 per hour
%umax sucrose = 0.26709 per hour
%umax lactose = 0.076945 per hour
%umax maltose = 0.18302 per hourc
	%in vivo	in silico	%ERROR
%Glucose	0.151	0.152	0.899
%Lactose	0.078	0.077	1.371
%Sucrose	0.247	0.274	8.133
%Maltose	0.199	0.183	8.731
