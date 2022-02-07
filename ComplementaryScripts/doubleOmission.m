%% Validation of iTN656 
% Purpose: to validate iTN656 by comparing in silico 
% and in vivo growth of L. reuteri KUB-AC5 on different carbon sources 
% i.e. glucose, sucrose, maltose and lactose
% 
% Written by Nachon Raethong, 5-DEC-2021
%% WORKSPACE
cd '/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM';
initCobraToolbox;
model_final=importModel('model/xml/iTN656.xml',false);
new7maxsOUT=model_final;
new7maxsOUT

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
nutrient_aa_l = union(amino_acids,lipid);
nutrients=union(nutrient_aa_l,vitamin);
carbon_sources = {'EXC_BOTH_melib_e', 'EXC_BOTH_raffin_e',...
    'EXC_BOTH_glc__D_e','EXC_BOTH_13ppd_e','EXC_BOTH_pyr_e',...
    'EXC_BOTH_glcn__D_e', 'EXC_BOTH_glyc_e', 'EXC_BOTH_lcts_e',...
    'EXC_BOTH_sucr_e', 'EXC_BOTH_tre_e', 'EXC_BOTH_acald_e',...
    'EXC_BOTH_rib__D_e',...
    'EXC_BOTH_gal_e', 'EXC_BOTH_malt_e', 'EXC_BOTH_arab__L_e',...
    'EXC_BOTH_tre_e', 'EXC_BOTH_sucr_e','EXC_BOTH_drib_e',...
    'EXC_BOTH_raffin_e', 'EXC_BOTH_melib_e',...
    'EXC_BOTH_malt_e', 'EXC_BOTH_lcts_e',...
    'EXC_BOTH_fuc__L_e', 'EXC_BOTH_arab__L_e', 'EXC_BOTH_glyc_e',...
    'EXC_BOTH_rib__D_e', 'EXC_BOTH_mnl_e',...
    'EXC_BOTH_glc__D_e', 'EXC_BOTH_glcn__D_e', 'EXC_BOTH_gal_e',...
    'EXC_BOTH_fru_e','EXC_BOTH_actn__R_e'};

for i = 1:numel(carbon_sources)
    new7maxsOUT=setParam(new7maxsOUT,'lb',carbon_sources(i),[0]);
end


products = {'EXC_BOTH_hxan_e' 'EXC_BOTH_xan_e',...
    'EXC_BOTH_succ_e','EXC_BOTH_gcald_e', 'EXC_BOTH_btd__RR_e',...
    'EXC_BOTH_Levan', 'EXC_BOTH_ac_e', 'EXC_BOTH_lac__D_e',...
    'EXC_BOTH_orot_e','EXC_BOTH_etoh_e', 'EXC_BOTH_diact_e',...
    'EXC_BOTH_mal__L_e', 'EXC_BOTH_lac__L_e'};
for i = 1:numel(products)
    new7maxsOUT=setParam(new7maxsOUT,'lb',products(i),[0]);
    new7maxsOUT=setParam(new7maxsOUT,'ub',products(i),[1000]);
end

new7maxsOUT2 = setParam(new7maxsOUT,'lb',{'EXC_BOTH_o2_e'},[-0.000000000001]); 


for i = 1:numel(nutrients)
    new7maxsOUT3=setParam(new7maxsOUT2,'lb',nutrients(i),[-1]);
    new7maxsOUT3=setParam(new7maxsOUT3,'ub',nutrients(i),[1000]);
end

new7maxsOUT3 = setParam(new7maxsOUT3,'lb',{'EXC_BOTH_glc__D_e'},[-25]); 

%% Omission test
omiModel=new7maxsOUT3;
solomiModel = solveLP(omiModel,1);
fprintf(['umax omiModel = ' num2str(solomiModel.f*-1) ' per hour' '\n']);

mVersatile = cell(numel(nutrients),numel(nutrients));
row=cell(numel(nutrients),1);
column=cell(numel(nutrients),1);
for i = 1:numel(nutrients)
    model = omiModel;
    row{i} = {nutrients{i}};
    model_i = setParam(model,'lb',nutrients{i},[0]);
    for j = 1:numel(nutrients)
        column{j} = {nutrients{j}};
        model_ij = setParam(model_i,'lb',nutrients{j},[0]);
        solmodel_ij = solveLP(model_ij,1);
        mVersatile(i,j) = {(solmodel_ij.f*-1)};
    end
end


