%% Simulation of L. reuteri KUB-AC5 growth on different carbon sources
% Purpose: to determine growth capability of L. reuteri KUB-AC5 
% on different carbon sources
% 
% Written by Nachon Raethong, 31-JAN-2019
% Updated by Nachon Raethong, 24-JULY-2019
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
new7maxsOUT = setParam(new7maxsOUT,'lb',{'EXC_BOTH_o2_e'},[-0.000000000001]); 

amino_acids = {'EXC_BOTH_arg__L_e', 'EXC_BOTH_cys__L_e', 'EXC_BOTH_ile__L_e', 'EXC_BOTH_leu__L_e', 'EXC_BOTH_lys__L_e', 'EXC_BOTH_met__L_e', 'EXC_BOTH_thr__L_e', 'EXC_BOTH_tyr__L_e', 'EXC_BOTH_val__L_e', 'EXC_BOTH_ala__L_e', 'EXC_BOTH_asn__L_e', 'EXC_BOTH_phe__L_e', 'EXC_BOTH_trp__L_e', 'EXC_BOTH_pro__L_e', 'EXC_BOTH_gln__L_e', 'EXC_BOTH_asp__L_e', 'EXC_BOTH_gly_e', 'EXC_BOTH_ser__L_e', 'EXC_BOTH_glu__L_e', 'EXC_BOTH_his__L_e'};
for i = 1:numel(amino_acids)
    new7maxsOUT=setParam(new7maxsOUT,'lb',amino_acids(i),[-1]);
end
lipid= {'EXC_BOTH_hdcea_e', 'EXC_BOTH_ocdcya_e', 'EXC_BOTH_ocdctr_e'};
for i = 1:numel(lipid)
    new7maxsOUT=setParam(new7maxsOUT,'lb',lipid(i),[-1]);
end

vitamin = {'EXC_BOTH_4abz_e', 'EXC_BOTH_pydam_e', 'EXC_BOTH_thm_e', 'EXC_BOTH_pydam_e', 'EXC_BOTH_pnto__R_e', 'EXC_BOTH_nac_e', 'EXC_BOTH_btn_e'};
for i = 1:numel(vitamin)
    new7maxsOUT=setParam(new7maxsOUT,'lb',vitamin(i),[-0.0001]);
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



%% The simulation
%  The growth simulation was carried out by constrainted the model different conditions
%  by seting the uptake rate of each carbon substrate to 1 mmol/gDW/h 
%  (uptake one carbon substrate at a time)


carbon_sources = {'EXC_BOTH_glyc_e',...
    'EXC_BOTH_lcts_e','EXC_BOTH_sucr_e', ...
    'EXC_BOTH_rib__D_e','EXC_BOTH_gal_e', 'EXC_BOTH_malt_e',...
    'EXC_BOTH_arab__L_e', 'EXC_BOTH_tre_e', ...
    'EXC_BOTH_raffin_e', 'EXC_BOTH_melib_e',...
    'EXC_BOTH_fuc__L_e', 'EXC_BOTH_mnl_e',...
    'EXC_BOTH_glc__D_e', 'EXC_BOTH_glcn__D_e', ...
    'EXC_BOTH_fru_e'};
for i = 1:numel(carbon_sources)
    new7maxsOUT=setParam(new7maxsOUT,'lb',carbon_sources{i},[0]);
    new7maxsOUT=setParam(new7maxsOUT,'ub',carbon_sources{i},[0]);
end
cs = cell(numel(carbon_sources),1);
mVersatile_0 = cell(numel(carbon_sources),1);
mVersatile_1 = cell(numel(carbon_sources),1);
mVersatile_2 = cell(numel(carbon_sources),1);
mVersatile_4 = cell(numel(carbon_sources),1);
mVersatile_8 = cell(numel(carbon_sources),1);
mVersatile_16 = cell(numel(carbon_sources),1);
mVersatile_20 = cell(numel(carbon_sources),1);
mVersatile_24 = cell(numel(carbon_sources),1);
model_mVersatile_III = new7maxsOUT; % initial model no carbon substrates
for i = 1:numel(carbon_sources)
    cs{i}=carbon_sources{i};
    model = model_mVersatile_III;
    model_0 = setParam(model,'lb',carbon_sources{i},0);
    model_1 = setParam(model,'lb',carbon_sources{i},-1);
    model_2 = setParam(model,'lb',carbon_sources{i},-2);
    model_4 = setParam(model,'lb',carbon_sources{i},-4);
    model_8 = setParam(model,'lb',carbon_sources{i},-8);
    model_16 = setParam(model,'lb',carbon_sources{i},-16);
    model_20 = setParam(model,'lb',carbon_sources{i},-20);
    model_24 = setParam(model,'lb',carbon_sources{i},-24);
    sol0 = solveLP(model_0,1);
    sol1 = solveLP(model_1,1);
    sol2 = solveLP(model_2,1);
    sol4 = solveLP(model_4,1);
    sol8 = solveLP(model_8,1);
    sol16 = solveLP(model_16,1);
    sol20 = solveLP(model_20,1);
    sol24 = solveLP(model_24,1);
    mVersatile_0{i} = (sol0.f*-1);
    mVersatile_1{i} = (sol1.f*-1);
    mVersatile_2{i} = (sol2.f*-1);
    mVersatile_4{i} = (sol4.f*-1);
    mVersatile_8{i} = (sol8.f*-1);
    mVersatile_16{i} = (sol16.f*-1);
    mVersatile_20{i} = (sol20.f*-1);
    mVersatile_24{i} = (sol24.f*-1);
end

variedResult = table(cs,mVersatile_0,mVersatile_1,mVersatile_2,...
    mVersatile_4,mVersatile_8,mVersatile_16,mVersatile_20,...
    mVersatile_24);

