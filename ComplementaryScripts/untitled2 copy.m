model_final=importModel('model/xml/iTN656.xml',false);
model = model_final;
carbonSourceRxns = {'EXC_BOTH_melib_e', 'EXC_BOTH_raffin_e',...
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

aminoSupplyRxns = {'EXC_BOTH_arg__L_e', 'EXC_BOTH_cys__L_e',...
    'EXC_BOTH_ile__L_e', 'EXC_BOTH_leu__L_e', 'EXC_BOTH_lys__L_e',...
    'EXC_BOTH_met__L_e', 'EXC_BOTH_thr__L_e', 'EXC_BOTH_tyr__L_e',...
    'EXC_BOTH_val__L_e', 'EXC_BOTH_ala__L_e', 'EXC_BOTH_asn__L_e',...
    'EXC_BOTH_phe__L_e', 'EXC_BOTH_trp__L_e', 'EXC_BOTH_pro__L_e',...
    'EXC_BOTH_gln__L_e', 'EXC_BOTH_asp__L_e', 'EXC_BOTH_gly_e',...
    'EXC_BOTH_ser__L_e', 'EXC_BOTH_glu__L_e', 'EXC_BOTH_his__L_e'};

carbonSources = struct();
carbonSources.rxns = carbonSourceRxns;
carbonSources.rxnNames = cell(numel(carbonSourceRxns),1);
[a, b]=ismember(carbonSources.rxns, model.rxns);
I=find(a);
carbonSources.rxnNames(I)=model.rxnNames(b(I));



aminoSupply = struct();
aminoSupply.rxns = aminoSupplyRxns;
aminoSupply.rxnNames = cell(numel(aminoSupplyRxns),1);
[a, b]=ismember(aminoSupply.rxns, model.rxns);
I=find(a);
aminoSupply.rxnNames(I)=model.rxnNames(b(I));

amino_acids=aminoSupplyRxns;
carbon_sources=carbonSourceRxns;
new7maxsOUT=model;

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

products = {'EXC_BOTH_hxan_e','EXC_BOTH_xan_e',...
    'EXC_BOTH_succ_e','EXC_BOTH_gcald_e', 'EXC_BOTH_btd__RR_e',...
    'EXC_BOTH_Levan', 'EXC_BOTH_ac_e', 'EXC_BOTH_lac__D_e',...
    'EXC_BOTH_orot_e','EXC_BOTH_etoh_e', 'EXC_BOTH_diact_e',...
    'EXC_BOTH_mal__L_e', 'EXC_BOTH_lac__L_e', 'EXC_BOTH_13ppd_e'...
    'EXC_BOTH_acald_e', 'EXC_BOTH_orn__L_e', 'EXC_BOTH_actn__R_e',...
    'EXC_BOTH_drib_e', 'EXC_BOTH_mal__L_e', 'EXC_BOTH_pyr_e',...
    'EXC_BOTH_succ_e'};
for i = 1:numel(products)
    new7maxsOUT=setParam(new7maxsOUT,'lb',products(i),[0]);
    new7maxsOUT=setParam(new7maxsOUT,'ub',products(i),[1000]);
end


for i = 1:numel(carbon_sources)
    new7maxsOUT=setParam(new7maxsOUT,'lb',carbon_sources{i},[0]);
    new7maxsOUT=setParam(new7maxsOUT,'ub',carbon_sources{i},[0]);
end

model_mVersatile_III = new7maxsOUT; % initial model no carbon substrates
for i = 1:numel(carbon_sources)
    model = model_mVersatile_III;
    model_i = setParam(model,'lb',carbon_sources{i},-1);
    for j = 1:numel(amino_acids)
        model_j = setParam(model_i,'lb',amino_acids{j},0);
        sol = solveLP(model_j,1);
        fprintf(['umax  = ' num2str(sol.f*-1) ' per hour' '\n'])
    end
end


variedResult = table(cs,mVersatile_III);



% create a structure and cells for the results
growthProfile = struct();
growthProfile.carbonSources = cell(numel(variedN),1);
growthProfile.aminoSupplyRxns = cell(numel(variedN),1);
growthProfile.umax = cell(numel(variedN),1);


