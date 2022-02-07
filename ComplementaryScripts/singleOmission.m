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
    new7maxsOUT=setParam(new7maxsOUT,'lb',amino_acids(i),[-1]);
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

vitaminB = {'EXC_BOTH_pydam_e', 'EXC_BOTH_thm_e',...
    'EXC_BOTH_pnto__R_e', 'EXC_BOTH_nac_e',...
    'EXC_BOTH_btn_e'};

nutrients = union(amino_acids,vitaminB);

for i = 1:numel(nutrients)
    new7maxsOUT3=setParam(new7maxsOUT2,'lb',nutrients(i),[-1]);
    new7maxsOUT3=setParam(new7maxsOUT3,'ub',nutrients(i),[1000]);
end

new7maxsOUT3 = setParam(new7maxsOUT3,'lb',{'EXC_BOTH_glc__D_e'},[0]); 

%% Omission test
omiModel=new7maxsOUT3;
omiModelr=setParam(omiModel,'lb',{'EXC_BOTH_pro__L_e'},[0]);
solomiModel = solveLP(omiModelr,1);

fprintf(['umax omiModel = ' num2str(solomiModel.f*-1) ' per hour' '\n']);

cs=cell(numel(nutrients),1);
mVersatile_1 = cell(numel(nutrients),1);
mVersatile_05 = cell(numel(nutrients),1);
mVersatile_025 = cell(numel(nutrients),1);
mVersatile_0125 = cell(numel(nutrients),1);
mVersatile_00625 = cell(numel(nutrients),1);
mVersatile_00312 = cell(numel(nutrients),1);
mVersatile_00156 = cell(numel(nutrients),1);
mVersatile_00078 = cell(numel(nutrients),1);
mVersatile_00039 = cell(numel(nutrients),1);
mVersatile_00019 = cell(numel(nutrients),1);
mVersatile_00001 = cell(numel(nutrients),1);
mVersatile_0 = cell(numel(nutrients),1);

for i = 1:numel(nutrients)
    model = omiModel;
    cs{i}=nutrients{i};
    model1 = setParam(model,'lb',nutrients{i},[-1]);
    model05 = setParam(model,'lb',nutrients{i},[-0.5]);
    model025 = setParam(model,'lb',nutrients{i},[-0.25]);
    model0125 = setParam(model,'lb',nutrients{i},[-0.125]);
    model00625 = setParam(model,'lb',nutrients{i},[-0.0625]);
    model00312 = setParam(model,'lb',nutrients{i},[-0.0312]);
    model00156 = setParam(model,'lb',nutrients{i},[-0.0156]);
    model00078 = setParam(model,'lb',nutrients{i},[-0.0078]);
    model00039 = setParam(model,'lb',nutrients{i},[-0.0039]);
    model00019 = setParam(model,'lb',nutrients{i},[-0.0019]);
    model00001 = setParam(model,'lb',nutrients{i},[-0.0001]);
    model0 = setParam(model,'lb',nutrients{i},[0]);
    
    sol1 = solveLP(model1,1);
    sol05 = solveLP(model05,1);
    sol025 = solveLP(model025,1);
    sol0125 = solveLP(model0125,1);
    sol00625 = solveLP(model00625,1);
    sol00312 = solveLP(model00312,1);
    sol00156 = solveLP(model00156,1);
    sol00078 = solveLP(model00078,1);
    sol00039 = solveLP(model00039,1);
    sol00019 = solveLP(model00019,1);
    sol00001 = solveLP(model00001,1);
    sol0 = solveLP(model0,1);

    mVersatile_1{i} = (sol1.f*-1);
    mVersatile_05{i} = (sol05.f*-1);
    mVersatile_025{i} = (sol025.f*-1);
    mVersatile_0125{i} = (sol0125.f*-1);
    mVersatile_00625{i} = (sol00625.f*-1);
    mVersatile_00312{i} = (sol00312.f*-1);
    mVersatile_00156{i} = (sol00156.f*-1);
    mVersatile_00078{i} = (sol00078.f*-1);
    mVersatile_00039{i} = (sol00039.f*-1);
    mVersatile_00019{i} = (sol00019.f*-1);
    mVersatile_00001{i} = (sol00001.f*-1);
    mVersatile_0{i} = (sol0.f*-1);
end


variedResult = table(cs,mVersatile_1,mVersatile_05,mVersatile_025,...
    mVersatile_0125,mVersatile_00625,mVersatile_00312,...
    mVersatile_00156,mVersatile_00078,mVersatile_00039,...
    mVersatile_00019,mVersatile_00001,mVersatile_0);
