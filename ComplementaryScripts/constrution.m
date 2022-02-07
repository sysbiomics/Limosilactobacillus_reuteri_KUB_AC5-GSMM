%% Metabolic network construction and development of L. reuteri KUB-AC5 genome-scale metabolic model (GSMM)
% Purposed: Here is a complete script used for constructing a metabolic network and
% developing a GSMM of L. reuteri KUB-AC5.
% The process was carried through the following steps:
% STEP 0: TEMPLATE MODEL PREPARATION
% STEP 1: STEP 2: SEQUENCE ANALYSIS FOR CONSTRUCTING A METABOLIC NETWORK
% STEP 2: TEMPLATE MODEL PREPARATION
% STEP 3: TEMPLATE MODEL PREPARATION

% This is an iterative process which ends up with a validated GSMM 
% that can represent L. reuteri KUB-AC5 growth on different substrates.
% Due to the close proximity between L. reuteri KUB-AC5 and strain JCM1112 genomes 


% Written by Nachon Raethong, 23-DEC-2019
% Updated by Nachon Raethong, 21-NOV-2021

%% WORKSPACE
cd '/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM';
%% STEP 1: TEMPLATE MODEL PREPARATION
% Initially, the GSMM of Limosilactobacillus reuteri JCM1112 was prepared and used as 
% an intial template network.
initCobraToolbox;

jcm1112Cobra = readCbModel('ComplementaryData/LbReuteri_jcm1112.xml');
jcm1112Raven = ravenCobraWrapper(jcm1112Cobra);
% For parsing Cobra model to RAVEN model, we created a new argument of metComps 
% into the template model.
for i = 1:numel(jcm1112Raven.metComps)
    if jcm1112Raven.metComps(i) == 0;
        jcm1112Raven.metComps(i) = 1;
    else
    end
end

%% STEP 2: SEQUENCE ANALYSIS FOR CONSTRUCTING A METABOLIC NETWORK
% The genome sequences of L. reuteri KUB-AC5 and JCM1112 were retrieved
% from Jatuponwiphat, T. et al., 2019 and Kristjansdottir, T. et al., 2019, respectively. 
% Then, the pairwise orthologous genes were performed between L. reuteri KUB-AC5 and 
% L. reuteri JCM 1112 using the best bidirectional best hits (BBHs) from BLASTP.

%blastAC5JCM1112=getBlast('AC5','ComplementaryData/AC5.fasta','LbReuteri','ComplementaryData/JCM1112.fasta');
%save('ComplementaryData/blastAC5JCM1112.mat','blastAC5JCM1112');

load('ComplementaryData/blastAC5JCM1112.mat');

% The identified orthologous genes were implemented to generate an initial metabolic 
% network of L. reuteri KUB-AC5 metabolism based on a set of existing genes, 
% reactions and metabolites in the template.

ac5DraftFromJCM1112 = getModelFromHomology({jcm1112Raven},blastAC5JCM1112,'AC5',{},2,false,10^-5,100);

% Non-gene-associated reactions, including spontaneous, transport and exchange reactions 
% and literature-based metabolic activities, were manually curated and
% added for optimizing global connectivity and feasibility of the metabolic network.

rxnToAdd.rxns =setdiff(jcm1112Raven.rxns,ac5DraftFromJCM1112.rxns);

ac5DraftPlusJCM1112 = addRxnsGenesMets(ac5DraftFromJCM1112,jcm1112Raven,rxnToAdd.rxns,...
    false,'additional rxns from JCM1112 model',2);


badrxnsToRemove = {'12PPDt','PPOHt'};
goodac5DraftPlusJCM1112 = removeReactions(ac5DraftPlusJCM1112,badrxnsToRemove,true,...
            true,true);
goodac5DraftPlusJCM1112 = setParam(goodac5DraftPlusJCM1112,'obj',{'BIOMASS'},1);
sol = solveLP(goodac5DraftPlusJCM1112,1);
printFluxes(goodac5DraftPlusJCM1112, sol.x,true);
fprintf(['umax = ' num2str(sol.f*-1) ' per hour' '\n']);

model_ac5=contractModel(goodac5DraftPlusJCM1112);

printModelStats(model_ac5,true,true);
%exportForGit(model_ac5,'model_ac5','');


%% STEP 3: GAP-FILLING AND POLISHING GSMM 
% 3.1 update grRules
model_ac5plus = model_ac5;
[~, newGrRules]=xlsread('ComplementaryData/updated.xlsx','updated_grRules');
rxnID = newGrRules(2:end,1);
geneAssoc = newGrRules(2:end,2);
for i = 1:numel(rxnID)
    rxn = rxnID(i);
    new = geneAssoc{i};
    model_ac5plus = changeGeneAssoc(model_ac5plus,rxn,new,true);
end

% 3.1 update metabolite names
model_ac5plus2=model_ac5plus;
[~, textData1]=xlsread('ComplementaryData/updated.xlsx','updated_MetNames');
metNames = struct();
metNames.old = textData1(2:end,1);
metNames.new = textData1(2:end,2);
[a, b]=ismember(model_ac5plus2.metNames,metNames.old);
I=find(a);
model_ac5plus2.metNames(I)=metNames.new(b(I));
printModelStats(model_ac5plus2,true,true);

% 3.1 update EC numbers
model_ac5plus3=model_ac5plus2;
[~, textData1]=xlsread('ComplementaryData/updated.xlsx','updated_RxnECs');
updated_RxnECs = struct();
updated_RxnECs.rxns = textData1(2:end,1);
updated_RxnECs.new = textData1(2:end,2);
[a, b]=ismember(model_ac5plus3.rxns,updated_RxnECs.rxns);
I=find(a);
model_ac5plus3.eccodes(I)=updated_RxnECs.new(b(I));
model_ac5plus3=contractModel(model_ac5plus3);

printModelStats(model_ac5plus3,true,true);

% 3.1 update updated_RxnNames
model_ac5plus4=model_ac5plus3;
[~, textData1]=xlsread('ComplementaryData/updated.xlsx','updated_RxnNames');
updated_RxnNames = struct();
updated_RxnNames.old = textData1(2:end,1);
updated_RxnNames.new = textData1(2:end,2);
[a, b]=ismember(model_ac5plus4.rxnNames,updated_RxnNames.old);
I=find(a);
model_ac5plus4.rxnNames(I)=updated_RxnNames.new(b(I));
model_ac5plus4=contractModel(model_ac5plus4);
printModelStats(model_ac5plus4,true,true);




sol = solveLP(model_ac5plus4,1);
printFluxes(model_ac5plus4, sol.x,true);
fprintf(['umax = ' num2str(sol.f*-1) ' per hour' '\n']);


%% STEP 4: ADD NEWLY-IDENTIFIED REACTIONS
% 4.1 Addition of new metabolites according to their presence in 
% newly-identified reactions to the network.

model_ac5new=model_ac5plus4;
[~, newMets]=xlsread('ComplementaryData/updated.xlsx','new_Metabolites');
metsToAdd = struct();
metsToAdd.mets = newMets(2:end,1);
metsToAdd.metNames = newMets(2:end,2);
metsToAdd.metFormulas = newMets(2:end,3);
metsToAdd.compartments = 'c';
new1=addMets(model_ac5new,metsToAdd);

% 4.2 Introduction of newly-identified reactions to the network.
[~, SheetS]=xlsread('ComplementaryData/updated.xlsx','new_Reactions');
ac5_newRxns = struct();
ac5_newRxns.rxns = SheetS(2:end,1);
ac5_newRxns.rxnNames = SheetS(2:end,2);
ac5_newRxns.equations = SheetS(2:end,3);
ac5_newRxns.eccodes = SheetS(2:end,4);
ac5_newRxns.grRules = SheetS(2:end,5);

new2 = addRxns(new1,ac5_newRxns,2,'c',true,true);
simplifyModel(new2);

% 4.2.HL Addition of iHL reactions to the network.

%blastAC5biogaia=getBlast('AC5','ComplementaryData/AC5.fasta','biogaia','ComplementaryData/Lreuteri_biogaia_v03.faa');
%save('ComplementaryData/blastAC5biogaia.mat','blastAC5biogaia');
%load('ComplementaryData/blastAC5biogaia.mat');

%iHL622c = readCbModel('ComplementaryData/iHL622.xml');
%iHL622r = ravenCobraWrapper(iHL622c);
%iHL622r.id = 'biogaia';
%[~, textData1]=xlsread('ComplementaryData/updated.xlsx','updated_MetNames');
%metNames = struct();
%metNames.old = textData1(2:end,1);
%metNames.new = textData1(2:end,2);
%[a, b]=ismember(iHL622r.metNames,metNames.old);
%I=find(a);
%iHL622r.metNames(I)=metNames.new(b(I));

%ac5DraftFrombiogaia = getModelFromHomology({iHL622r},blastAC5biogaia,'AC5',{},2,false,10^-5,100);
%ac5DraftFrombiogaia.equations = constructEquations(ac5DraftFrombiogaia);
%exportForGit(ac5DraftFrombiogaia,'ac5DraftFrombiogaia','','xlsx');
%rxns=setdiff(ac5DraftFrombiogaia.rxns,new7maxsOUT22.rxns);

[~, SheetS]=xlsread('ComplementaryData/updated.xlsx','new_Reactions2');
ac5_newRxns2 = struct();
ac5_newRxns2.rxns = SheetS(2:end,1);
ac5_newRxns2.rxnNames = SheetS(2:end,2);
ac5_newRxns2.equations = SheetS(2:end,3);
ac5_newRxns2.eccodes = SheetS(2:end,4);
ac5_newRxns2.grRules = SheetS(2:end,5);

new22 = addRxns(new2,ac5_newRxns2,3,'',true,true);
simplifyModel(new22);
% 4.3 Addition of Levan biosysthesis

[~, newMets_Levan]=xlsread('ComplementaryData/updated.xlsx','new_Metabolites_Levan');
metsToAdd_Levan = struct();
metsToAdd_Levan.mets = newMets_Levan(2:end,1);
metsToAdd_Levan.metNames = newMets_Levan(2:end,2);
metsToAdd_Levan.metFormulas = newMets_Levan(2:end,3);
metsToAdd_Levan.compartments = 'e';
new3=addMets(new22,metsToAdd_Levan);

[~, SheetS]=xlsread('ComplementaryData/updated.xlsx','new_Reactions_Levan');
LevanRxns = struct();
LevanRxns.rxns = SheetS(2:end,1);
LevanRxns.rxnNames = SheetS(2:end,2);
LevanRxns.equations = SheetS(2:end,3);
LevanRxns.eccodes = SheetS(2:end,4);
LevanRxns.grRules = SheetS(2:end,5);

new4 = addRxns(new3,LevanRxns,2,'e',true,true);
simplifyModel(new4);
new5=contractModel(new4);


%% STEP 5: ADD BIOMASS REACTIONS OF KUB-AC5
[~, SheetS]=xlsread('ComplementaryData/updated.xlsx','add_BiomassAC5');
add_BiomassAC5 = struct();
add_BiomassAC5.rxns = SheetS(2:end,1);
add_BiomassAC5.rxnNames = SheetS(2:end,2);
add_BiomassAC5.equations = SheetS(2:end,3);
new6 = addRxns(new5,add_BiomassAC5,2,'c',true,false);
for i = 1:numel(add_BiomassAC5.rxns)
    new6=setParam(new6,'lb',add_BiomassAC5.rxns(i),[0]);
    new6=setParam(new6,'ub',add_BiomassAC5.rxns(i),[1000]);
end
rxnsToRemove = {'BIOMASS','PROTS_LRE','DNAS_LRE','RNAS_LRE'};
new7max = removeReactions(new6,rxnsToRemove,true,...
            true,true);
        
new7maxa = setParam(new7max,'obj','AC5_GROWTH',1);

sol = solveLP(new7maxa);
printFluxes(new7maxa, sol.x);
fprintf(['umax = ' num2str(sol.f*-1) ' per hour' '\n']);
%% STEP 6: SUBSYSTEMS CHARACTERIZATION
%model=importExcelModel('ComplementaryData/iTN607.xlsx');
%characterizedSubsystems = struct();
%characterizedSubsystems.rxns = model.rxns;
%characterizedSubsystems.subSystems =model.subSystems;
%save('ComplementaryData/characterizedSubsystems.mat','characterizedSubsystems');

load('ComplementaryData/characterizedSubsystems.mat');
[a, b]=ismember(new7maxa.rxns,characterizedSubsystems.rxns);
I=find(a);
new7maxa.subSystems(I)=characterizedSubsystems.subSystems(b(I));
printModelStats(new7maxa,true,true);

%exportForGit(ac5DraftFrombiogaias,'ac5DraftFrombiogaias','','xlsx');


exchangeRxns = setdiff(getExchangeRxns(new7maxa,'both'),'AC5_GROWTH');
rxnsToRemove = exchangeRxns;
new7maxs = removeReactions(new7maxa,rxnsToRemove,true,...
            true,true);

metNames=getMetsInComp(new7maxs,'e');
[new7maxsIN, addedRxns] = addExchangeRxns(new7maxs,'both',new7maxs.mets(metNames));
new7maxsIN.equations = constructEquations(new7maxsIN);
for i = 1:numel(addedRxns)
    new7maxsIN=setParam(new7maxsIN,'lb',addedRxns(i),[-1000]);
    new7maxsIN=setParam(new7maxsIN,'ub',addedRxns(i),[1000]);
end

new7maxsIN2 = removeReactions(new7maxsIN,...
    {'EXC_BOTH_m_0002' 'EXC_BOTH_m_0003' 'EXC_BOTH_succ_e' 'EXC_BOTH_m_0004' 'EXC_BOTH_hista_e' 'EXC_BOTH_nh4_e' 'EXC_BOTH_Levan_L'},true,...
            true,true);
        
new30 = setParam(new7maxsIN2,'obj','AC5_GROWTH',1); 
new301 = setParam(new30,'eq',{'ATPM'},[0.37]); 
solx = solveLP(new301);
printFluxes(new301, solx.x);
fprintf(['umax = ' num2str(solx.f*-1) ' per hour' '\n']);


new7maxsOUT=new301;
products = {'EXC_BOTH_gcald_e', 'EXC_BOTH_btd__RR_e',...
    'EXC_BOTH_Levan', 'EXC_BOTH_ac_e', 'EXC_BOTH_lac__D_e',...
    'EXC_BOTH_orot_e','EXC_BOTH_etoh_e', 'EXC_BOTH_diact_e',...
    'EXC_BOTH_mal__L_e', 'EXC_BOTH_lac__L_e'};
for i = 1:numel(products)
    new7maxsOUT=setParam(new7maxsOUT,'lb',products(i),[0]);
    new7maxsOUT=setParam(new7maxsOUT,'ub',products(i),[1000]);
end
amino_acids = {'EXC_BOTH_arg__L_e', 'EXC_BOTH_cys__L_e', 'EXC_BOTH_ile__L_e', 'EXC_BOTH_leu__L_e', 'EXC_BOTH_lys__L_e', 'EXC_BOTH_met__L_e', 'EXC_BOTH_thr__L_e', 'EXC_BOTH_tyr__L_e', 'EXC_BOTH_val__L_e', 'EXC_BOTH_ala__L_e', 'EXC_BOTH_asn__L_e', 'EXC_BOTH_phe__L_e', 'EXC_BOTH_trp__L_e', 'EXC_BOTH_pro__L_e', 'EXC_BOTH_gln__L_e', 'EXC_BOTH_asp__L_e', 'EXC_BOTH_gly_e', 'EXC_BOTH_ser__L_e', 'EXC_BOTH_glu__L_e', 'EXC_BOTH_his__L_e'};
for i = 1:numel(amino_acids)
    new7maxsOUT=setParam(new7maxsOUT,'lb',amino_acids(i),[-20]);
end
lipid= {'EXC_BOTH_hdcea_e', 'EXC_BOTH_ocdcea_e', 'EXC_BOTH_ocdcya_e', 'EXC_BOTH_ocdctr_e'};
for i = 1:numel(lipid)
    new7maxsOUT=setParam(new7maxsOUT,'lb',lipid(i),[-0.0001]);
end

vitamin = {'EXC_BOTH_m_0006', 'EXC_BOTH_thm_e', 'EXC_BOTH_pydam_e', 'EXC_BOTH_pnto__R_e', 'EXC_BOTH_nac_e', 'EXC_BOTH_btn_e'};
for i = 1:numel(vitamin)
    new7maxsOUT=setParam(new7maxsOUT,'lb',vitamin(i),[-20]);
end
nucleotide= {'EXC_BOTH_gua_e', 'EXC_BOTH_ura_e', 'EXC_BOTH_ade_e', 'EXC_BOTH_thymd_e', 'EXC_BOTH_ins_e'};
for i = 1:numel(nucleotide)
    new7maxsOUT=setParam(new7maxsOUT,'lb',nucleotide(i),[0]);
end
carbon_sources = {'EXC_BOTH_succ_e', 'EXC_BOTH_melib_e', 'EXC_BOTH_raffin_e', 'EXC_BOTH_pyr_e', 'EXC_BOTH_drib_e', 'EXC_BOTH_glc__D_e', 'EXC_BOTH_glcn__D_e', 'EXC_BOTH_glyc_e', 'EXC_BOTH_lcts_e', 'EXC_BOTH_orn__L_e', 'EXC_BOTH_sucr_e', 'EXC_BOTH_tre_e', 'EXC_BOTH_rib__D_e', 'EXC_BOTH_gal_e', 'EXC_BOTH_malt_e', 'EXC_BOTH_13ppd_e', 'EXC_BOTH_acald_e', 'EXC_BOTH_actn__R_e', 'EXC_BOTH_arab__L_e'};
for i = 1:numel(carbon_sources)
    new7maxsOUT=setParam(new7maxsOUT,'lb',carbon_sources(i),[0]);
end

new7maxsOUT2 = setParam(new7maxsOUT,'lb',{'EXC_BOTH_o2_e'},[-0.000000000001]); 
new7maxsOUT2 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_glc__D_e'},[-6.6]); 
%new7maxsOUT2 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_sucr_e'},[-5.2]); 
%new7maxsOUT2 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_lcts_e'},[-0.9]); 
%new7maxsOUT2 = setParam(new7maxsOUT2,'lb',{'EXC_BOTH_malt_e'},[-3.3]); 

sol = solveLP(new7maxsOUT2);
printFluxes(new7maxsOUT2, sol.x);
fprintf(['umax = ' num2str(sol.f*-1) ' per hour' '\n']);



%EXC_BOTH_glc__D_e 6.634	0.684 > umax 0.151 ± 0.004 [umax = 0.15128]
%EXC_BOTH_lcts_e 0.940	0.322 > umax 0.078 ± 0.005 [umax = 0.10966]
%EXC_BOTH_sucr_e 5.192	0.204 > umax 0.247 ± 0.003 [umax = 0.23871]
%EXC_BOTH_malt_e 3.310	0.764 > umax 0.199 ± 0.009 [umax = 0.1793]
