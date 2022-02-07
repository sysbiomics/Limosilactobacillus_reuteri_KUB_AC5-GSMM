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
for i = 1:numel(ac5_newRxns)
    new2=setParam(new2,'ub',ac5_newRxns.rxns(i),[1000]);
end

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
ac5_newRxns2.subSystems = SheetS(2:end,6);


new22 = addRxns(new2,ac5_newRxns2,3,'',true,true);
for i = 1:numel(ac5_newRxns2.rxns)
    new22=setParam(new22,'ub',ac5_newRxns2.rxns(i),[1000]);
end
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
rxnsToRemove = {'MNLabc','FRUabc','BIOMASS','PROTS_LRE',...
    'DNAS_LRE','RNAS_LRE'};
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

for i = 1:numel(new7maxsIN.rxns)
    if new7maxsIN.rev(i) == 1
        new7maxsIN=setParam(new7maxsIN,'lb',new7maxsIN.rxns(i),[-1000]);
    else
        new7maxsIN=setParam(new7maxsIN,'lb',new7maxsIN.rxns(i),[0])
    end
end

new7maxsIN2 = removeReactions(new7maxsIN,...
    {'CYSTL' 'HL_ALATA_L' 'EXC_BOTH_for_e' 'LTA_total'...
    'THRt2' 'EXC_BOTH_ribflv_e'...
    'EXC_BOTH_gua_e' 'EXC_BOTH_ura_e'...
    'EXC_BOTH_ade_e' 'EXC_BOTH_thymd_e' 'EXC_BOTH_ins_e'...
    'HL_PPK2' 'THMDt2r' 'HL_GALt2' 'HL_ALAALAr'...
    'EXC_BOTH_ocdcea_e' 'EXC_BOTH_m_0006' 'EXC_BOTH_m_0002'...
    'EXC_BOTH_m_0003' 'EXC_BOTH_m_0004' 'EXC_BOTH_hista_e'...
    'EXC_BOTH_nh4_e' 'EXC_BOTH_Levan_L'},...
    true,true,true);
        
c={'ATPS3r''HL_ADK2' 'HL_ADK4' 'CBMKr'...
    'HL_ALATA_L' 'DTMPK' 'CYTK1' 'CYTK2' 'GK1'...
    'NDPK9' 'NDPK8' 'NDPK7' 'NDPK6' 'NDPK5' 'NDPK1'...
    'NDPK4' 'NDPK3' 'NDPK2' 'PPK2r' 'UMPK' 'ADK1' 'ASPK'...
    'PFK' 'R08572' 'R02971' 'R04391' 'R01547' 'HL_ADADir' 'HL_GALKr'};
[a, b]=ismember(new7maxsIN2.rxns,c);
I=find(a);
new7maxsIN2.rev(I)=0;
model=simplifyModel(new7maxsIN2,true,true,true,false);
reducedModel=deleteUnusedGenes(model);

for i = 1:numel(reducedModel.rxns)
    if reducedModel.rev(i) == 1
        reducedModel=setParam(reducedModel,'lb',reducedModel.rxns(i),[-1000])
        reducedModel=setParam(reducedModel,'ub',reducedModel.rxns(i),[1000])
    else
        reducedModel=setParam(reducedModel,'lb',reducedModel.rxns(i),[0])
        reducedModel=setParam(reducedModel,'ub',reducedModel.rxns(i),[1000])
    end
end

new30 = setParam(reducedModel,'obj','AC5_GROWTH',1); 
new301 = setParam(new30,'eq',{'ATPM'},[0.37]); 
solx = solveLP(new301);
printFluxes(new301, solx.x);
fprintf(['umax = ' num2str(solx.f*-1) ' per hour' '\n']);

exchangeRxns=getExchangeRxns(new301);
man={'SHCHD2' 'UPP3MT' 'COPREC2MT' 'COPRECT' 'ARBt2' 'R01547'};
[a, b]=ismember(new301.rxns,exchangeRxns);
I=find(a);
new301.subSystems(I)={'exchange'};

[a, b]=ismember(new301.rxns,man);
I=find(a);
new301.subSystems(I)={'unclassified'};

new301.rxnNotes=new301.subSystems;
new301model=new301;
new301model.subSystems=cell(numel(new301model.rxns),1);
for i=1:numel(new301model.rxns)
    new301model.subSystems{i,1}=cellstr(new301model.rxnNotes{i})
end
new301model.rxnNotes=cell(numel(new301model.rxns),1);



new7maxsOUT=new301model;

rxnsToRemove = {'ACTNdiff', 'ALOX','3HPPt','EXC_BOTH_3hpp_e',...
    'ACALDt','AACPS4','AACPS5','AACPS182','NH4DISex','Cuabc','RBLK2',...
    'DIACTt','EXC_BOTH_diact_e','EXC_BOTH_m_0005','EXC_BOTH_fe2_e',...
    'tr_fe2','FORt','PYRt2','SUCCt2r','EXC_BOTH_k_e',...
    'EXC_BOTH_Levan','EXC_BOTH_mg2_e','EXC_BOTH_mn2_e','EXC_BOTH_m_0007',...
    'FA181tr','EXC_BOTH_pyr_e','EXC_BOTH_na1_e','EXC_BOTH_succ_e',...
    'EXC_BOTH_so4_e','EXC_BOTH_actn__R_e','EXC_BOTH_acald_e',...
    'EXC_BOTH_drib_e','EXC_BOTH_dha_e','EXC_BOTH_glyb_e',...
    'HISTAap','EXC_BOTH_ptrc_e'};

remodel_AC5 = removeReactions(new7maxsOUT,rxnsToRemove,true,...
            true,true);
%%
model_final=remodel_AC5;
model_final.id= 'iTN656';
model_final.name = 'Limosilactobacillus reuteri KUB-AC5 GSMM';
exportForGit(model_final,'iTN656','',{'mat', 'txt', 'xlsx', 'xml'});
