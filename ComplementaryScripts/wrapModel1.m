% Written by Nachon Raethong, 23-DEC-2019
cd 'C:\Users\User\Documents\GitHub\reuteri'
% Prepare template model from jcm1112 
jcm1112Cobra = readCbModel('Additional_file_6\Additional_file_6a.xml');
jcm1112Raven = ravenCobraWrapper(jcm1112Cobra);

for i = 1:numel(jcm1112Raven.metComps)
    if jcm1112Raven.metComps(i) == 0;
        jcm1112Raven.metComps(i) = 1;
    else
    end
end


blastAC5JCM1112=getBlast('AC5','AC5.fasta','LbReuteri','JCM1112.fasta');
ac5DraftFromJCM1112 = getModelFromHomology({jcm1112Raven},blastAC5JCM1112,'AC5',{},2,false,10^-5,100);
%get the remaining reactions without gene annotation from JCM1112

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

model_ac5 = goodac5DraftPlusJCM1112;
model_ac5.id = 'AC5';



[~, SheetS]=xlsread('addBiomass.xlsx','Sheet1');
ac5Biomass = struct();
ac5Biomass.rxns = SheetS(2:end,1);
ac5Biomass.rxnNames = SheetS(2:end,2);
ac5Biomass.equations = SheetS(2:end,3);

optimizedBiomass = addRxns(model_ac5,ac5Biomass,2,'c',true,false);
optimizedBiomass = setParam(optimizedBiomass,'lb',ac5Biomass.rxns,0);
optimizedBiomass = setParam(optimizedBiomass,'ub',ac5Biomass.rxns,1000);
optimizedBiomass = setParam(optimizedBiomass,'eq',{'BIOMASS','PROTS_LRE','DNAS_LRE','RNAS_LRE'},0);

jcm1112biomass = {'BIOMASS','PROTS_LRE','DNAS_LRE','RNAS_LRE'};
model_ac5b = removeReactions(optimizedBiomass,jcm1112biomass,true,...
            true,true);

model_ac5b = setParam(model_ac5b,'obj','AC5_BIOMASS',1); 
sol = solveLP(model_ac5b);
printFluxes(model_ac5b, sol.x);
fprintf(['umax = ' num2str(sol.f*-1) ' per hour' '\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%D%O%N%E%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
iTN535 = model_ac5b;
iTN535.id = 'iTN535';
iTN535.description = 'Lactobacillus reuteri KUB-AC5 GEM [23-DEC-2019]';
iTN535.annotation = [];
......................................................................
cd D:/Github/reuteri/ModelFiles/xml
model=importModel('D:/GitHub/reuteri/ModelFiles/xml/iTN535.xml')
modelAC5 = model;
model=importExcelModel('D:/GitHub/reuteri/ModelFiles/xlsx/iTN607.xlsx')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%D%O%N%E%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelAC5.S(121,712)  =  47
modelAC5.S(151,712)  = -47
modelAC5.S(318,712)  = -47
modelAC5.S(320,712)  =  47
modelAC5.S(429,712)  =  47;

model=importExcelModel('D:/GitHub/reuteri/ModelFiles/xlsx/iTN607.xlsx')
glucoseAC5 = setParam(model,'lb',{'EX_glc__D_e'},-6.625); %Uptake of rate of glucose
glucoseAC5 = setParam(glucoseAC5,'obj','AC5_BIOMASS',1); %Maximize AC5 biomass/growth
sol_glucoseAC5 = solveLP(glucoseAC5);
printFluxes(glucoseAC5, sol_glucoseAC5.x,false);fprintf(['umax = ' num2str(sol_glucoseAC5.f*-1) ' per hour in glucose' '\n']);


sucroseAC5 = setParam(model,'lb',{'EX_sucr_e'},-4.292); %Uptake of rate of sucrose
sucroseAC5 = setParam(sucroseAC5,'lb',{'EX_glc__D_e'},0); %Uptake of rate of glucose
sucroseAC5 = setParam(sucroseAC5,'obj','AC5_BIOMASS',1); %Maximize AC5 biomass/growth
sol_sucroseAC5 = solveLP(sucroseAC5);
printFluxes(sucroseAC5, sol_sucroseAC5.x,false);
fprintf(['umax = ' num2str(sol_sucroseAC5.f*-1) ' per hour in sucrose' '\n']);

riboseAC5 = setParam(modelAC5,'lb',{'EX_rib__D_e'},-3.824); %Uptake of rate of ribose
riboseAC5 = setParam(riboseAC5,'lb',{'EX_glc__D_e'},0); %Uptake of rate of glucose
riboseAC5 = setParam(riboseAC5,'lb',{'EX_lcts_e'},0); %Uptake of rate of lactose
riboseAC5 = setParam(riboseAC5,'obj','AC5_BIOMASS',1); %Maximize AC5 biomass/growth
sol_riboseAC5 = solveLP(riboseAC5);print
printFluxes(riboseAC5, sol_riboseAC5.x);
fprintf(['umax = ' num2str(sol_riboseAC5.f*-1) ' per hour in ribose' '\n']);


lactoseAC5 = setParam(modelAC5,'lb',{'EX_lcts_e'},-1.205); %Uptake of rate of lactose
lactoseAC5 = setParam(lactoseAC5,'lb',{'EX_glc__D_e'},0); %Uptake of rate of glucose
lactoseAC5 = setParam(lactoseAC5,'obj','AC5_BIOMASS',1); %Maximize AC5 biomass/growth
sol_lactoseAC5 = solveLP(lactoseAC5);
printFluxes(lactoseAC5, sol_lactoseAC5.x);
fprintf(['umax = ' num2str(sol_lactoseAC5.f*-1) ' per hour in lactose' '\n']);

maltoseAC5 = setParam(modelAC5,'lb',{'EX_malt_e'},-2.544); %Uptake of rate of maltose
maltoseAC5 = setParam(maltoseAC5,'lb',{'EX_glc__D_e'},0); %Uptake of rate of glucose
maltoseAC5 = setParam(maltoseAC5,'obj','AC5_BIOMASS',1); %Maximize AC5 biomass/growth
sol_maltoseAC5 = solveLP(maltoseAC5);
printFluxes(maltoseAC5, sol_maltoseAC5.x);
fprintf(['umax = ' num2str(sol_maltoseAC5.f*-1) ' per hour in maltose' '\n']);

fructoseAC5 = setParam(modelAC5,'lb',{'EX_fru_e'},-3.123); %Uptake of rate of fructose
fructoseAC5 = setParam(fructoseAC5,'lb',{'EX_glc__D_e'},0); %Uptake of rate of glucose
fructoseAC5 = setParam(fructoseAC5,'obj','AC5_BIOMASS',1); %Maximize AC5 biomass/growth
sol_fructoseAC5 = solveLP(fructoseAC5);
printFluxes(fructoseAC5, sol_fructoseAC5.x,false);
fprintf(['umax = ' num2str(sol_fructoseAC5.f*-1) ' per hour in fructose' '\n']);

galactoseAC5 = setParam(modelAC5,'lb',{'EX_gal_e'},-8.398); %Uptake of rate of galactose
galactoseAC5 = setParam(galactoseAC5,'lb',{'EX_glc__D_e'},0); %Uptake of rate of glucose
galactoseAC5 = setParam(galactoseAC5,'obj','AC5_BIOMASS',1); %Maximize AC5 biomass/growth
sol_galactoseAC5 = solveLP(galactoseAC5);
printFluxes(galactoseAC5, sol_galactoseAC5.x);
fprintf(['umax = ' num2str(sol_galactoseAC5.f*-1) ' per hour in galactose' '\n']);

%add rxn %ห้ามเว้นวรรคคำ
model=addReaction(model, 'ชื่อrxnหใม่ที่จะadd','reactionFormula','%ใส่สมการrxn');
model=addReaction(model, 'AC713','reactionFormula','g6p_B_c => f6p_c');

%add met %addเฉพาะที่ไม่มี
model = addmetabolites(model, 'AC714','bfg_c'); 

%add gene %addเฉพาะที่ไม่มี ถ้าซ้ำจะไม่เพิ่มขึ้นให้
model = changeGeneAssociation(model, 'AC713', 'AC5u0009GL000150');

%remove metabolite
model = removeMetabolites(model,{'beta-D-Glucose'}, false);

%remove reaction
model = removeRxns(model,'', false);

%covert to cobra
model = ravenCobraWrapper(modelAC);

%Converting COBRA structure to RAVEN..
model_r2 = ravenCobraWrapper(model_c);

%export
exportForGit(model, 'cobramodel', 'D:\GitHub\reuteri\ModelFiles');

%phpp
model=importExcelModel('D:/GitHub/reuteri/ModelFiles/xlsx/iTN607add_met_name.xlsx')
changeCobraSolver ('gurobi', 'all', 1) %change solver to gurobi
modelphpp = model
ATPphppRates = zeros(0);
for i = 0:40
    for j = 0:40
        modelphpp = changeRxnBounds(modelphpp, 'EX_o2_e', -i, 'b'); %O2
        modelphpp = changeRxnBounds(modelphpp, 'EX_glc__D_e', -j, 'b'); %glu
        modelphpp = changeObjective(modelphpp, 'AC5_BIOMASS'); %objective
        FBAphpp = optimizeCbModel(modelphpp, 'max');
        ATPphppRates(i+1,j+1) = FBAphpp.f;
    end
end
showprogress(0,'generating PhPP');
surfl(ATPphppRates) % 3d plot
xlabel('oxygen uptake (mmol/gcDW/h)')
ylabel('glucose uptake (mmol/gcDW/h)')
zlabel('growth rate (/h)')  

%readcbmodel
model = readcbModel ('%ชื่อโมเดล')%cobramodel.xml

%model feature
printModelStats(model,true,true);

%add reaction ratio
model = addRatioReaction (model, {'EX_glc_D[c]', 'EX_glc_D[e]'}, [1; 2])

% this displays an array with reaction names and flux bounds.
{'AC714', '-1000', '1000'}
    https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorialModelManipulation.html
 
%print all flux
printFluxBounds(model);

model=importExcelModel('D:/GitHub/reuteri/ModelFiles/xlsx/AorKEGGHomology.xlsx')

%arabinose non uptake
model=importExcelModel('D:/GitHub/reuteri/ModelFiles/xlsx/iTN535.xlsx')
glucoseAC5 = setParam(model,'eq',{'EX_arab__L_e'},-6); %Uptake of rate of glucose
glucoseAC5 = setParam(model,'lb',{'EX_glc__D_e'},0); %Uptake of rate of glucose
glucoseAC5 = setParam(glucoseAC5,'obj','AC5_BIOMASS',1); %Maximize AC5 biomass/growth
sol_glucoseAC5 = solveLP(glucoseAC5);
printFluxes(glucoseAC5, sol_glucoseAC5.x,false);fprintf(['umax = ' num2str(sol_glucoseAC5.f*-1) ' per hour in glucose' '\n']);


%iTN535
%add 88 reactions
%add_ec_number
%add_subsystems
%add_biomass
