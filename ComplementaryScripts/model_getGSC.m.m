%% get GSC from iTN656 
% Purpose: to extract gene-set from iTN656
% 
% Written by Nachon Raethong, 10-DEC-2021
%% WORKSPACE
cd '/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM';
model=importModel('model/xml/iTN656.xml',false);
model=readCbModel('model/xlsx/iTN656.xlsx');
model=importExcelModel('model/xlsx/iTN656.xlsx',false);
extractSubsystemGSC(model,'','GSC/iTN656_subSystems.txt')
extractMetaboliteGSC(model,'','GSC/iTN656_mets.txt')
submodel=model;
submodel.subSystems=submodel.rxnNames;
extractSubsystemGSC(submodel,'','GSC/iTN656_rxnNames.txt')
submodel.equations = constructEquations(submodel);
submodel.subSystems=submodel.equations;
extractSubsystemGSC(submodel,'','GSC/iTN656_equations.txt')