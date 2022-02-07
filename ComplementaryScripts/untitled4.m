initCobraToolbox;
cd '/Users/nachonase/Documents/GitHub/Limosilactobacillus_reuteri_KUB_AC5-GSMM';
model_t=readCbModel('ComplementaryData/LbReuteri_jcm1112.xml');
jcm1112Raven=ravenCobraWrapper(model_t);
% For parsing Cobra model to RAVEN model, we created a new argument of metComps 
% into the template model.
for i = 1:numel(jcm1112Raven.metComps)
    if jcm1112Raven.metComps(i) == 0;
        jcm1112Raven.metComps(i) = 1;
    else
    end
end


model_g=readCbModel('ComplementaryData/iHL622.xml');
model_iHL622=ravenCobraWrapper(model_g);
[~, textData1]=xlsread('ComplementaryData/updated.xlsx','updated_MetNames');
metNames = struct();
metNames.old = textData1(2:end,1);
metNames.new = textData1(2:end,2);
[a, b]=ismember(model_iHL622.metNames,metNames.old);
I=find(a);
model_iHL622.metNames(I)=metNames.new(b(I));


model_i=readCbModel('model/xml/iTN656.xml');
model_iTN656=ravenCobraWrapper(model_i);
models = {model_iTN656, model_iHL622, jcm1112Raven};
printResults=true;
%compStruct=compareModels(models,true);

compStruct.modelIDs={};
for i=1:numel(models)
    compStruct.modelIDs=[models{i}.id];
    models{i}.equ=constructEquations(models{i},models{i}.rxns,true,true,true);
    models{i}.uEqu=constructEquations(models{i},models{i}.rxns,false,true,true);
end

field='rxns';
compStruct.rxns.comparison=getToCheck(models,field);
compStruct.rxns.nElements=checkStuff(getElements(models,field),compStruct.rxns.comparison);
if printResults==true
    fprintf('*** Comparison of reaction IDs:\n');
    printList(models,compStruct.rxns.comparison,compStruct.rxns.nElements);
    fprintf('\n\n');
end

field='mets';
compStruct.mets.comparison=getToCheck(models,field);
compStruct.mets.nElements=checkStuff(getElements(models,field),compStruct.mets.comparison);
if printResults==true
    fprintf('*** Comparison of metabolite IDs:\n');
    printList(models,compStruct.mets.comparison,compStruct.mets.nElements);
    fprintf('\n\n');
end

field='genes';
compStruct.genes.comparison=getToCheck(models,field);
compStruct.genes.nElements=checkStuff(getElements(models,field),compStruct.genes.comparison);
if printResults==true
    fprintf('*** Comparison of gene IDs:\n');
    printList(models,compStruct.genes.comparison,compStruct.genes.nElements);
    fprintf('\n\n');
end

field='eccodes';
compStruct.eccodes.comparison=getToCheck(models,field);
compStruct.eccodes.nElements=checkStuff(getElements(models,field),compStruct.eccodes.comparison);
if printResults==true
    fprintf('*** Comparison of ec-numbers:\n');
    printList(models,compStruct.eccodes.comparison,compStruct.eccodes.nElements);
    fprintf('\n\n');
end

field='metNames';
compStruct.metNames.comparison=getToCheck(models,field);
compStruct.metNames.nElements=checkStuff(getElements(models,field),compStruct.metNames.comparison);
if printResults==true
    fprintf('*** Comparison of metabolite names:\n');
    printList(models,compStruct.metNames.comparison,compStruct.metNames.nElements);
    fprintf('\n\n');
end

field='equ';
compStruct.equ.comparison=getToCheck(models,field);
compStruct.equ.nElements=checkStuff(getElements(models,field),compStruct.equ.comparison);
if printResults==true
    fprintf('*** Comparison of equations with compartment:\n');
    printList(models,compStruct.equ.comparison,compStruct.equ.nElements);
    fprintf('\n\n');
end

field='uEqu';
compStruct.uEqu.comparison=getToCheck(models,field);
compStruct.uEqu.nElements=checkStuff(getElements(models,field),compStruct.uEqu.comparison);
if printResults==true
    fprintf('*** Comparison of equations without compartment:\n');
    printList(models,compStruct.uEqu.comparison,compStruct.uEqu.nElements);
    fprintf('\n\n');
end




printModelStats(model_iTN656,false,false);

printModelStats(jcm1112Raven,false,false);

printModelStats(model_iHL622,false,false);



jcm1112Raven.id = 'jcm1112Raven';
jcm1112Raven.name = 'jcm1112Raven';
