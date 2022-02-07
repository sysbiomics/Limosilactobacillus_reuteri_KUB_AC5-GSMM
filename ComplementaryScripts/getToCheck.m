function toCheck=getToCheck(models,field)
%Get all the combinations that should be checked for overlap (including the
%single ones)
toCheckA=[];
I=find(cellfun(@checkField,models));
nI=numel(I);
for i=nI:-1:1
    combs=combnk(1:nI,i);
    toAdd=false(size(combs,1),nI);
    for j=1:size(combs,1)
        toAdd(j,combs(j,:))=true;
    end
    toCheckA=[toCheckA;toAdd];
end

%If not all of the models have the required field
toCheck=false(size(toCheckA,1),numel(models));
toCheck(:,I)=toCheckA;

%Ugly thing to get around parameters
    function I=checkField(A)
        I=isfield(A,field);
    end
end
