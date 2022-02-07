function nElements=checkStuff(A,toCheck)
%Now loop through the toCheck matrix, starting with the combination with
%the most models. Only elements that weren't in iteration n are considered
%in iteration n+1.
nElements=zeros(size(toCheck,1),1);
alreadyChecked=[];
for i=1:size(toCheck,1)
    I=find(toCheck(i,:));
    inCommon=setdiff(A{I(1)},alreadyChecked);
    for j=2:numel(I)
        inCommon=intersect(inCommon,A{I(j)});
    end
    alreadyChecked=union(alreadyChecked,inCommon);
    nElements(i)=numel(inCommon);
end
end
