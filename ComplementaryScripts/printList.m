
function printList(models,toCheck,nElements)
%To guess how many spaces that are needed to align
firstLen=[];
for i=1:size(toCheck,1)
    label=[];
    I=find(toCheck(i,:));
    for j=1:numel(I)
        label=[label models{I(j)}.id '/'];
    end
    if i==1
        firstLen=numel(label);
    end
    nSpaces=firstLen-numel(label);
    fprintf([label(1:end-1) '  ' repmat(sprintf(' '),1,nSpaces) num2str(nElements(i)) '\n']);
end
end

