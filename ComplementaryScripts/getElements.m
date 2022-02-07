
function A=getElements(models,field)
A={};
for i=1:numel(models)
    if isfield(models{i},field)
        A=[A;{models{i}.(field)}];
    end
end
end