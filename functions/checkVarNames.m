%% check whether variable names match
ClassicRestrictions = SR;
tmp = (ClassicRestrictions{1,1}{2});
for i = 2:length(ClassicRestrictions)
    aux = unique(ClassicRestrictions{1,i}{2});
    for j = 1:length(aux)
        if sum(strcmp(tmp,aux(j))) == 0
            tmp = [tmp, aux(j)];
        end
    end
end

for i = 1:length(ClassicRestrictions)
    aux{i} = ClassicRestrictions{i}{1};
end
spec.shocks = unique(aux,'stable');


for i = 1:length(tmp)
    check(i) = ismember(tmp(i),varNames);
end
    
if sum(check) ~= length(check)
    error('Wrong variable name in shock specification, check whether names match to varNames');
end
