f1 = fieldnames(data1);
f2 = fieldnames(data2);

for i = 1:length(f1)
    if isfield(data2,(f1{i}))
        if ~isequal(data1.(f1{i}),data2.(f1{i}))
            f1{i}
            fprintf(" is not equal\n")
        end
    else
        f1{i}
        fprintf(" DNE in data2\n")
    end
end