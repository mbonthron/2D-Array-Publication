function copy = deepCopyStruct(original)
    if isstruct(original)
        copy = struct();  % Initialize empty struct
        fields = fieldnames(original);
        for i = 1:numel(fields)
            field = fields{i};
            value = original.(field);
            if isstruct(value)
                copy.(field) = deepCopyStruct(value);  % Recursive copy
            elseif iscell(value)
                copy.(field) = cellfun(@(x) deepCopyStruct(x), value, 'UniformOutput', false);
            else
                copy.(field) = value;  % Direct assignment
            end
        end
    else
        copy = original;  % For non-struct input, return input itself
    end
end