function [a1stab,a1unst,a3stab,a3unst,b] = get_data_from_coco(run_name)
%GET_DATA_FROM_COCO Summary of this function goes here
%   Detailed explanation goes here
bd = coco_bd_read(run_name);

a_data = bd(2:end,19);
b_data = bd(2:end,8);
eigs   = bd(2:end,15);

b = cellfun(@(x) x(1),b_data);

a1 = cellfun(@(x) x(1),a_data);
a2 = cellfun(@(x) x(2),a_data);
a3 = cellfun(@(x) x(3),a_data);


lambda = cellfun(@(x) max(real(x)),eigs);

a1 = a1 - a3;

[a1stab,a1unst] = separate(a1,lambda);
[a2stab,a2unst] = separate(a2,lambda);
[a3stab,a3unst] = separate(a3,lambda);

end

