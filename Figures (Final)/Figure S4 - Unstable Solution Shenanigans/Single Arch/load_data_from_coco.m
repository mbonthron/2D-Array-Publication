function [b,xstab,xunst] = load_data_from_coco(run_name)
bd = coco_bd_read(run_name);

b = cellfun(@(x) x(1),bd(2:end,8));

x = cellfun(@(x) x(1),bd(2:end,16));

eigs = bd(2:end,15);
lambda1 = cellfun(@(x) max(real(x)),eigs);

[xstab,xunst] = separate(x,lambda1);

end

