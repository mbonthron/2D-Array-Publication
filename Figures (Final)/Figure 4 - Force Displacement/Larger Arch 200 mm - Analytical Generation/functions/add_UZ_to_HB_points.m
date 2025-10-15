function [] = add_UZ_to_HB_points(root_name,prob,new_name,data)
%ADD_UZ_TO_HB_POINTS Summary of this function goes here
%   Detailed explanation goes here
bd = coco_bd_read(root_name);
HBlbls = coco_bd_labs(root_name, 'HB');

bcrits = zeros(1,length(HBlbls));

for k = 1:length(HBlbls)
    bcrits(k) = coco_bd_val(bd,HBlbls(k),'b');
end


bcrits = unique_tol(bcrits,1e-3);

prob = coco_add_event(prob,'UZ','b',bcrits);

coco(prob,new_name,[],1,data.parameter_names,data.computational_domain)

end

function v_out = unique_tol(v_in, tol)
    % Remove repeated values in v_in within tolerance tol
    % Keeps the first occurrence of each cluster
    
    % Sort the input
    [v_sorted, idx] = sort(v_in(:));
    
    % Initialize keep mask
    keep = true(size(v_sorted));
    
    % Compare consecutive differences
    for i = 2:length(v_sorted)
        if abs(v_sorted(i) - v_sorted(i-1)) <= tol
            keep(i) = false; % remove repeated within tolerance
        end
    end
    
    % Recover original order
    v_out = v_sorted(keep);
    v_out = v_out(ismember(v_in, v_out, 'legacy')); % keep order from input
end