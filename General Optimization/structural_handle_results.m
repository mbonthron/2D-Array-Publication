index_idx       = find(strcmp(results(1,:),'Index'));
binary_idx      = find(strcmp(results(1,:),'Binary'));
connected_idx   = find(strcmp(results(1,:),'Connected'));
stable_idx      = find(strcmp(results(1,:),'Stable Found'));
asym_idx        = find(strcmp(results(1,:),'Asymmetric Wells'));
vhigh_idx       = find(strcmp(results(1,:),'VHigh'));
vlow_idx        = find(strcmp(results(1,:),'VLow'));
thigh_idx       = find(strcmp(results(1,:),'thetaHigh'));
tlow_idx        = find(strcmp(results(1,:),'thetaLow'));
archd_idx       = find(strcmp(results(1,:),'Arch Displacement'));
prohib_idx      = find(strcmp(results(1,:),'Prohibiting Walls'));


[m,n] = size(results);



%%
for count = 0:(m-2)
    % Check to see  if any of the values are empty
    if ~any(cellfun(@isempty, results(count+2,:)))
        if results{count+2,stable_idx} && results{count+2,asym_idx} && ~results{count+2,prohib_idx} && results{count+2,connected_idx}
             disp(count)
        end
    end
end
