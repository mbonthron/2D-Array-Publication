function [Ahat,A,b] = get_Ahat_from_coco(run_name,b_val,data)
%%
N = data.N;
N_modes = data.N_modes;

%% Load the data from coco for the mode shapes
bd = coco_bd_read(run_name);
UZ = coco_bd_labs(run_name, 'UZ');

% Constraint length
C = data.constraint_count;


%% Load the values of 'b' and Ahat from COCO
bcrits      = zeros(1,length(UZ));
Ahat        = zeros(2*(N*N_modes-C),length(UZ));
stability   = zeros(1,length(UZ));

for k = 1:length(UZ)
    bcrits(k)    = coco_bd_val(bd,UZ(k),'b');
    Ahat(:,k)    = coco_bd_val(bd,UZ(k),'x');
    maxeig       = max(real(coco_bd_val(bd,UZ(k),'eigs')));
    if maxeig > 0
        stability(k) = 1;   % Unstable Solution
    else
        stability(k) = -1;  % Stable Solution
    end
end

% keep only stable solution
bcrits(stability>0) = nan;

%% Recover the missing modes from the system
A = determine_A_from_Ahat(Ahat',data)';

%% Find the b value that is closest to the b value request
[~,idx] = min(abs(b_val-bcrits));

Ahat = Ahat(:,idx);
A    = A(:,idx);
b    = bcrits(idx);
end

