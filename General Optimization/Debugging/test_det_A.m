clearvars;
load("./det_A_test.mat");

A_1 = A;
tic
for i = 1:C % Incredibly slow
    mode = modes_to_skip(i);
    A_1 = [A_1(1:mode-1,:); missingvals(i,:) ; A_1(mode:end,:)];
end
shift_modes = N*N_modes;    % Do the same for the derivative terms
for i = 1:C
    mode = modes_to_skip(i);
    A_1 = [A_1(1:shift_modes+mode-1,:); Dmissingvals(i,:) ; A_1(shift_modes+mode:end,:)];
end
toc

tic
% Preallocate
shift_modes = N*N_modes;
[M_size, N_size] = size(A);
C = length(modes_to_skip);
A_new = zeros(M_size + 2*C, N_size);

insert_locs = false(M_size + 2*C, 1);
insert_locs(modes_to_skip) = true;
insert_locsD = false(M_size + 2*C, 1);
insert_locsD(modes_to_skip+shift_modes) = true;

% Create final matrix
A_new(insert_locs, :) = missingvals;
A_new(insert_locsD, :) = Dmissingvals;
A_new(~(insert_locs| insert_locsD), :) = A;

% Derivatives
% [M_size, N_size] = size(A_new);
% C = length(modes_to_skip);
% A_new_new = zeros(M_size + C, N_size);
% 
% insert_locs = false(M_size + C, 1);
% insert_locs(modes_to_skip+shift_modes) = true;
% 
% % Create final matrix
% A_new_new(insert_locs, :) = Dmissingvals;
% A_new_new(~insert_locs, :) = A_new;
toc

isequal(A_1(1:shift_modes),A_new(1:shift_modes))
isequal(A_1,A_new)