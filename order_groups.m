% Provides the order in v3_summary to use for copy-pasteable listing
% ie clump each exemplar's C and A slice groups together

% Note that the numbers here are based on the v3_summary struct that is
% created when the program is run. To update the numbers, see
% prepare_slices.m for slice numbers for each pattern

function order = order_groups(type)
    if strcmp(type,'pos')
        order = [1,2,3,4,16,5,6,17,18,19,20,21,22,7,8,23,24,9,25,26,10,11,27,28,32,33,34,35,36,15,37];
    elseif strcmp(type,'neg')
        order = [12,13,14,29,30,31];
    end
end
