function ind = interval_union(ind)
% Finds the union of intervals
% Parameters
% ----------
% ind (float):
%   A matrix of size [n 2], with n being the number of intervals. ind(:,1)
%   is the start point of the interval and ind(:,2) the end point. It might
%   have intersecting intervals.
%
% Returns
% -------
% ind (float)
%   A matrix of size [m 2], with m <= n. It contains non-intersecting
%   intervals.

if isempty(ind)
    return;
end

% Union of intersecting intervals
%% Keep the longest range in repeating indices
ind = sortrows(ind);
[~, ia, ~] = unique(ind(:,1), 'last');
ind = ind(ia,:);
[~, ia, ~] = unique(ind(:,2), 'first');
ind = ind(ia,:);
ind = sortrows(ind);

%% Delete intermediate intersecting ranges
% Keep row with largest end point: keep largest interval and delete smaller
% subintervals that are inside. A n B = A.
I = [true; diff(ind(:,2))>0];
ind = ind(I,:);

% Unite intersecting intervals. A n B ~= A.
% Find sets of rows (intervals) that are intersecting. Pick the start point of 
% the first row of each set as the start point of the intersection of the set:
I = [true; diff([ind(1:end-1,2) ind(2:end,1)],1,2)>0];
c1 = ind(I, 1); % Start indices of new intervals
% identify sets composed of only one interval. their end point doesn't
% change:
one_ind = find(I);
if one_ind
    if one_ind(end) == length(I)
        J = [diff(one_ind)==1; true];
    else
        J = [diff(one_ind)==1; false];
    end
    one_ind = one_ind(J); % rows that do not intersect with any other row
end
% for sets with more than one interval, pick the end point of the last row
% as the end point of the intersection of the set:
intersect_ind = find(~I); % rows that intersect with other rows
if intersect_ind
    K = [diff(intersect_ind)>1; true]; % pick the last row of each set
    intersect_ind = intersect_ind(K);
end
c2 = sort([one_ind; intersect_ind]);
c2 = ind(c2,2);
ind = [c1 c2];