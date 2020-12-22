function new_intervals = get_good_segment_indices(interval, bad_segments)

if isempty(bad_segments)
    new_intervals = interval;
    return;
end

n_bad_segments = size(bad_segments, 1);

new_intervals = zeros(n_bad_segments + 1, 2);

% Good segment before first bad segment?
if bad_segments(1,1) > interval(1)  % An interval before 1st bad segment
        new_intervals(1,:) = [interval(1), bad_segments(1,1)-1];
end

% define good segments between bad segments
for ii = 1:n_bad_segments-1
    new_intervals(ii+1,:) = [bad_segments(ii,2)+1, bad_segments(ii+1,1)-1];         
end

% Good segment after last bad segment?
if bad_segments(end,2) < interval(2)
    new_intervals(end,:) = [bad_segments(end,2)+1, interval(2)];
end
end