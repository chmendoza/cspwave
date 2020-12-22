function varargout = count_unique(x)
% Wrapper of unique() to count unique values
%
% Parameters
% ----------
% x : 1D array
%
% Returns
% -------
% vals : 1D array
%   Array with unique values
% counts: 1D array
%   Array with frequency of each unique value
%

[vals, ia, ~] = unique(sort(x), 'stable');

ia = ia';
counts = diff(ia);

% Add count of last unique value
counts = [counts length(x) + 1 - ia(end)];

if iscolumn(x)
    counts = counts';       
end

switch nargout
    case 1
        varargout{1} = counts;
    case 2
        varargout{1} = counts;
        varargout{2} = vals;
    otherwise
        error('Wrong number of output arguments\n')
end

end