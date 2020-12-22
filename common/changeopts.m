function opts = changeopts(opts, key, value)
% Change key-value pairs on a cell array
%
% Examples:
% opts = changeopts(opts, 'sym', true)
% opts = changeopts(opts, {'max_iter','sym'},{500,true})

if ~iscell(key)
    key = {key};
    value = {value};
end

ind = cellfun(@(x)find(strcmpi(opts,x)),key);
opts(ind+1) = value;

end