function chunk_indices = find_chunk_indices(indices, min_chunk_length)
% Find chunks of consecutive indices
%
% Parameters
% ----------
% indices (1D array):
%   The indices to be chunked.
% min_chunk_length (int):
%   Minimum length of a chunk.
%
% Return
% -------
% zero_chunk_indices (N x 2 - array):
%   zero_chunk_indices(i,:) are the start and the end indices of the i-th
%   chunk of zero samples.


for ii = 1:min_chunk_length - 1
    indices = indices(diff(indices) == 1);
end

n_indices = length(indices);

if n_indices == 0
    chunk_indices = reshape([], [], 2); % A 0x2 empty matrix
    return
elseif n_indices == 1
    chunk_indices = [indices(1) indices(1) + min_chunk_length - 1];
    return
end
% Each index is the start of a new chunk
indLim = find(diff(indices) > 1) + 1;

% Only one chunk. FIXME: is this still needed? See line 27
if isempty(indLim)
    chunk_indices = [indices(1) indices(end)];
    return
end

n_chunks = length(indLim) + 1; % >= 2
chunk_indices = zeros(n_chunks,2);

% First chunk
% Consecutive indices in indLim means that each index is the start of a
% chunk of minimum length.
if indLim(1) - 1 == 1   % Chunk of length min_chunk_lenght
    chunk_indices(1,:) = ...
        [indices(1) indices(1)+min_chunk_length-1];
else
    chunk_indices(1,:) = [indices(1) indices(indLim(1)-1)];
end

% n_chunks-2 chunks in the middle
for ii = 1:n_chunks-2
    if indLim(ii+1) - indLim(ii) == 1 % Chunk of length min_chunk_length
        chunk_indices(ii+1,:) = ...
            [indices(indLim(ii)) indices(indLim(ii))+min_chunk_length-1];
    else
        chunk_indices(ii+1,:) = ...
            [indices(indLim(ii)) indices(indLim(ii+1)-1)];
    end
end

% Last chunk
if  indLim(end) == n_indices % Chunk of length min_chunk_length
    chunk_indices(n_chunks,:) = ...
        [indices(end) indices(end)+min_chunk_length-1];
else
    chunk_indices(n_chunks,:) = ...
        [indices(indLim(end)) indices(end)+min_chunk_length-1];
end
end