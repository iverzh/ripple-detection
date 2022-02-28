function B = padToLength(A, targ_len, padval)
%padToLength: this function symmetrically pads an array A with padval
%(default 0) OR symmetrically trims array A to reach exactly length
%targ_len. If targ_len-length(A) is odd, the extra pad or trim goes on the
%end of the array. Useful after resampling something to get an exact
%required final length.
%   INPUTS
%   A: n-dimensional array to be padded (but only dimension 1 is padded)
%   targ_len: length in dimension 1 to be reached; if less than length(A),
%   A is trimmed, if more, A is padded
%   padval: value used for padding
%   OUTPUTS
%   B: padded array
%TO DO: add option for padding only at beginning or end

if nargin < 3
    padval = false;
end

if mod(targ_len,1) || targ_len < 1
    error('targ_len must be a positive integer.')
end

A_sz = size(A);
len_diff = targ_len - A_sz(1);
if len_diff == 0
    B = A;
    return
end

if ~mod(len_diff,2)
    start_shift = len_diff/2;
    end_shift = len_diff/2;
else
    end_shift = round(len_diff/2);
    start_shift = len_diff-end_shift;
end
if start_shift > 0
    B = [repmat(padval, [start_shift A_sz(2:end)]);...
        A;...
        repmat(padval, [end_shift A_sz(2:end)])];
elseif numel(A_sz) == 2
    B = A(1-start_shift:end+end_shift,:);
else
    A_msk = [false([-start_shift A_sz(2:end)]);...
        true([targ_len A_sz(2:end)]);...
        false([-end_shift A_sz(2:end)])];
    B = reshape(A(A_msk), [targ_len A_sz(2:end)]);
end

end

