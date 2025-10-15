function A = convnfft(A, B, shape, dims, Power2Flag)
% CONVNFFT  FFT-BASED N-dimensional convolution.
%   C = CONVNFFT(A, B) performs the N-dimensional convolution of
%   matrices A and B. If nak = size(A,k) and nbk = size(B,k), then
%   size(C,k) = max([nak+nbk-1,nak,nbk]);
% 
%   C = CONVNFFT(A, B, SHAPE) controls the size of the answer C:
%       'full'   - (default) returns the full N-D convolution
%       'same'   - returns the central part of the convolution that
%                  is the same size as A.
%       'valid'  - returns only the part of the result that can be
%                  computed without assuming zero-padded arrays.
%                  size(C,k) = max([nak-max(0,nbk-1)],0).
%
%   C = CONVNFFT(..., SHAPE, DIMS) with DIMS is vector of dimensions where
%       the convolution will be carried out. By default DIMS is
%       [1:max(ndims(A),ndims(B))] (all dimensions). A and B must have the
%       same lengths on other dimensions.
% 
%   'Power2Flag', boolean. If it is TRUE, use FFT with length rounded
%       to the next power-two. It is faster but requires more memory.
%       Default value is FALSE.
%
% Class support for inputs A,B:
% float: double, single
%
% METHOD: CONVNFFT uses Fourier transform (FT) convolution theorem, i.e.
%         FT of the convolution is equal to the product of the FTs of the
%         input functions.
%         In 1-D, the complexity is O((na+nb)*log(na+nb)), where na/nb are
%         respectively the lengths of A and B.
%
% Usage recommendation:
%         In 1D, this function is faster than CONV for nA, nB > 1000.
%         In 2D, this function is faster than CONV2 for nA, nB > 20.
%         In 3D, this function is faster than CONVN for nA, nB > 5.
% 
% See also conv, conv2, convn.
% 
% derived from: convnfft (Bruno Luong <brunoluong@yahoo.com>)

if nargin<3 || isempty(shape)
    shape = 'full';
end
nd = max(ndims(A),ndims(B));
% work on all dimensions by default
if nargin<4 || isempty(dims)
    dims = 1:nd;
end
dims = reshape(dims, 1, []); % row (needed for for-loop index)
if nargin<5 || isempty(Power2Flag)
    Power2Flag = 0;
end
% IFUN function will be used later to truncate the result
% M and N are respectively the length of A and B in some dimension
switch lower(shape)
    case 'full',
        ifun = @(m,n) 1:m+n-1;
    case 'same',
        ifun = @(m,n) ceil((n-1)/2)+(1:m);
    case 'valid',
        ifun = @(m,n) n:m;
    otherwise
        error('convnfft: unknown shape %s', shape);
end
classA = class(A);
classB = class(B);
ABreal = isreal(A) && isreal(B);
% Special case, empty convolution, try to follow MATLAB CONVN convention
if any(size(A)==0) || any(size(B)==0)
    szA = zeros(1,nd); szA(1:ndims(A))=size(A);
    szB = zeros(1,nd); szB(1:ndims(B))=size(B);
    % Matlab wants these:
    szA = max(szA,1); szB = max(szB,1);
    szC = szA;
    for dim=dims
        szC(dim) = length(ifun(szA(dim),szB(dim)));
    end
    A = zeros(szC,classA); % empty -> return zeros
    return
end
if Power2Flag
    % faster FFT if the dimension is power of 2
    lfftfun = @(l) 2^nextpow2(l);
else
    % slower, but smaller temporary arrays
    lfftfun = @(l) l;
end

% Matlab FFT
subs(1:ndims(A)) = {':'};
for dim=1:nd
    m = size(A,dim);
    n = size(B,dim);
    % compute the FFT length
    l(dim) = lfftfun(m+n-1);
    subs{dim} = ifun(m,n);
end
if nd==numel(unique(dims));
    A = fftn(A,l);
    B = fftn(B,l);
else
  for dim=dims
    A = fft(A,l(dim),dim);
    B = fft(B,l(dim),dim);
  end
end

A = A.*B;
clear B

% Back to the non-Fourier space
% Matlab IFFT  
if nd==numel(unique(dims));
    A = ifftn(A);
else
  for dim=dims
    A = ifft(A,[],dim);
  end
end
% Truncate the results
if ABreal
    % Make sure the result is real
    A = real(A(subs{:}));
else
    A = A(subs{:});
end
end % convnfft