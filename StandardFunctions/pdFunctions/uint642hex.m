function h = uint642hex(d,n)
%uint642hex Convert decimal integer (uint64) to hexadecimal string.
%   Example
%       uint642hex(2748,16) returns '0000000000000ABC'.
%       uint642hex(2748,4)  returns '0ABC'.
%       uint642hex(2748)    returns 'ABC'.

if nargin<1
    error(nargchk(1,2,nargin,'struct'));
end
if nargin<2
    n=16;
end

d = d(:); % Make sure d is a column vector.

if any(~isreal(d)) || any(d < 0) || any(d ~= fix(d))
    error('PD:uint642hex:FirstArgIsInvalid','FirstArgIsInvalid')
end
if any(d > uint64(2)^64)
    warning('PD:uint642hex:FirstArgTooLarge','FirstArgTooLarge');
end

d = uint64(d);

numD = numel(d);

nMin=0;
if nargin==1 || any(d > uint64(2)^(n*4)),
    for t=0:4:63
        if any(d > uint64(2)^t);
            nMin=t/4+1;
        else
            break;
        end
    end
end

if nargin==2 && nMin>n
    warning('PD:uint642hex:SecondArgIsTooSmall','SecondArgIsTooSmall');
    n=nMin;
elseif nargin==1
    n=nMin;
end

if numD>1
    n = n*ones(numD,1);
end

h = sprintf('%0*X',[n,d]');

h = reshape(h,n(1),numD)';
