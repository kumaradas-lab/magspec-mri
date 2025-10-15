function d = hex2uint64(h)
%hex2uint64 Convert hexadecimal string to hex2uint64.

if iscellstr(h), h = char(h); end
if isempty(h), d = []; return, end

% Work in upper case.
h = upper(h);

[m,n]=size(h);

if n>16
   error('PD:hex2uint64:TooLargeHexadecimal', 'Hexadecimal number is too large');
end

% Right justify strings and form 2-D character array.
if ~isempty(find((h==' ' | h==0),1))
  h = strjust(h);

  % Replace any leading blanks and nulls by 0.
  h(cumsum(h ~= ' ' & h ~= 0,2) == 0) = '0';
else
  h = reshape(h,m,n);
end

% Check for out of range values
if any(any(~((h>='0' & h<='9') | (h>='A'&h<='F'))))
   error('PD:hex2uint64:IllegalHexadecimal', 'Illegal Hexadecimal');
end

p = uint64(16).^uint64(n-1:-1:0);
p = p(ones(m,1),:);

d = h <= 64; % Numbers
h(d) = h(d) - 48;

d =  h > 64; % Letters
h(d) = h(d) - 55;

d = sum(uint64(h).*p,2,'native');
