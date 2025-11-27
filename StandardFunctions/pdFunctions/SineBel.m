function window = SineBel(windowSize)
%%
if nargin < 1
  windowSize = [];
end
if isempty(windowSize)
  windowSize = 2^10;
end

window = sin(linspace(0, pi, windowSize)).';

end
