function p = unwrapVector(p,cutoff)
%LocalUnwrap   Unwraps column vector of phase values.

m = length(p);

% Unwrap phase angles.  Algorithm minimizes the incremental phase variation 
% by constraining it to the range [-pi,pi]
dp = diff(p,1,1);                % Incremental phase variations
dps = mod(dp+pi,2*pi) - pi;      % Equivalent phase variations in [-pi,pi)
dps(dps==-pi & dp>0,:) = pi;     % Preserve variation sign for pi vs. -pi
dp_corr = dps - dp;              % Incremental phase corrections
dp_corr(abs(dp)<cutoff,:) = 0;   % Ignore correction when incr. variation is < CUTOFF

% Integrate corrections and add to P to produce smoothed phase values
p(2:m,:) = p(2:m,:) + cumsum(dp_corr,1);