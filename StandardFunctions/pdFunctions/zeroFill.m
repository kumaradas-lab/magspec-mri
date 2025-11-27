function filled = zeroFill(toFill, outputSize, winsize)
%% Apply zero-filling and filter to k-space
%
%   filled = zeroFill(toFill, outputSize, winsize)
%
% Zero filling is used to essentially enhance the resolution of an acquired
% image. The complex data allows to locate borders of a measured structure on a
% sub-pixel level. For this, the k-space is expanded to a higher frequency range
% and filled with zeros before calculating the image via FFT.
% Additionally, this function allows to smooth out images by reducing the
% amplitude of high frequency portions (Gibbs ringing, truncation artifact).
%
%
% INPUT:
%
%   toFill
%         complex double multi-dimensional array containing the k-space data
%         (frequency 0 at the center).
%   outputSize
%         vector containing the dimensions of the zero filled output k-space.
%   winsize
%         scalar or vector containing the relative size of the k-space filter or
%         structure with the following fields:
%     size
%           scalar or vector containing the relative size of the k-space filter.
%           The highest k-space frequency in each dimension corresponds to 1.
%           (no default, mandatory)
%     winexp
%           scalar with the exponent used in the filter function (default: 1)
%     used
%           Boolean value indicating whether a filter should be applied
%           (default: true if nargin >=3; false otherwise)
%     winType
%           string indicating the type of the k-space filter function (default:
%           'RaisedCos'). All filter functions are symmetric with respect to the
%           center of the k-space and only depend on the kartesian distance (R)
%           from the k-space center. The multiplicative filter is applied before
%           zero-filling.
%           The following strings are supported:
%       'RaisedCos':
%             Raised cosine function (cos(R*pi).^(2.*winsize.winexp) where the
%             relative kartesian distance to the k-space center is smaller
%             than 0.5 (see winsize.size) and 0 otherwise.
%       'Gaussian':
%             Gaussian bell curve:
%                 exp(- ((R*exp(1/2)*2).^2 .* winsize.winexp)
%
%
% OUTPUT:
%
%   filled
%         The zero-filled and filtered k-space with the same number of
%         dimensions as input toFill and size according to outputSize.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2022 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
if nargin < 3
  winsize.used = 0;
end
if ~isa(winsize, 'struct')
  t = winsize;
  clear winsize
  winsize.size = t;
else
  % backwards compatibility with deprecated interface
  if isfield(winsize, 'i1'), winsize.size(1) = winsize.i1; end
  if isfield(winsize, 'i2'), winsize.size(2) = winsize.i2; end
  if isfield(winsize, 'i3'), winsize.size(3) = winsize.i3; end
  if isfield(winsize, 'i4'), winsize.size(4) = winsize.i4; end
  if isfield(winsize, 'i5'), winsize.size(5) = winsize.i5; end
end

if isfield(winsize, 'size') && isscalar(winsize.size)
  winsize.size(1:5) = winsize.size;
end

if ~isfield(winsize, 'winexp'), winsize.winexp = 1; end
if ~isfield(winsize, 'used'), winsize.used = 1; end
if ~isfield(winsize, 'winType'), winsize.winType = 'RaisedCos'; end


%% apply filter
if winsize.used
  % use input precision
  if isa(toFill, 'single')
    castFcn = @single;
    dataClass = 'single';
  else
    castFcn = @double;
    dataClass = 'double';
  end

  % get "index vectors" with distance to center for each dimension
  r = cell(1, 5);
  sizeIn = size(toFill);
  for t = 1:numel(sizeIn)
    if sizeIn(t) == 1
      r{t} = 0;
    elseif mod(sizeIn(t), 2)
      r{t} = linspace((-0.5+0.5/size(toFill,t)) /winsize.size(t), (0.5-0.5/size(toFill,t)) /winsize.size(t), sizeIn(t)).';
    else
      r{t} = linspace((-0.5)                    /winsize.size(t), (0.5-1.0/size(toFill,t)) /winsize.size(t), sizeIn(t)).';
    end
    r{t} = reshape(castFcn(r{t}), [], 1);  % make sure orientation is in dim 1
  end

  % calculate n-dimensional Eucledian distance from center
  dist = r{1}.^2;
  for iDim = 2:numel(sizeIn)
    if sizeIn(iDim) ~= 0
      permVector = 1:5;
      permVector(1:iDim) = permVector(1:iDim) + 1;
      permVector(iDim) = 1;
      dist = bsxfun(@plus, dist, permute(r{iDim}.^2, permVector));
    end
  end
  dist = sqrt(dist);

  % calculate filter window
  switch winsize.winType
    case 'RaisedCos'
      dist(dist>0.5) = 0.5;
      filter_win = cosd(dist*180) .^ (2 .* winsize.winexp);
    case 'Gaussian'
      filter_win = exp(- ((dist*exp(0.5)*2).^2) .* winsize.winexp);
    otherwise
      error('PD:zeroFill:UnknownFilter', ...
        'Use one of the filter window functions: "RaisedCos" or "Gaussian".');
  end
  clear dist;

  filter_win(filter_win<eps(dataClass)/1e6) = 0;

  % apply filter window
  toFill = toFill .* filter_win;

  clear filter_win;
end


%% pad with zeros
filled = zeros(outputSize, 'like', toFill);
filled(...
  (floor(end/2)-floor(size(toFill,1)/2))+(1:size(toFill,1)),...
  (floor(end/2)-floor(size(toFill,2)/2))+(1:size(toFill,2)),...
  (floor(end/2)-floor(size(toFill,3)/2))+(1:size(toFill,3)),...
  (floor(end/2)-floor(size(toFill,4)/2))+(1:size(toFill,4)),...
  (floor(end/2)-floor(size(toFill,5)/2))+(1:size(toFill,5))) = toFill;


end
