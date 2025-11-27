function data = plot_Image(varargin)
%% Plot image data
%
%   data = plot_Image(data, AQSlice, img_tag)
%   data = plot_Image(hax, ...)
%
% Plot the image (in data.ImageZ) with matching axes decorations.
% This function can be used with the output of the function get_kSpaceAndImage.
% (Currently,) this function only works for 2-D images.
%
%
% INPUT:
%
%   hax
%       Handle to an axes graphics object. If omitted, the current axes (gca)
%       are used.
%
%   data
%       Structure with image data as returned by the function get_kSpaceAndImage
%       (in particular the data in data.ImageZ is displayed).
%
%   AQSlice
%       Structure with the image encoding information.
%
%   img_tag
%       String that is used to identify the image graphics object. This tag is
%       used to optimize repeated plotting to the same axes parents.
%       (Default: Image_2d_image')
%
%
% OUTPUT:
%
%   data
%       Same as input argument "data" but potentially with the additional field
%       "RoI" (if it didn't exist or didn't match the size of data.ImageZ.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% check input and set default values

if nargin > 0 && ~isempty(varargin{1}) && ishghandle(varargin{1}(1), 'axes')
  hax = varargin{1}(1);
  haveAx = 1;
else
  hax = gca();
  haveAx = 0;
end

if nargin < 2 + haveAx
  error('PD:arrow3d:Nargin', ...
    'Function must be called with at least two input arguments.');
end

data = varargin{1+haveAx};
AQSlice = varargin{2+haveAx};
if nargin < 3+haveAx || isempty(varargin{3+haveAx})
  img_tag = 'Image_2d_image';
else
  img_tag = varargin{3+haveAx};
end

szImage = arrayfun(@(x) size(data.ImageZ, x), 1:6);
if sum(szImage(1:4)>1) ~= 2
  error('PD:plot_Image:UnsupportedDimension', ...
    'Only 2-D images are supported.');
end


%% find Region of Interest
if isemptyfield(data, 'RoI') || ...
    (~(isscalar(data.RoI) && data.RoI == 1) && ...
    ~all(size(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice))) == size(data.RoI)))
  ZeroFillFactorSqueezed = data.ZeroFillFactor(size(data.ImageZ(:,:,:,:,AQSlice.iSlice))>1);
  % FIXME: Is the image always along the first dimension?
  data = smoothRoI(data, AQSlice, data.ImageZ(:,:,:,:,AQSlice.iSlice), 2*ZeroFillFactorSqueezed+1);
  data.RoI = double(data.RoI);  % There is no logical NaN.
  data.RoI(~data.RoI) = NaN;
end


%% plot image
hKids = get(hax, 'Children');
hImg = findobj(hKids, 'flat', 'Tag', img_tag);

if sum(AQSlice.nPhase==1)==2
  % read-phase encoded
  t = find(AQSlice.nPhase>1, 1, 'first');
  aspect = [AQSlice.nRead    / AQSlice.nPhase(t), ...
            AQSlice.sizeRead / AQSlice.sizePhase(t), 1];

  xData = data.Ticks(1).ReadZ;
  if ~isinf(AQSlice.sizeRead)
    xData = xData / AQSlice.LengthUnitScale;
  end
  yData = data.Ticks(t).PhaseZ;
  if ~isinf(AQSlice.sizePhase(t))
    yData = yData / AQSlice.LengthUnitScale;
  end
  cData = abs(permute(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice)), [2,1]) / ...
    AQSlice.AmplitudeUnitScale);
  if isempty(hImg)
    imagesc(xData, yData, cData, 'Parent', hax, 'Tag', img_tag);
  else
    set(hImg, 'XData', xData, 'YData', yData, 'CData', cData);
    axis(hax, 'tight');
  end
  set(hax, 'YDir', 'normal');
  % colorbar
  colormap(hax, gray)
  title(hax, ['Image Amplitude in ' AQSlice.AmplitudeUnit])
  if isinf(AQSlice.sizeRead)
    lengthUnit = 'Hz';
  else
    lengthUnit = AQSlice.LengthUnit;
  end
  xlabel(hax, [AQSlice.ReadCartesianAxis{1} 'read in ' lengthUnit]);
  if isinf(AQSlice.sizePhase(t))
    % FIXME: Correct unit?
    lengthUnit = 'px';
  else
    lengthUnit = AQSlice.LengthUnit;
  end
  ylabel(hax, [AQSlice.PhaseCartesianAxis{t} 'phase(' num2str(t) ') in ' lengthUnit]);
else
  % phase-phase encoded
  t = find(AQSlice.nPhase>1, 2, 'first');
  aspect = [AQSlice.nPhase(t(1))    / AQSlice.nPhase(t(2)), ...
            AQSlice.sizePhase(t(1)) / AQSlice.sizePhase(t(2)), 1];

  xData = data.Ticks(t(2)).PhaseZ;
  if ~isinf(AQSlice.sizePhase(t(2)))
    xData = xData / AQSlice.LengthUnitScale;
  end
  yData = data.Ticks(t(1)).PhaseZ;
  if ~isinf(AQSlice.sizePhase(t(1)))
    yData = yData / AQSlice.LengthUnitScale;
  end
  cData = abs(abs(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice)) / ...
    AQSlice.AmplitudeUnitScale));
  if isempty(hImg)
    imagesc(xData, yData, cData, 'Parent', hax, 'Tag', img_tag);
  else
    set(hImg, 'XData', xData, 'YData', yData, 'CData', cData);
    axis(hax, 'tight');
  end
  set(hax, 'YDir', 'normal');
  % colorbar
  colormap(hax, gray)
  title(hax, ['Image Amplitude in ' AQSlice.AmplitudeUnit]);
  if isinf(AQSlice.sizePhase(t(2)))
    % FIXME: Correct unit?
    lengthUnit = 'px';
  else
    lengthUnit = AQSlice.LengthUnit;
  end
  xlabel(hax, [AQSlice.PhaseCartesianAxis{t(2)} 'phase(' num2str(t(2)) ') in ' lengthUnit]);
  if isinf(AQSlice.sizePhase(t(2)))
    % FIXME: Correct unit?
    lengthUnit = 'px';
  else
    lengthUnit = AQSlice.LengthUnit;
  end
  ylabel(hax, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' lengthUnit]);
end
if any([isnan(aspect), isinf(aspect), aspect==0])
  set(hax, 'DataAspectRatioMode', 'auto');
else
  set(hax, 'DataAspectRatio', [1 1 1]);
end

end
