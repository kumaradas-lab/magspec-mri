function data = smoothRoI(data, AQSlice, image, kernel_sz)
%% Calculate Boolean array "RoI", remove isolated outliers and fill isolated holes
%
%   data = smoothRoI(data, AQSlice, image, kernel_sz)
%
% Calculate a Boolean array with the same size as "image" which is true in the
% "region of interest" (RoI).
%
%
% INPUT:
%
%   data
%       data structure as returned, e.g., by get_kSpaceAndImage.
%
%   AQSlice
%       Structure with the following (optional) fields:
%     RoiAbsoluteValue
%         If set, this value is used as the absolute cut-off value for the RoI.
%     RoiCutOffPercentile, RoiRelativeValue
%         If RoiAbsoluteValue is not set, these values are used to calculate the
%         absolute cut of value from the absolute value of the data in "image".
%
%   image
%       ND-array with the image data.
%
%   kernel_sz
%       Size of the convolution kernel that is used to remove isolated outliers
%       and fill isolated holes.
%
%
% OUTPUT:
%
%   data
%       Same as the input structure data. But with the following additional
%       fields:
%     RoICutOff
%         The selected cut-off value (see the documentation for the fields in
%         AQSlice above).
%     RoI
%         Logical ND-array with the size of non-singleton dimensions of the
%         input argument "image" where the region of interest is marked with
%         true.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% prepare input
if isemptyfield(AQSlice, 'RoiAbsoluteValue')
  % determine value of lower percentile
  if AQSlice.RoiCutOffPercentile == 0 || AQSlice.RoiRelativeValue == 0
    data.RoICutOff = 0;
  else
    imageSorted = sort(abs(image(:)));
    data.RoICutOff = imageSorted(round(numel(image)*AQSlice.RoiCutOffPercentile)) ...
      * AQSlice.RoiRelativeValue;
  end
else
  data.RoICutOff = AQSlice.RoiAbsoluteValue;
end

if data.RoICutOff == 0
  data.RoI = 1;
  return;
end

% permute singleton dimensions of image to end
img_sz = size(image);
image = permute(image, [find(img_sz~=1), find(img_sz==1)]);
img_sz = size(image);
if numel(img_sz) == 2
  img_sz=[img_sz,1];
end
% make sure that the convolution kernel is 3d
if isempty(kernel_sz)
  kernel_sz = [1, 1, 1];
elseif isscalar(kernel_sz)
  kernel_sz = [kernel_sz, 1, 1];
elseif numel(kernel_sz) == 2
  kernel_sz = [kernel_sz, 1];
elseif numel(kernel_sz) > 3
  kernel_sz = kernel_sz(1:3);
end

zffp = (kernel_sz-1)/2;  % ZeroFillFactorPermuted

RoI = double(abs(image) >= data.RoICutOff);


%% pad with the values mirrored at the border
% create RoI with space for mirrored values
data.RoI = NaN(size(RoI,1)+zffp(1)*2, size(RoI,2)+zffp(2)*2, size(RoI,3)+zffp(3)*2);
% mirror at borders in first dim
data.RoI(:,zffp(2)+(1:size(RoI,2)),zffp(3)+(1:size(RoI,3))) = ...
  RoI([1+(zffp(1):-1:1),1:size(RoI,1),size(RoI,1)-(1:zffp(1))],:,:);
if numel(kernel_sz) > 1 && zffp(2) >= 1
  % mirror at borders in second dim
  data.RoI(:,:,zffp(3)+(1:size(RoI,3))) = ...
    data.RoI(:,zffp(2)+[1+(zffp(2):-1:1),1:size(RoI,2),size(RoI,2)-(1:zffp(2))],zffp(3)+(1:size(RoI,3)));
end
if numel(kernel_sz) > 2 && zffp(3) >= 1
  % mirror at borders in third dim
  data.RoI(:,:,:) = ...
    data.RoI(:,:,zffp(3)+[1+(zffp(3):-1:1),1:size(RoI,3),size(RoI,3)-(1:zffp(3))]);
end


%% smooth the image and get RoI
% convolution with kernel of weight 1
data.RoI = fftn(data.RoI);
kernel_roi = zeros(size(data.RoI));  % pad with zeros
kernel_idx = arrayfun(@(x) 1:x, kernel_sz, 'UniformOutput', false);
kernel_roi(kernel_idx{:}) = 1/prod(kernel_sz);
kernel_roi = circshift(kernel_roi, -zffp);  % shift center of kernel to idx 1
kernel_roi = fftn(kernel_roi);
data.RoI = data.RoI .* kernel_roi;
data.RoI = ifftn(data.RoI);

% select center of RoI (size of image)
img_idx = arrayfun(@(x,y) x+(1:y), zffp, img_sz, 'UniformOutput', false);
data.RoI = data.RoI(img_idx{1:ndims(data.RoI)});

% finally create RoI array
data.RoI = data.RoI > 1/2;

end
