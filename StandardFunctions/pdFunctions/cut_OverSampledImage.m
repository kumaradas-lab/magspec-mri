function cutImage = cut_OverSampledImage(OSimage, ImageSize)
%% Get center part of input image
%
%     cutImage = cut_OverSampledImage(OSimage, ImageSize)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if isscalar(ImageSize)
  cutImage = zeros(ImageSize, 1);
else
  cutImage = zeros(ImageSize);
end
cutImage = OSimage(...
  (floor(end/2)-floor(size(cutImage,1)/2))+(1:size(cutImage,1)),...
  (floor(end/2)-floor(size(cutImage,2)/2))+(1:size(cutImage,2)),...
  (floor(end/2)-floor(size(cutImage,3)/2))+(1:size(cutImage,3)),...
  (floor(end/2)-floor(size(cutImage,4)/2))+(1:size(cutImage,4)));

end
