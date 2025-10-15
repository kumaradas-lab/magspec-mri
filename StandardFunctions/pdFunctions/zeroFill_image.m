function [Filled_image, kSpaceZ, kSpace] = zeroFill_image(ImageToFill, outputSize, winsize)
%% Apply zero-filling and k-space filter to image
%
%   [Filled_image, kSpaceZ, kSpace] = zeroFill_image(ImageToFill, outputSize, winsize)
%
% A k-space is created from the input image with a inverse FFT. The function
% "zeroFill" is used to apply zero-filling and filter to that k-space. The
% resulting image is the FFT of the zero-filled and filtered k-space.
% See also the documentation for function "zeroFill".
%
%
% INPUT:
%
%   ImageToFill
%         complex double multi-dimensional array containing the image data
%   outputSize
%         vector containing the dimensions of the zero filled output image.
%   winsize
%         scalar, vector, or structure containing the relative size of the
%         k-space filter. See documentation for function "zeroFill" for further
%         details.
%
%
% OUTPUT:
%
%   Filled_image
%         The image after zero-filling and filtering the k-space. It has the
%         same number of dimensions as the input image ImageToFill and a size
%         according to outputSize.
%   kSpaceZ
%         The zero-filled and filtered k-space with the same number of
%         dimensions as input ImageToFill and size according to outputSize.
%   kSpace
%         The Fourier transform of the input image ImageToFill used as an input
%         for zero-filling and filtering with the function "zeroFill".
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2020 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if nargin < 3
  winsize.used = 0;
end

if ~isa(winsize, 'struct')
  t = winsize;
  clear winsize
  winsize.size = t;

  if nargin < 3
    winsize.used = 0;
  else
    winsize.used = 1;
  end
end

Filled_image = zeros([outputSize, size(ImageToFill, 5)]);
for iImage = 1:size(ImageToFill, 5)
  kSpace = fftshift(ifftn(ifftshift(ImageToFill(:,:,:,:,iImage))));
  kSpaceZ = zeroFill(kSpace, outputSize, winsize);
  Filled_image(:,:,:,:,iImage) = fftshift(fftn(ifftshift(kSpaceZ)));
end

end
