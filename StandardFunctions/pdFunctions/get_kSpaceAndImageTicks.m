function [data, AQSlice] = get_kSpaceAndImageTicks(data, AQSlice)
%% Calculate k-ticks and ticks for k-space and image
%
%       [data, AQSlice] = get_kSpaceAndImageTicks(data, AQSlice)
%
%
% INPUT:
%
%   data        structure with image and k-space data (as returned from
%               get_kSpaceAndImage)
%     ZeroFillFactor
%                 The factors that will be used for zero-filling the k-space
%                 [read phase1 phase2 phase3].
%     ZeroFillFactorOS
%                 The factor that will be used for zero-filling the over-sampled
%                 k-space [read phase1 phase2 phase3].
%
%   AQSlice     structure with settings
%     nRead       Number of pixels/voxels in read direction.
%     ReadOS      Over-sampling factor in read direction.
%     nPhase      Number of pixels/voxels in phase (1,2,3) direction.
%     PhaseOS     Over-sampling factor in phase (1,2,3) direction.
%     CenterRot   Center of rotation in meter [read phase1 phase2 phase3]
%                 (default: [0 0 0]).
%
%
% OUTPUT:
%
%   data        Structure with image and k-space data (see corresponding input
%               variable) with added fields "Ticks" and "kTicks".
%               Both of this fields are three element vectors of structures. The
%               first element contains the ticks (or k-ticks) in read and
%               phase(1) direction. The second and third elements correspond to
%               phase(2) and phase(3) directions, accordingly. The ticks and
%               k-ticks are stored in the following fields:
%       Read
%                   Vector with the (k-)ticks in read direction corresponding to
%                   Image or kSpace.
%       ReadZ
%                   Vector with the (k-)ticks in read direction corresponding to
%                   the zero-padded kSpaceZ or the corresponding ImageZ.
%       ReadOs
%                   Vector with the (k-)ticks in the over-sampled read direction
%                   corresponding to ImageOs or kSpaceOs.
%       ReadOsZ
%                   Vector with the (k-)ticks in the over-sampled read direction
%                   corresponding to the zero-padded kSpaceOsZ or the
%                   corresponding ImageOsZ.
%       Phase
%                   Vector with the (k-)ticks in phase(1, 2, or 3) direction
%                   corresponding to Image or kSpace.
%       PhaseZ
%                   Vector with the (k-)ticks in phase(1, 2, or 3) direction
%                   corresponding to the zero-padded kSpaceZ or the
%                   corresponding ImageZ.
%       PhaseOs
%                   Vector with the (k-)ticks in the over-sampled phase(1, 2, or
%                   3) direction corresponding to ImageOs or kSpaceOs.
%       PhaseOsZ
%                   Vector with the (k-)ticks in the over-sampled phase(1, 2, or
%                   3) direction corresponding to the zero-padded kSpaceOsZ or
%                   the corresponding ImageOsZ.
%
%               Additionally, the following field is added to data:
%     PermuteOrder
%                 For 3d images, this is the permute order for the image and
%                 k-space arrays to correctly display the images in sliceomatic.
%
%   AQSlice
%               Same as the input AQSlice with the additional fields:
%     SliceCartesianAxis
%                 If the image was rotated by a multiple of 90 degrees around
%                 all rotation axes, this is the label for the slice direction
%                 (in gradient coordinate sytem). It is empty, otherwise.
%     ReadCartesianAxis
%                 If the image was rotated by a multiple of 90 degrees around
%                 all rotation axes, this is the label for the read direction
%                 (in gradient coordinate sytem). It is empty, otherwise.
%     PhaseCartesianAxis
%                 If the image was rotated by a multiple of 90 degrees around
%                 all rotation axes, this is 1x3 cell string with the labels for
%                 the phase directions (in gradient coordinate sytem). It is a
%                 1x3 cell string with empty strings, otherwise.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2021 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------
%
% See also: get_kSpaceAndImage, plot_kSpaceAndImage

%% check input
if numel(AQSlice) > 1
  warning('PD:plot_kSpaceAndImage:MultipleAQSlice', ...
    'Multiple AQSlices are not supported. First slice is used.');
  AQSlice = AQSlice(1);
end

%% default input
if isemptyfield(AQSlice, 'CenterRot'), AQSlice.CenterRot = [0, 0, 0]; end
if isemptyfield(AQSlice, 'angle2Turns'), AQSlice.angle2Turns = 1/(2*pi); end
if isemptyfield(data, 'ZeroFillFactor')
  if isemptyfield(AQSlice, 'ZeroFillFactor')
    error('PD:get_kSpaceAndImageTicks:ZeroFillFactorUndefined', ...
      '"ZeroFillFactor" must be set in "data" or "AQSlice".');
  else
    MySize = [AQSlice.nRead, AQSlice.nPhase(1), AQSlice.nPhase(2), AQSlice.nPhase(3)];
    MySizeZ = ceil(MySize .* AQSlice.ZeroFillFactor);
    data.ZeroFillFactor = MySizeZ ./ MySize;
  end
end
if isemptyfield(data, 'ZeroFillFactorOS')
  MySizeOS = [AQSlice.nRead*AQSlice.ReadOS, AQSlice.nPhase(1)*AQSlice.PhaseOS(1), ...
    AQSlice.nPhase(2)*AQSlice.PhaseOS(2), AQSlice.nPhase(3)*AQSlice.PhaseOS(3)];
  MySizeOSZ = ceil(MySizeOS .* data.ZeroFillFactor);
  data.ZeroFillFactorOS = MySizeOSZ ./ MySizeOS;
end

%% find CartesianAxis names
Angle2Deg = AQSlice.angle2Turns * 360;
[Rx, Ry, Rz] = get_aptDegRotationMatrix(AQSlice.alfa(1)*Angle2Deg, ...
  AQSlice.phi(1)*Angle2Deg, AQSlice.theta(1)*Angle2Deg);
R = Rz * Ry * Rx;
Coordinates = {'-Z ', '-Y ', '-X ', '', 'X ', 'Y ', 'Z '};
SliceCoordinate = zeros(3, 1);
SliceCoordinate(AQSlice.SliceCoordinate) = 1;
SliceCoordinate = R * SliceCoordinate;
SliceCoordinate = [find(SliceCoordinate>0.9), -find(-SliceCoordinate>0.9)];
if isempty(SliceCoordinate), SliceCoordinate = 0; end
AQSlice.SliceCartesianAxis = Coordinates(SliceCoordinate+4);

ReadCoordinate = zeros(3,1);
ReadCoordinate(AQSlice.ReadCoordinate) = 1;
ReadCoordinate = R * ReadCoordinate;
ReadCoordinate = [find(ReadCoordinate>0.9), -find(-ReadCoordinate>0.9)];
if isempty(ReadCoordinate), ReadCoordinate = 0; end
AQSlice.ReadCartesianAxis = Coordinates(ReadCoordinate+4);

for iPhase = 1:3
  temp = zeros(3,1);
  temp(AQSlice.PhaseCoordinate(iPhase)) = 1;
  temp = R * temp;
  temp = [find(temp>0.9), -find(-temp>0.9)];
  if isempty(temp), temp = 0; end
  AQSlice.PhaseCartesianAxis(iPhase) = Coordinates(temp+4);
end


%% grids
if isinf(AQSlice.sizeRead)
  % k-Space label of axis x,y and z tick
  % Spatial frequency read for k-Space axis
  data.kTicks(1).Read = get_FFTGrid(AQSlice.nRead, AQSlice.nRead).';
  data.kTicks(1).ReadZ = ...
    get_FFTGrid(AQSlice.nRead*data.ZeroFillFactor(1), AQSlice.nRead*data.ZeroFillFactor(1)).';
  % Spatial frequency read for kSpaceOs axis
  data.kTicks(1).ReadOs = get_FFTGrid(AQSlice.nRead*AQSlice.ReadOS, AQSlice.nRead*AQSlice.ReadOS).';
  data.kTicks(1).ReadOsZ = ...
    get_FFTGrid(AQSlice.nRead*data.ZeroFillFactorOS(1)*AQSlice.ReadOS, ...
    AQSlice.nRead*data.ZeroFillFactorOS(1)*AQSlice.ReadOS).';
  % data.kTicks(1).ReadHzZ = get_FFTGrid(AQSlice.nRead*ReadZ/AQSlice.sizeRead, AQSlice.nRead*ReadZ).';
  % Ticks of Image
  % data.Ticks(1).Read = AQSlice.CenterRot(1) + get_FFTGrid(AQSlice.nRead, AQSlice.nRead).';
  % data.Ticks(1).ReadZ = AQSlice.CenterRot(1) + ...
  %   get_FFTGrid(AQSlice.nRead*data.ZeroFillFactor(1), AQSlice.nRead*data.ZeroFillFactor(1)).';
  % data.Ticks(1).ReadOs = AQSlice.CenterRot(1) + ...
  %   get_FFTGrid(AQSlice.nRead*AQSlice.ReadOS, AQSlice.nRead*AQSlice.ReadOS).';
  % data.Ticks(1).ReadOsZ = AQSlice.CenterRot(1) + ...
  %   get_FFTGrid(AQSlice.nRead*AQSlice.ReadOS*data.ZeroFillFactorOS(1), ...
  %               AQSlice.nRead*AQSlice.ReadOS*data.ZeroFillFactorOS(1)).';
  data.Ticks(1).Read = get_FFTGrid(AQSlice.nRead/AQSlice.AcquisitionTime, AQSlice.nRead).';
  data.Ticks(1).ReadZ = get_FFTGrid(AQSlice.nRead/AQSlice.AcquisitionTime, ...
                                    AQSlice.nRead*data.ZeroFillFactor(1)).';
  data.Ticks(1).ReadOs = get_FFTGrid(AQSlice.nRead*AQSlice.ReadOS/AQSlice.AcquisitionTime, ...
                                     AQSlice.nRead*AQSlice.ReadOS).';
  data.Ticks(1).ReadOsZ = get_FFTGrid(AQSlice.nRead*AQSlice.ReadOS/AQSlice.AcquisitionTime, ...
                                      AQSlice.nRead*AQSlice.ReadOS*data.ZeroFillFactorOS(1)).';
else
  % k-Space label of axis x,y and z tick
  % Spatial frequency read for k-Space axis
  data.kTicks(1).Read = get_FFTGrid(AQSlice.nRead/AQSlice.sizeRead, AQSlice.nRead).';
  data.kTicks(1).ReadZ = ...
    get_FFTGrid(AQSlice.nRead*data.ZeroFillFactor(1)/AQSlice.sizeRead, ...
    AQSlice.nRead*data.ZeroFillFactor(1)).';
  % Spatial frequency read for kSpaceOs axis
  data.kTicks(1).ReadOs = get_FFTGrid(AQSlice.nRead/AQSlice.sizeRead, ...
    AQSlice.nRead*AQSlice.ReadOS).';
  data.kTicks(1).ReadOsZ = ...
    get_FFTGrid(AQSlice.nRead*data.ZeroFillFactorOS(1)/AQSlice.sizeRead, ...
    AQSlice.nRead*data.ZeroFillFactorOS(1)*AQSlice.ReadOS).';
  % data.kTicks(1).ReadHzZ=get_FFTGrid(AQSlice.nRead*ReadZ/AQSlice.sizeRead,AQSlice.nRead*ReadZ).';
  % Ticks of Image
  data.Ticks(1).Read = AQSlice.CenterRot(1) + ...
    get_FFTGrid(AQSlice.sizeRead, AQSlice.nRead).';
  data.Ticks(1).ReadZ = AQSlice.CenterRot(1) + ...
    get_FFTGrid(AQSlice.sizeRead, AQSlice.nRead*data.ZeroFillFactor(1)).';
  data.Ticks(1).ReadOs = AQSlice.CenterRot(1) + ...
    get_FFTGrid(AQSlice.sizeRead*AQSlice.ReadOS, AQSlice.nRead*AQSlice.ReadOS).';
  data.Ticks(1).ReadOsZ = AQSlice.CenterRot(1) + ...
    get_FFTGrid(AQSlice.sizeRead*AQSlice.ReadOS, AQSlice.nRead*AQSlice.ReadOS*data.ZeroFillFactor(1)).';
end
for iPhase = 1:numel(AQSlice.nPhase)
  if isinf(AQSlice.sizePhase(iPhase))
    % Spatial frequency phase for k-Space axis
    data.kTicks(iPhase).Phase = ...
      get_FFTGrid(AQSlice.nPhase(iPhase), AQSlice.nPhase(iPhase)).';
    data.kTicks(iPhase).PhaseZ = ...
      get_FFTGrid(AQSlice.nPhase(iPhase).*data.ZeroFillFactor(iPhase+1), ...
      AQSlice.nPhase(iPhase).*data.ZeroFillFactor(iPhase+1)).';
    % Spatial frequency phase for kSpaceOs axis
    data.kTicks(iPhase).PhaseOs = ...
      get_FFTGrid(AQSlice.nPhase(iPhase).*AQSlice.PhaseOS(iPhase), ...
      AQSlice.nPhase(iPhase).*AQSlice.PhaseOS(iPhase)).';
    data.kTicks(iPhase).PhaseOsZ = ...
      get_FFTGrid(AQSlice.nPhase(iPhase).*data.ZeroFillFactorOS(iPhase+1).*AQSlice.PhaseOS(iPhase), ...
      AQSlice.nPhase(iPhase).*data.ZeroFillFactorOS(iPhase+1).*AQSlice.PhaseOS(iPhase)).';
    data.Ticks(iPhase).Phase = AQSlice.CenterRot(iPhase) + ...
      get_FFTGrid(AQSlice.nPhase(iPhase), AQSlice.nPhase(iPhase)).';
    data.Ticks(iPhase).PhaseZ = AQSlice.CenterRot(iPhase) + ...
      get_FFTGrid(AQSlice.nPhase(iPhase).*data.ZeroFillFactor(iPhase+1), ...
      AQSlice.nPhase(iPhase).*data.ZeroFillFactor(iPhase+1)).';
    data.Ticks(iPhase).PhaseOs = AQSlice.CenterRot(iPhase) + ...
      get_FFTGrid(AQSlice.nPhase(iPhase)*AQSlice.PhaseOS(iPhase), ...
      AQSlice.nPhase(iPhase)*AQSlice.PhaseOS(iPhase)).';
    data.Ticks(iPhase).PhaseOsZ = AQSlice.CenterRot(iPhase) + ...
      get_FFTGrid(AQSlice.nPhase(iPhase)*AQSlice.PhaseOS(iPhase).*data.ZeroFillFactor(iPhase+1), ...
      AQSlice.nPhase(iPhase)*AQSlice.PhaseOS(iPhase).*data.ZeroFillFactor(iPhase+1)).';
  else
    % Spatial frequency phase for k-Space axis
    data.kTicks(iPhase).Phase = ...
      get_FFTGrid(AQSlice.nPhase(iPhase)./AQSlice.sizePhase(iPhase), ...
      AQSlice.nPhase(iPhase)).';
    data.kTicks(iPhase).PhaseZ = ...
      get_FFTGrid(AQSlice.nPhase(iPhase).*data.ZeroFillFactor(iPhase+1)./AQSlice.sizePhase(iPhase), ...
      AQSlice.nPhase(iPhase).*data.ZeroFillFactor(iPhase+1)).';
    % Spatial frequency phase for kSpaceOs axis
    data.kTicks(iPhase).PhaseOs = ...
      get_FFTGrid(AQSlice.nPhase(iPhase)./AQSlice.sizePhase(iPhase), ...
      AQSlice.nPhase(iPhase).*AQSlice.PhaseOS(iPhase)).';
    data.kTicks(iPhase).PhaseOsZ = ...
      get_FFTGrid(AQSlice.nPhase(iPhase).*data.ZeroFillFactorOS(iPhase+1)./AQSlice.sizePhase(iPhase), ...
      AQSlice.nPhase(iPhase).*data.ZeroFillFactorOS(iPhase+1).*AQSlice.PhaseOS(iPhase)).';
    data.Ticks(iPhase).Phase = AQSlice.CenterRot(iPhase) + ...
      get_FFTGrid(AQSlice.sizePhase(iPhase), AQSlice.nPhase(iPhase)).';
    data.Ticks(iPhase).PhaseZ = AQSlice.CenterRot(iPhase) + ...
      get_FFTGrid(AQSlice.sizePhase(iPhase), ...
      AQSlice.nPhase(iPhase).*data.ZeroFillFactor(iPhase+1)).';
    data.Ticks(iPhase).PhaseOs = AQSlice.CenterRot(iPhase) + ...
      get_FFTGrid(AQSlice.sizePhase(iPhase)*AQSlice.PhaseOS(iPhase), ...
      AQSlice.nPhase(iPhase)*AQSlice.PhaseOS(iPhase)).';
    data.Ticks(iPhase).PhaseOsZ = AQSlice.CenterRot(iPhase) + ...
      get_FFTGrid(AQSlice.sizePhase(iPhase)*AQSlice.PhaseOS(iPhase), ...
      AQSlice.nPhase(iPhase)*AQSlice.PhaseOS(iPhase).*data.ZeroFillFactor(iPhase+1)).';
  end
end

%% permute Order (for sliceomatic)
% default PermuteOrder for non-3D images (no permutation)
data.PermuteOrder = [1, 2, 3, 4];

if sum([AQSlice.nPhase(:)>1; AQSlice.nRead>1]) >= 3
  % 3d image
  if AQSlice.nRead > 1
    % read encoding for first dimension
    data.PermuteOrder = [2, 1 ,3];
  else
    % phase encoding for all 3 dimensions (CSI)
    data.PermuteOrder = [1, 3, 2];
  end
end

end
