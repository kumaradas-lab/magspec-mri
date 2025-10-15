function [data, AQSlice] = correct_read_B0(data, AQSlice, dataB0, GradAmpRead)
%% Correct shift in read direction with B0 map
%
%   [data, AQSlice] = correct_read_B0(data, AQSlice, dataB0, GradAmpRead)
%
% The function deskews the image in data.ImageZ to correct the read out scaling
% error due to the deviations in the B0 field (dataB0).
%
% INPUT:
%
%   data
%         Structure with measurement data (as returned in the field SeqLoop.data
%         by e.g. sequence_Flash or sequence_SpinEcho).
%
%   AQSlice
%         AQSlice settings of the measurement to correct. Fields that are
%         specific to this function include:
%     CorrectAmplitude
%           Correct the image amplitude according to how the pixels/voxels are
%           stretched or shrunk in read direction. (Default: true)
%
%   dataB0
%         Structure with the measured B0 map (see Seq.CorrectB0Read.Get in e.g.
%         sequence_Flash or sequence_SpinEcho).
%
%   GradAmpRead
%         Gradient amplitude for readout encoding in T/m.
%
%
% OUTPUT:
%
%   data
%         Similar to input data with the corrected image. The original image
%         (before the correction) is returned in data.ImageZNoB0Corr. The
%         corrected image is returned in data.ImageZ.
%
%   AQSlice
%         Same as input AQSlice but with potentially added fields.
%
%
% SEE ALSO: sequence_Flash, sequence_SpinEcho
%
% ------------------------------------------------------------------------------
% (C) Copyright 2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

if isemptyfield(AQSlice, 'CorrectAmplitude')
  AQSlice.CorrectAmplitude = true;
end

if isempty(dataB0)
  warning('No data for correction of read offset with B0 map available. Skipping correction.');
  return;
end

% copy original image
data.ImageZNoB0Corr = data.ImageZ;

% get ticks for map and image
dataB0 = get_kSpaceAndImageTicks(dataB0, dataB0.AQSlice);
[data, AQSlice] = get_kSpaceAndImageTicks(data, AQSlice);

% calculate deformation in read direction due to B0 map
ReadZOffsetMeter = -dataB0.ImageZ / GradAmpRead;  % dataB0.ImageZ is B0 map in T

% rotate from B0 coordinate system to magnet coordinate system
% For alpha=phi=theta=0 -> read : z , phase(1) : x , phase(2) : y
[RzB0, RxB0, RyB0] = get_aptDegRotationMatrix(dataB0.AQSlice.theta/pi*180, ...
  dataB0.AQSlice.alfa/pi*180, dataB0.AQSlice.phi/pi*180);
% rotate from magnet coordinate system to image coordinate sytem
[RzIm, RxIm, RyIm] = get_aptDegRotationMatrix(AQSlice.theta/pi*180, ...
  AQSlice.alfa/pi*180, AQSlice.phi/pi*180);
% total rotation matrix (passive rotation) of B0 coordinates to image
% coordinates
% R = (RzIm * RyIm * RxIm)';  % <- rotation matrix if B0 map is not rotated
% R = (RxB0' * RyB0' * RzB0')';  % <- rotation matrix if image is not rotated
R = (RzIm * RyIm * RxIm)' * (RxB0' * RyB0' * RzB0')';
% FIXME: Is B0 map always 3d?
[rB0, p1B0, p2B0] = ndgrid(...
  dataB0.Ticks(1).ReadZ + dataB0.AQSlice.Center2OriginImage(3), ...
  dataB0.Ticks(1).PhaseZ + dataB0.AQSlice.Center2OriginImage(1), ...
  dataB0.Ticks(2).PhaseZ + dataB0.AQSlice.Center2OriginImage(2));
CoorB0 = R * [rB0(:).'; p1B0(:).'; p2B0(:).'];
rB0 = reshape(CoorB0(1,:), size(dataB0.ImageZ));
p1B0 = reshape(CoorB0(2,:), size(dataB0.ImageZ));
p2B0 = reshape(CoorB0(3,:), size(dataB0.ImageZ));

ReadZMeterDeformed = rB0 + ReadZOffsetMeter;
scInterp = scatteredInterpolant(rB0(:), p1B0(:), p2B0(:), ...
  double(ReadZMeterDeformed(:)), 'natural', 'none');
[rIm, p1Im, p2Im] = ndgrid(...
  data.Ticks(1).ReadZ + AQSlice.Center2OriginImage(3), ...
  data.Ticks(1).PhaseZ + AQSlice.Center2OriginImage(1), ...
  data.Ticks(2).PhaseZ + AQSlice.Center2OriginImage(2));
ReadZMeterDeformed = scInterp(rIm, p1Im, p2Im);

% correct each line in read direction individually
for t4 = 1:size(data.ImageZ, 4)
  for t3 = 1:size(data.ImageZ, 3)
    for t2 = 1:size(data.ImageZ, 2)
      isInRoi = ~isnan(ReadZMeterDeformed(:,t2,t3));
      if sum(isInRoi) > 2
        % correct position
        data.ImageZ(:,t2,t3,t4) = ...
          interp1(ReadZMeterDeformed(isInRoi,t2,t3), data.ImageZ(isInRoi,t2,t3,t4), data.Ticks(1).ReadZ);
      else
        data.ImageZ(:,t2,t3,t4) = NaN;
      end
    end
  end
end

if AQSlice.CorrectAmplitude
  % gradient of corrected coordinates in read direction
  if isvector(ReadZMeterDeformed)
    ReadZOffsetMeterGrad = gradient(ReadZMeterDeformed);
  else
    % gradient on dimension 1 is 2nd output of gradient
    [~, ReadZOffsetMeterGrad] = gradient(ReadZMeterDeformed);
  end
  ReadZOffsetMeterGrad = ReadZOffsetMeterGrad ./ diff(data.Ticks(1).ReadZ(1:2));

  % Correct image amplitude for stretched or squashed pixel/voxel size
  data.ImageZ = bsxfun(@rdivide, data.ImageZ, ReadZOffsetMeterGrad);
end

end
