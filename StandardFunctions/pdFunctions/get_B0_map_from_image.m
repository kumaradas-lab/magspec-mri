function dataB0 = get_B0_map_from_image(HW, Seq, data1, data2)
%% Get (read corrected) B0 map from image data
%
%   dataB0 = get_B0_map_from_image(HW, Seq, data1, data2)
%
% INPUT:
%
%   HW
%       HW object or structure
%
%   Seq
%       Structure with the following mandatory fields:
%
%     AQSlice
%         Structure with the following mandatory fields:
%             CenterRot, sizeRead, sizePhase, nRead, nPhase, ZeroFillFactor
%
%     CorrectB0Read
%         Structure with the following fields:
%       tReadoutDiff or tEchoIncr
%           If called with 4 input arguments, this is the difference in the
%           readout shift (spin echo) or gradient echo times (FLASH) of the two
%           images.
%       MaxFreqOffset
%           Maximum frequency offset (for selection of RoI) of the B0 map in Hz
%           (default: 2000).
%       MinRelAmp
%           Minimum amplitude (relative to maximum amplitude in image) that is
%           part of the RoI (default: 0.2).
%       MaxRelAmpDiff
%           Maximum relative amplitude deviation (from 1.0) between the images
%           from the two echo times (default: 0.25).
%       RoIExtension
%           Extend the region of interest in all directions by approximately
%           the set value in meter. If this is a scalar, that value applies to
%           all directions. Otherwise, it applies to the dimensions as sorted
%           in data1.ImageZ (see below).
%           (Default: 0.15*image size)
%
%     tEcho
%         Mandatory gradient echo time when called with four arguments. In this
%         case, it is assumed that the image in data1 is the result of a
%         gradient echo measurement (sequence_Flash).
%
%   data1
%       Structure with the following mandatory field:
%
%     ImageZ
%         Array with the complex data of the (first) image as returned by
%         get_kSpaceAndImage.
%
%   data1
%       Oprional structure with the following mandatory field:
%
%     ImageZ
%         Array with the complex data of the second image as returned by
%         get_kSpaceAndImage.
%
%
% OUTPUT:
%
%   dataB0
%       Structure with (among others) the following fields:
%
%     AQSlice
%         Same as Seq.AQSlice(1).
%
%     ImageZ
%         B0 map derived from the input. This B0 map is itself corrected by the
%         read shift caused by the B0 deviations.
%
%
% SEE ALSO: sequence_Flash, sequence_Spin_Echo, correct_read_B0
%
% ------------------------------------------------------------------------------
% (C) Copyright 2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% input check
if nargin < 3
  error('PD:get_B0_map_from_image:TooFewInputArguments', ...
    ['The function "get_B0_map_from_image" must be called at least ', ...
    'with arguments "HW", "Seq" and "data1".']);
end

if ~isfield(Seq, 'AQSlice') || isempty(Seq.AQSlice)
  error('PD:get_B0_map_from_image:MissingFieldSeqAQSlice', ...
    'The input argument "Seq" must have the field "AQSlice".');
end

if isemptyfield(Seq, 'CorrectB0Read')
  error('PD:get_B0_map_from_image:MissingFieldSeqCorrectB0Read', ...
    'The input argument "Seq" must have the field "CorrectB0Read".');
end

if isemptyfield(Seq, 'Read') || isemptyfield(Seq.Read(1), 'GradAmp')
  error('PD:get_B0_map_from_image:MissingFieldSeqReadGradAmp', ...
    'The input argument "Seq" must contain "Read(1).GradAmp".');
end

if ~isfield(data1, 'ImageZ')
  error('PD:get_B0_map_from_image:MissingFieldData1ImageZ', ...
    'The input argument "data1" must have the field "ImageZ".');
end

if nargin > 3
  if ~isfield(data2, 'ImageZ') || any(size(data1.ImageZ) ~= size(data2.ImageZ))
    error('PD:get_B0_map_from_image:InvalidData1ImageZ', ...
    'The input argument "data2" must have the field "ImageZ" with the same size as "data1.ImageZ".');
  end
end


%% default parameters
% maximum frequency offset of B0 map in Hz
if isemptyfield(Seq.CorrectB0Read, 'MaxFreqOffset')
  Seq.CorrectB0Read.MaxFreqOffset = 2000;
end
% minimum amplitude relative to maximum amplitude
if isemptyfield(Seq.CorrectB0Read, 'MinRelAmp')
  Seq.CorrectB0Read.MinRelAmp = 0.2;
end
% maximum relative amplitude deviation from 1.0 between measurements 1 and 2
if isemptyfield(Seq.CorrectB0Read, 'MaxRelAmpDiff')
  Seq.CorrectB0Read.MaxRelAmpDiff = 0.25;
end
% (approximate) extension of the RoI for B0 map in meters
if isemptyfield(Seq.CorrectB0Read, 'RoIExtension')
  Seq.CorrectB0Read.RoIExtension ...
    = 0.15 * [Seq.AQSlice(1).sizeRead, Seq.AQSlice(1).sizePhase(1:2)];
end
if isemptyfield(Seq.CorrectB0Read, 'unwrap'), Seq.CorrectB0Read.unwrap = true; end

dataB0.AQSlice = Seq.AQSlice(1);


%% calculate frequency offset of each voxel
if nargin < 4
  % gradient echo measurement with a single image
  % assume that all phase variation is due to B0 inhomogeneities

  tDephase = Seq.tEcho;

  % phase difference in rad
  if Seq.CorrectB0Read.unwrap
    dataB0.phaseDiff = unwrap3Dmiddle(angle(data1.ImageZ));
  else
    dataB0.phaseDiff = angle(data1.ImageZ);
  end

  % relative amplitude
  dataB0.ampRel = abs(data1.ImageZ./max(data1.ImageZ(:)));

else
  if isfield(Seq.CorrectB0Read, 'tReadoutDiff')
    % spin echo measurements with shifted read-out
    tDephase = Seq.CorrectB0Read.tReadoutDiff;
  elseif isfield(Seq.CorrectB0Read, 'tEchoIncr')
    % gradient echo measurements with different echo times
    tDephase = Seq.CorrectB0Read.tEchoIncr;
  else
    error('PD:get_B0_map_from_image:UnknownDephaseTime', ...
      'Either Seq.CorrectB0Read.tReadoutDiff or Seq.CorrectB0Read.tEchoIncr must be set.');
  end

  % phase difference in rad
  if Seq.CorrectB0Read.unwrap
    dataB0.phaseDiff = unwrap3Dmiddle(angle(data1.ImageZ ./ data2.ImageZ));
  else
    dataB0.phaseDiff = (angle(data1.ImageZ ./ data2.ImageZ));
  end

  % amplitude relation between tEcho1 and tEcho2
  dataB0.ampRel = abs(data2.ImageZ ./ data1.ImageZ);

end

dataB0.frequencyOffset = dataB0.phaseDiff ./ (2*pi*tDephase);
% B0 map in Tesla
% (saved as ImageZ so that correct_read_B0 can be applied to the map)
dataB0.ImageZ = dataB0.phaseDiff / tDephase / HW.GammaDef;


%% get region of interest
% FIXME: Improve selection of RoI
% get ROI
dataB0.RoI = ones(size(dataB0.ImageZ));
% maximum frequency offset
dataB0.RoI(abs(dataB0.frequencyOffset) > Seq.CorrectB0Read.MaxFreqOffset) = 0;
if nargin > 3
  % minimum amplitude measurement tEcho1
  dataB0.RoI(abs(data2.ImageZ) ...
             < max(abs(data2.ImageZ(:))).*Seq.CorrectB0Read.MinRelAmp) = 0;
end
% minimum amplitude measurement tEcho2
dataB0.RoI(abs(data1.ImageZ) ...
           < max(abs(data1.ImageZ(:))).*Seq.CorrectB0Read.MinRelAmp) = 0;
% min relative amplitude difference (tEcho2 > tEcho1)
dataB0.RoI(( dataB0.ampRel - 1) > Seq.CorrectB0Read.MaxRelAmpDiff) = 0;
% min relative amplitude difference (tEcho1 > tEcho2)
dataB0.RoI((1./dataB0.ampRel-1) > Seq.CorrectB0Read.MaxRelAmpDiff) = 0;

% Include points that have many neighbors in ROI
kernel_sz = 3.*(size(dataB0.RoI)>1)+(size(dataB0.RoI)==1);
kernel = ones(kernel_sz);
kernel = kernel / sum(kernel(:));
dataB0.RoI(convn(dataB0.RoI, kernel, 'same') > 0.8) = 1;
% Exclude points that have many neighbors not in ROI
dataB0.RoI(convn(dataB0.RoI, kernel, 'same') < 0.3) = 0;

% figure(220), imagesc(permute(abs(SeqLoop.dataLoop(Loop).ImageZ), [3 1 2]))
% figure(221), imagesc(permute(abs(SeqLoop.dataLoop(Loop-1).ImageZ), [3 1 2]))
% figure(222), imagesc(permute(dataB0.RoI, [3 1 2]))

dataB0.RoI = dataB0.RoI > 0; % convert to logical


%% compensate read shift due to B0 deviation in region of interest
% B0 deviation leads to read offset
% from frequency to meter
dataB0.ReadZOffsetMeter = dataB0.ImageZ./(max(sum(Seq.Read(1).GradAmp.^2,1)).^0.5);
dataB0.ReadZOffsetMeter(~dataB0.RoI) = NaN; % remove points not in ROI
% corrected coordinates in read direction (zero filled)
ticksReadZ = Seq.AQSlice(1).CenterRot(1) + ...
  get_FFTGrid(Seq.AQSlice(1).sizeRead, Seq.AQSlice(1).nRead*Seq.AQSlice(1).ZeroFillFactor(1));
dataB0.ReadZMeterCorr = bsxfun(@minus, ticksReadZ, dataB0.ReadZOffsetMeter);

% copy measured B0 map before readout correction
dataB0.ImageZNoB0Corr = dataB0.ImageZ;

% correct each line in read direction individually with read offset due to B0 error
for t2 = 1:size(dataB0.ImageZ, 2)
  for t3 = 1:size(dataB0.ImageZ, 3)
    isInRoi = ~isnan(dataB0.ReadZMeterCorr(:,t2,t3));
    if sum(isInRoi) > 2
      % correct position
      dataB0.ImageZ(:,t2,t3) = ...
        interp1(dataB0.ReadZMeterCorr(isInRoi,t2,t3), dataB0.ImageZ(isInRoi,t2,t3), ticksReadZ);
    else
      dataB0.ImageZ(:,t2,t3) = NaN;
    end
  end
end


%% optionally extent RoI
if Seq.CorrectB0Read.RoIExtension > 0
  % extend size of ROI by approximately requested value in each direction
  kernel_sz = ceil(Seq.CorrectB0Read.RoIExtension*1.1 ./ ...
                   [Seq.AQSlice(1).sizeRead, Seq.AQSlice(1).sizePhase(1:2)] .* ...
                   [Seq.AQSlice(1).nRead, Seq.AQSlice(1).nPhase(1:2)]*4) .* ...
              (size(dataB0.RoI)>1) + ...
              (size(dataB0.RoI)==1);
  [center_idx{1:3}] = ndgrid(linspace(-1,1,kernel_sz(1)), linspace(-1,1,kernel_sz(2)), linspace(-1,1,kernel_sz(3)));
  kernel = 1 - sqrt(center_idx{1}.^2 + center_idx{2}.^2 + center_idx{3}.^2);
  kernel(kernel < 0) = 0;
  origRoI = ~isnan(dataB0.ImageZ);
  largerRoI = origRoI;
  largerRoI(convn(largerRoI, kernel, 'same') > 0.025*sum(kernel(:))) = 1;
  largerRoI = largerRoI & ~origRoI;

  % extrapolate correction map
  dataB0 = get_kSpaceAndImageTicks(dataB0, dataB0.AQSlice);
  [R, P1, P2] = ndgrid(dataB0.Ticks(1).ReadZ, dataB0.Ticks(1).PhaseZ, dataB0.Ticks(2).PhaseZ);
  scInt = scatteredInterpolant(R(origRoI), P1(origRoI), P2(origRoI), ...
    double(dataB0.ImageZ(origRoI)), 'natural', 'linear');

  dataB0.ImageZ(largerRoI) = scInt(R(largerRoI), P1(largerRoI), P2(largerRoI));
end


% FIXME: Rotate B0 map to magnet coordinate system?


end
