function [data, AQSlice] = plot_kSpaceAndImage(data, AQSlice)
%% Display plots of k-space and image space
%
%       [data, AQSlice] = plot_kSpaceAndImage(data, AQSlice)
%
%
% INPUT:
%
%   data        Structure with image and k-space data (see get_kSpaceAndImage)
%
%   AQSlice     Structure with settings including:
%     iSlice
%                 Index of image to display (fifth linear index of image and
%                 k-space fields in data). (Default: 1)
%
%     RoiRelativeValue, RoiCutOffPercentile
%                 To differentiate between regions with and without sample
%                 inside the image, the region of interest (RoI) is calculated.
%                 This is done by determining the amplitude of the lower
%                 percentile set by AQSlice.RoiCutOffPercentile (default: 0.95).
%                 All pixels/voxels with values >= AQSlice.RoiRelativeValue
%                 (default: 1/5) times that percentile are included inside the
%                 RoI. Isolated outliers, edge and corner pixels/voxels are
%                 removed from the RoI (and similarly for holes inside the RoI).
%                 The distributing effects of zero-filling are taken into
%                 account.
%
%     RoiAbsoluteValue
%                 Amplitude that is used to differentiate between pixels/voxels
%                 inside and outside the region of interest (RoI). This setting
%                 takes precedence over AQSlice.RoiRelativeValue and
%                 AQSlice.RoiCutOffPercentile and is only used if it is defined.
%                 All pixels/voxels with values >= AQSlice.RoiAbsoluteValue are
%                 part of the RoI (no default value).
%
%     sliceomaticProps
%                 Structure with additional argument that are passed to
%                 sliceomatic (default: struct() with no fields).
%
%
% OUTPUT:
%
%   data        Structure with image and k-space data (see corresponding input
%               variable) with added fields. Amongst the added fields are:
%     RoI
%                 Matrix corresponding to the selected slice of data.ImageZ
%                 (or data.ImageSliceomaticZ for 3d images) that is 1 inside the
%                 region of interest and NaN otherwise (see also input
%                 AQSlice.RoiCutOffPercentile, AQSlice.RoiRelativeValue and
%                 AQSlice.RoiAbsoluteValue).
%                 If data.RoICutOff (see below) calculates to 0, data.RoI is
%                 scalar 1.
%                 If data.RoI is already part of the input (with the correct
%                 dimensions or the scalar 1), the RoI is NOT calculated again.
%
%     RoICutOff
%                 Minimum amplitude that is considered to be inside the region
%                 of interest (RoI). Isolated outliers, edge and corner pixels/
%                 voxels are removed. This may lead to pixels/voxels with higher
%                 amplitude lying outside the RoI or pixels/voxels with lower
%                 amplitude lying inside it.
%
%   AQSlice     See corresponding input variable with added fields:
%     plotImagehAxes
%                 1x4 cell with handles to the axes (1: k-space amplitude,
%                 2: image amplitude, 3: k-space phase, 4: image phase)
%
%     plotImageOshAxes
%                 The same as before but for the over-sampled image
%
%     plotB0ppmhAxes
%                 1x3 cell with handles to the axes of the B0 map in ppm
%                 (1: amplitude, 2: phase, 3: gradient)
%
%     plotB0HzhAxes
%                 1x3 cell with handles to the axes of the B0 map in Hz
%                 (1: amplitude, 2: phase, 3: gradient)
%
%     plotB0GradienthAxes
%                 1x3 cell with handles to the axes of the B0 Gradient map in
%                 T/m
%                 (1: tEcho, 2: tAQ)
%
%     plotB1hAxes
%                 Cell with handle to the axes of the B1 map in percent.
%
%     plotDataAxes
%                 1x4 cell with handles to the axes (1: data amplitude,
%                 2: fft1_data amplitude, 3: data phase, 4: fft1_data phase)
%
%   See "get_kSpaceAndImageTicks" for info on additional fields that are added
%   to these structures.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2025 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

%% check input
if numel(AQSlice) > 1
  warning('PD:plot_kSpaceAndImage:MultipleAQSlice', ...
    'Plotting of multiple AQSlices is not supported. First slice is plotted.');
  AQSlice = AQSlice(1);
end

%% default input
AQSlice = set_EmptyField(AQSlice, 'iSlice', 1);  % index of image to display
AQSlice = set_EmptyField(AQSlice, 'CenterRot', [0, 0, 0]);
AQSlice = set_EmptyField(AQSlice, 'plotImageHandle', 3030);
AQSlice = set_EmptyField(AQSlice, 'plotPhase', 0);

AQSlice = set_EmptyField(AQSlice, 'plotImage', 50);
% if AQSlice.plotImage == 1,  AQSlice.plotImage = 50; end
if AQSlice.plotImage == 1,  AQSlice.plotImage = AQSlice.plotImageHandle; end
AQSlice = setPhaseFigure(AQSlice, 'plotImage', 'plotImagePhase', 51);

AQSlice = set_EmptyField(AQSlice, 'plotkSpace', 0);
if AQSlice.plotkSpace == 1, AQSlice.plotkSpace = 52; end
AQSlice = setPhaseFigure(AQSlice, 'plotkSpace', 'plotkSpacePhase', 53);

AQSlice = set_EmptyField(AQSlice, 'plotImageOs', 0);
if AQSlice.plotImageOs == 1,  AQSlice.plotImageOs = 54; end
AQSlice = setPhaseFigure(AQSlice, 'plotImageOs', 'plotImageOsPhase', 55);

AQSlice = set_EmptyField(AQSlice, 'plotB0ppm', 0);
if AQSlice.plotB0ppm == 1,  AQSlice.plotB0ppm = 56; end
AQSlice = set_EmptyField(AQSlice, 'plotB0PpmPhase', 0);
if AQSlice.plotB0PpmPhase == 1
  if ishghandle(AQSlice.plotB0ppm, 'uipanel') || AQSlice.plotB0ppm == 0
    AQSlice.plotB0PpmPhase = 57;
  else
    AQSlice.plotB0PpmPhase = double(AQSlice.plotB0ppm) + 1;
  end
end
AQSlice = set_EmptyField(AQSlice, 'plotB0ppmGradient', 0);
if AQSlice.plotB0ppmGradient == 1
  if ishghandle(AQSlice.plotB0ppm, 'uipanel') || AQSlice.plotB0ppm == 0
    AQSlice.plotB0ppmGradient = 58;
  else
    AQSlice.plotB0ppmGradient = double(AQSlice.plotB0ppm) + 2;
  end
end

AQSlice = set_EmptyField(AQSlice, 'plotB0Hz', 0);
if AQSlice.plotB0Hz == 1,  AQSlice.plotB0Hz = 59; end
AQSlice = set_EmptyField(AQSlice, 'plotB0HzPhase', 0);
if AQSlice.plotB0HzPhase == 1
  if ishghandle(AQSlice.plotB0Hz, 'uipanel') || AQSlice.plotB0Hz == 0
    AQSlice.plotB0HzPhase = 60;
  else
    AQSlice.plotB0HzPhase = double(AQSlice.plotB0Hz) + 1;
  end
end
AQSlice = set_EmptyField(AQSlice, 'plotB0HzGradient', 0);
if AQSlice.plotB0HzGradient == 1
  if ishghandle(AQSlice.plotB0Hz, 'uipanel') || AQSlice.plotB0Hz == 0
    AQSlice.plotB0HzGradient = 60;
  else
    AQSlice.plotB0HzGradient = double(AQSlice.plotB0Hz) + 2;
  end
end

AQSlice = set_EmptyField(AQSlice, 'plotB1percent', 0);
if AQSlice.plotB1percent == 1,  AQSlice.plotB1percent = 62; end

AQSlice = set_EmptyField(AQSlice, 'plotB0Gradient', 0);
if AQSlice.plotB0Gradient == 1,  AQSlice.plotB0Gradient = 63; end
AQSlice = set_EmptyField(AQSlice, 'plotB0GradientPhase', 0);
if AQSlice.plotB0GradientPhase == 1,  AQSlice.plotB0GradientPhase = 64; end
AQSlice = set_EmptyField(AQSlice, 'plotData', 0);
if AQSlice.plotData == 1,  AQSlice.plotData = 65; end
AQSlice = set_EmptyField(AQSlice, 'plotFft1_data', 0);
if AQSlice.plotFft1_data == 1, AQSlice.plotFft1_data = 65; end

if isemptyfield(AQSlice, 'angle2Turns'), AQSlice.angle2Turns = 1/(2*pi); end

AQSlice = set_EmptyField(AQSlice, 'raiseFigures', true);
AQSlice = set_EmptyField(AQSlice, 'AmplitudeUnit', 'T');
AQSlice = set_EmptyField(AQSlice, 'AmplitudeUnitScale', 1);
AQSlice = set_EmptyField(AQSlice, 'LengthUnit', 'mm');
AQSlice = set_EmptyField(AQSlice, 'LengthUnitScale', 1e-3);
AQSlice = set_EmptyField(AQSlice, 'ZeroFillFactor', 1);
AQSlice = set_EmptyField(AQSlice, 'ZeroFillWindowSize', Inf);
AQSlice = set_EmptyField(AQSlice, 'RoiRelativeValue', 1/5);
AQSlice = set_EmptyField(AQSlice, 'RoiCutOffPercentile', .95);

AQSlice = set_EmptyField(AQSlice, 'plotImagehAxes', cell(1, 4));
AQSlice = set_EmptyField(AQSlice, 'plotDatahAxes', cell(1, 4));
AQSlice = set_EmptyField(AQSlice, 'plotImageOshAxes', cell(1, 4));
AQSlice = set_EmptyField(AQSlice, 'plotB0ppmhAxes', cell(1, 3));
AQSlice = set_EmptyField(AQSlice, 'plotB0HzhAxes', cell(1, 3));
AQSlice = set_EmptyField(AQSlice, 'plotB0GradienthAxes', cell(1,3));
AQSlice = set_EmptyField(AQSlice, 'plotB1hAxes', cell(1));

AQSlice = set_EmptyField(AQSlice, 'sliceomaticProps', struct());  % structure that is passed to sliceomatic

data = set_EmptyField(data, 'ImageZ', data.Image);
% data = set_EmptyField(data, 'kSpaceZ', data.kSpace);
% data = set_EmptyField(data, 'kSpaceOsZ', data.kSpaceOsRaw);

if AQSlice.iSlice < 1
  AQSlice.iSlice = 1;
  warning('PD:plot_kSpaceAndImage:iSlice', ...
    'iSlice < 1. First slice is plotted.');
end
szImage = arrayfun(@(x) size(data.ImageZ, x), 1:6);
if AQSlice.iSlice > prod([szImage(5), szImage(6)])
  AQSlice.iSlice = prod([szImage(5), szImage(6)]);
  warning('PD:plot_kSpaceAndImage:iSlice', ...
    ['iSlice > prod([size(data.ImageZ,5), size(data.ImageZ,6)]). ', ...
    'Last slice (iSlice=%d) is plotted.'], ...
    AQSlice.iSlice);
end


%% grids
[data, AQSlice] = get_kSpaceAndImageTicks(data, AQSlice);


%% display image and k-space
if sum(szImage(1:4)>1) >= 3
  %% 3d image
  % (FIXME: This does not handle 4d images completely)

  % clean up parent before plotting axes
  if ishghandle(AQSlice.plotImageHandle, 'uipanel')
    deleteChildren = true;
    if isappdata(AQSlice.plotImageHandle, 'sliceomatic')
      d = getappdata(AQSlice.plotImageHandle, 'sliceomatic');
      if isvalid(d.axmain_obj)
        % keep the sliceomatic and re-use it
        deleteChildren = false;
      end
    end
    if deleteChildren
      delete(get(AQSlice.plotImageHandle, 'Children'));
    end
  end

  if sum(AQSlice.nPhase(:) > 1) < 3
    % read encoding for first dimension
    ZeroFillFactorPermuted = data.ZeroFillFactor(data.PermuteOrder);
  else
    % phase encoding for all 3 dimensions (CSI)
    ZeroFillFactorPermuted = data.ZeroFillFactor(data.PermuteOrder+1);
  end

  % data.ImageSliceomatic = permute(squeeze(data.Image), data.PermuteOrder);
  % data.kSpaceSliceomatic = permute(squeeze(data.kSpace), data.PermuteOrder);
  data.ImageSliceomaticZ = permute(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice)), data.PermuteOrder);
  % data.kSpaceSliceomaticZ = permute(squeeze(data.kSpaceZ), data.PermuteOrder);
  % data.kSpaceSliceomaticOsZ = permute(squeeze(data.kSpaceOsZ), data.PermuteOrder);
  % data.kSpaceSliceomaticOs = permute(squeeze(data.kSpaceOs), data.PermuteOrder);
  if isfield(data, 'ImageZFrequency')
    data.ImageFrequencySliceomaticZ = permute(squeeze(data.ImageZFrequency(:,:,:,:,AQSlice.iSlice)), data.PermuteOrder); %KPM CSI
  end

  % find Region of Interest
  if isemptyfield(data, 'RoI') || ...
      (~(isscalar(data.RoI) && data.RoI == 1) && ...
      ~all(size(data.ImageSliceomaticZ) == size(data.RoI)))
    % FIXME: Would it be better if this were a Boolean array (memory footprint)
    data = smoothRoI(data, AQSlice, data.ImageSliceomaticZ, 2*ZeroFillFactorPermuted+1);
    data.RoI = double(data.RoI);  % There is no logical NaN.
    data.RoI(~data.RoI) = NaN;
  end

  if ~isemptyfield(data, 'f_fft1_data')
    firstkLineIdx = find(AQSlice.UseAQWindow(:,1,1)~=0 & AQSlice.UsetRep(:,1,1)~=0, 1, 'first');
    data.fCenter = data.f_fft1_data(floor(AQSlice.nRead*AQSlice.ReadOS/2)+1,AQSlice.UseAQWindow(firstkLineIdx),AQSlice.UsetRep(firstkLineIdx));
  else
    data.fCenter = [];
  end

  if AQSlice.plotImage ~= 0
    AQSlice.plotImagehAxes{2} = plot3Ddata(AQSlice.plotImage, data, AQSlice, ...
      abs(data.ImageSliceomaticZ)/AQSlice.AmplitudeUnitScale, true, 'ReadZ', 'PhaseZ');
    title(AQSlice.plotImagehAxes{2}.hAxes, ['Image Amplitude in ' AQSlice.AmplitudeUnit ]);
  end

  if AQSlice.plotImagePhase ~= 0
    AQSlice.plotImagehAxes{4} = plot3Ddata(AQSlice.plotImagePhase, data, AQSlice, ...
      angle(data.ImageSliceomaticZ), true, 'ReadZ', 'PhaseZ');
    title(AQSlice.plotImagehAxes{4}.hAxes, 'Image Phase in rad');
  end

  if AQSlice.plotImageOs ~= 0
    AQSlice.plotImageOshAxes{2} = plot3Ddata(AQSlice.plotImageOs, data, AQSlice, ...
      abs(permute(squeeze(data.ImageOsZ(:,:,:,:,AQSlice.iSlice)), data.PermuteOrder)), true, 'ReadOsZ', 'PhaseOsZ');
    title(AQSlice.plotImageOshAxes{2}.hAxes, ['Over-Sampled Image Amplitude in ' AQSlice.AmplitudeUnit ]);
  end

  if AQSlice.plotImageOsPhase ~= 0
    AQSlice.plotImageOshAxes{4} = plot3Ddata(AQSlice.plotImageOsPhase, data, AQSlice, ...
      angle(permute(squeeze(data.ImageOsZ(:,:,:,:,AQSlice.iSlice)), data.PermuteOrder)), true, 'ReadOsZ', 'PhaseOsZ');
    title(AQSlice.plotImageOshAxes{4}.hAxes, 'Over-Sampled Image Phase in rad');
  end

  if AQSlice.plotB0ppm ~= 0
    AQSlice.plotB0ppmhAxes{1} = plot3Ddata(AQSlice.plotB0ppm, data, AQSlice, ...
      data.RoI.*unwrap3Dmiddle(-angle(data.ImageSliceomaticZ))/pi/2/AQSlice.tEcho/data.fCenter*1e6, true, 'ReadZ', 'PhaseZ');
    title(AQSlice.plotB0ppmhAxes{1}.hAxes, ...
      ['Offset to = ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT, in ppm']);
  end

  if AQSlice.plotB0Hz ~= 0
    AQSlice.plotB0HzhAxes{1} = plot3Ddata(AQSlice.plotB0Hz, data, AQSlice, ...
      data.RoI.*unwrap3Dmiddle(-angle(data.ImageSliceomaticZ))/pi/2/AQSlice.tEcho, true, 'ReadZ', 'PhaseZ');
    title(AQSlice.plotB0HzhAxes{1}.hAxes, ...
      ['Offset to = ' num2str(data.fCenter*1e-6, '%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT, in Hz']);
  end

  if AQSlice.plotB0HzPhase ~= 0 && isfield(data, 'ImageZFrequency')
    AQSlice.plotB0HzhAxes{2} = plot3Ddata(AQSlice.plotB0HzPhase, data, AQSlice, ...
      data.RoI.*data.ImageFrequencySliceomaticZ, true, 'ReadZ', 'PhaseZ');
    title(AQSlice.plotB0HzhAxes{2}.hAxes, ...
      ['Offset to = ' num2str(data.fCenter*1e-6, '%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT, in Hz']);
  end

  if AQSlice.plotB0PpmPhase ~= 0 && isfield(data, 'ImageZFrequency')
    AQSlice.plotB0ppmhAxes{2} = plot3Ddata(AQSlice.plotB0PpmPhase, data, AQSlice, ...
      data.RoI.*data.ImageFrequencySliceomaticZ./data.fCenter.*1e6, true, 'ReadZ', 'PhaseZ');
    title(AQSlice.plotB0ppmhAxes{2}.hAxes, ...
      ['Offset to = ' num2str(data.fCenter*1e-6, '%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT, in ppm']);
  end

  if AQSlice.plotB0HzGradient ~= 0
    B0Hz = data.RoI.*unwrap3Dmiddle(-angle(data.ImageSliceomaticZ))/pi/2/AQSlice.tEcho;
    hParent = AQSlice.plotB0HzGradient;
    if ~ishghandle(hParent, 'uipanel') && ...
        (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
      hParent = figure(hParent);
    end
    if AQSlice.nRead>1
      [B0HzGradient_xyz{[2,1,3]}] = gradient(B0Hz, diff(data.Ticks(1).ReadZ(1:2)), diff(data.Ticks(1).PhaseZ(1:2)), diff(data.Ticks(2).PhaseZ(1:2)));
      B0HzGradient_xyz = cat(4, B0HzGradient_xyz{:});
      B0HzGradient = sqrt(sum(B0HzGradient_xyz.^2, 4));
    else
      [B0HzGradient_xyz{[2,1,3]}] = gradient(B0Hz, diff(data.Ticks(3).PhaseZ(1:2)), diff(data.Ticks(1).PhaseZ(1:2)), diff(data.Ticks(2).PhaseZ(1:2)));
      B0HzGradient_xyz = cat(4, B0HzGradient_xyz{:});
      B0HzGradient = sqrt(sum(B0HzGradient_xyz.^2, 4));
    end
    AQSlice.plotB0HzhAxes{3} = plot3Ddata(hParent, data, AQSlice, B0HzGradient, true, 'ReadZ', 'PhaseZ');
    clear B0Hz B0HzGradient_xyz B0HzGradient
    title(AQSlice.plotB0HzhAxes{3}.hAxes, ...
      ['Gradient at = ' num2str(data.fCenter*1e-6, '%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT, in Hz/m']);
  end

  if AQSlice.plotB0ppmGradient ~= 0
    B0ppm = data.RoI .* unwrap3Dmiddle(-angle(data.ImageSliceomaticZ))/pi/2/AQSlice.tEcho/data.fCenter*1e6;
    hParent = AQSlice.plotB0ppmGradient;
    if ~ishghandle(hParent, 'uipanel') && ...
        (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
      hParent = figure(hParent);
    end
    if AQSlice.nRead>1
      [B0ppmGradient_xyz{[2,1,3]}] = gradient(B0ppm, diff(data.Ticks(1).ReadZ(1:2)), diff(data.Ticks(1).PhaseZ(1:2)), diff(data.Ticks(2).PhaseZ(1:2)));
      B0ppmGradient_xyz = cat(4, B0ppmGradient_xyz{:});
      B0ppmGradient = sqrt(sum(B0ppmGradient_xyz.^2, 4));
    else
      [B0ppmGradient_xyz{[2,1,3]}] = gradient(B0ppm, diff(data.Ticks(3).PhaseZ(1:2)), diff(data.Ticks(1).PhaseZ(1:2)), diff(data.Ticks(2).PhaseZ(1:2)));
      B0ppmGradient_xyz = cat(4, B0ppmGradient_xyz{:});
      B0ppmGradient = sqrt(sum(B0ppmGradient_xyz.^2, 4));
    end
    AQSlice.plotB0ppmhAxes{3} = plot3Ddata(hParent, data, AQSlice, B0ppmGradient, true, 'ReadZ', 'PhaseZ');
    clear B0ppm B0ppmGradient_xyz B0ppmGradient
    title(AQSlice.plotB0ppmhAxes{3}.hAxes, ...
      ['Gradient at = ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT, in ppm/m']);
  end

  if AQSlice.plotB0Gradient ~= 0
    B0 = data.RoI ...
      .* unwrap3Dmiddle(-angle(data.ImageSliceomaticZ)) ...
      / AQSlice.tEcho / AQSlice.Gamma;
    hParent = AQSlice.plotB0Gradient;
    if ~ishghandle(hParent, 'uipanel') && ...
        (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
      hParent = figure(hParent);
    end
    if AQSlice.nRead>1
      [B0Gradient_xyz{[2,1,3]}] = gradient(B0, diff(data.Ticks(1).ReadZ(1:2)), diff(data.Ticks(1).PhaseZ(1:2)), diff(data.Ticks(2).PhaseZ(1:2)));
      B0Gradient_mean(3) = mean(reshape(B0Gradient_xyz{3}, [], 1), 'omitnan');
      B0Gradient_mean(2) = mean(reshape(B0Gradient_xyz{2}, [], 1), 'omitnan');
      B0Gradient_mean(1) = mean(reshape(B0Gradient_xyz{1}, [], 1), 'omitnan');
      B0Gradient_xyz = cat(4, B0Gradient_xyz{:});
      B0Gradient = sqrt(sum(B0Gradient_xyz.^2, 4));
    else
      [B0Gradient_xyz{[2,1,3]}] = gradient(B0, diff(data.Ticks(3).PhaseZ(1:2)), diff(data.Ticks(1).PhaseZ(1:2)), diff(data.Ticks(2).PhaseZ(1:2)));
      B0Gradient_mean(3) = mean(reshape(B0Gradient_xyz{3}, [], 1), 'omitnan');
      B0Gradient_mean(2) = mean(reshape(B0Gradient_xyz{2}, [], 1), 'omitnan');
      B0Gradient_mean(1) = mean(reshape(B0Gradient_xyz{1}, [], 1), 'omtinan');
      B0Gradient_xyz = cat(4, B0Gradient_xyz{:});
      B0Gradient = sqrt(sum(B0Gradient_xyz.^2, 4));
    end
    AQSlice.plotB0GradienthAxes{1} = plot3Ddata(hParent, data, AQSlice, B0Gradient*1e3, true, 'ReadZ', 'PhaseZ');
    clear B0 B0Gradient_xyz B0Gradient
    if (AQSlice.nPhase(:)>1) == 3  % 3d image
      xlabel(AQSlice.plotB0GradienthAxes{1}.hAxes, ...
        [AQSlice.PhaseCartesianAxis{3} 'phase(3) in ' AQSlice.LengthUnit ' (mean ' num2str(B0Gradient_mean(2)*1e3,'%10.3f') ' mT/m)']);
    else
      xlabel(AQSlice.plotB0GradienthAxes{1}.hAxes, ...
        [AQSlice.ReadCartesianAxis{1} 'read in ' AQSlice.LengthUnit ' (mean ' num2str(B0Gradient_mean(2)*1e3,'%10.3f') ' mT/m)']);
    end
    ylabel(AQSlice.plotB0GradienthAxes{1}.hAxes, ...
      [AQSlice.PhaseCartesianAxis{1} 'phase(1) in ' AQSlice.LengthUnit ' (mean ' num2str(B0Gradient_mean(1)*1e3,'%10.3f') ' mT/m)']);
    zlabel(AQSlice.plotB0GradienthAxes{1}.hAxes, ...
      [AQSlice.PhaseCartesianAxis{2} 'phase(2) in ' AQSlice.LengthUnit ' (mean ' num2str(B0Gradient_mean(3)*1e3,'%10.3f') ' mT/m)']);
    title(AQSlice.plotB0GradienthAxes{1}.hAxes, ...
      ['Gradient at = ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT, in mT/m']);
  end

  if AQSlice.plotB0GradientPhase ~= 0 && isfield(data, 'ImageZFrequency')
    B0 = data.RoI .* data.ImageFrequencySliceomaticZ*2*pi/AQSlice.Gamma;
    hParent = AQSlice.plotB0GradientPhase;
    if ~ishghandle(hParent, 'uipanel') && ...
        (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
      hParent = figure(hParent);
    end
    [B0Gradient_xyz{[2,1,3]}] = gradient(B0, diff(data.Ticks(3).PhaseZ(1:2)), diff(data.Ticks(1).PhaseZ(1:2)), diff(data.Ticks(2).PhaseZ(1:2)));
    B0Gradient_mean(3) = mean(reshape(B0Gradient_xyz{3}, [], 1), 'omitnan');
    B0Gradient_mean(2) = mean(reshape(B0Gradient_xyz{2}, [], 1), 'omitnan');
    B0Gradient_mean(1) = mean(reshape(B0Gradient_xyz{1}, [], 1), 'omitnan');
    B0Gradient_xyz = cat(4, B0Gradient_xyz{:});
    B0Gradient = sqrt(sum(B0Gradient_xyz.^2, 4));
    AQSlice.plotB0GradienthAxes{2} = plot3Ddata(hParent, data, AQSlice, B0Gradient*1e3, true, 'ReadZ', 'PhaseZ');
    clear B0 B0Gradient_xyz B0Gradient
    xlabel(AQSlice.plotB0GradienthAxes{2}.hAxes, ...
      [AQSlice.PhaseCartesianAxis{3} 'phase(3) in ' AQSlice.LengthUnit ' (mean ' num2str(B0Gradient_mean(2)*1e3,'%10.3f') ' mT/m)']);
    ylabel(AQSlice.plotB0GradienthAxes{2}.hAxes, ...
      [AQSlice.PhaseCartesianAxis{1} 'phase(1) in ' AQSlice.LengthUnit ' (mean ' num2str(B0Gradient_mean(1)*1e3,'%10.3f') ' mT/m)']);
    zlabel(AQSlice.plotB0GradienthAxes{2}.hAxes, ...
      [AQSlice.PhaseCartesianAxis{2} 'phase(2) in ' AQSlice.LengthUnit ' (mean ' num2str(B0Gradient_mean(3)*1e3,'%10.3f') ' mT/m)']);
    title(AQSlice.plotB0GradienthAxes{2}.hAxes, ...
      ['Gradient at = ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT, in mT/m']);
  end


  if AQSlice.plotkSpace ~= 0
    % relative amplitude to 1uVeff of the ADC of k-Space in dB
    AQSlice.plotImagehAxes{1} = plot3Ddata(AQSlice.plotkSpace, data, AQSlice, ...
      max(-40, 20*log10(abs(permute(squeeze(data.kSpaceOs(:,:,:,:,AQSlice.iSlice)), data.PermuteOrder)) * ...
                        AQSlice.VoxelVolume / AQSlice.AreaCoil * data.Amplitude2Uin(AQSlice.UseAQWindow(1),1)/2^0.5*1e6)), ...
      false, 'ReadOs', 'PhaseOs');
    title(AQSlice.plotImagehAxes{1}.hAxes, ['k-Space in dB' char(181) 'Veff']);
  end

  if AQSlice.plotkSpacePhase ~= 0
    AQSlice.plotImagehAxes{3} = plot3Ddata(AQSlice.plotkSpacePhase, data, AQSlice, ...
      angle(permute(squeeze(data.kSpaceOs(:,:,:,:,AQSlice.iSlice)), data.PermuteOrder)), false, 'ReadOs', 'PhaseOs');
    title(AQSlice.plotImagehAxes{3}.hAxes, 'k-Space Phase in rad');
  end

  if AQSlice.plotB1percent ~= 0
    % At very small excitation angles (sin(x) ~ x), the combined effect of
    % sending and receiving at a B1 offset is ~ (sqrt(abs(x))/abs(center)*100
    dataB1percent = (sqrt(abs(data.ImageSliceomaticZ)) / sqrt(abs(data.ImageSliceomaticZ(round(end/2),round(end/2),round(end/2))))*100);
    % B1percent_min = 50;
    % B1percent_max = 150;
    % dataB1percent(dataB1percent < B1percent_min) = NaN;
    % dataB1percent(dataB1percent > B1percent_max) = NaN;
    AQSlice.plotB1hAxes{1} = plot3Ddata(AQSlice.plotB1percent, data, AQSlice, dataB1percent, true, 'ReadZ', 'PhaseZ');
    title(AQSlice.plotB1hAxes{1}.hAxes, 'B1+ intensity relative to center in percent');
    set(AQSlice.plotB1hAxes{1}, 'CLim', [50 150]);
  end

elseif sum(szImage(1:4)>1) == 2
  %% 2d image
  if (AQSlice.plotImage ~= 0) || (AQSlice.plotImageOs ~= 0) || (AQSlice.plotkSpace ~= 0) || ...
      (AQSlice.plotB0ppm ~= 0) || (AQSlice.plotB0Hz ~= 0)
    % find or create parents
    hParent{1} = AQSlice.plotImageHandle;
    if ~ishghandle(hParent{1}, 'uipanel') && ...
        (AQSlice.raiseFigures || ~ishghandle(hParent{1}, 'figure'))
      hParent{1} = figure(hParent{1});
    end
    if (AQSlice.plotImage ~= 0) && (AQSlice.plotImageOs ~= 0)
      % open second figure
      hParent{2} = AQSlice.plotImageOs;
      if ~ishghandle(hParent{2}, 'uipanel')
        if AQSlice.raiseFigures || ~ishghandle(hParent{2}, 'figure')
          hParent{2} = figure(hParent{2});
        end
        % FIXME: Implement re-using axes in second figure
        hParent{2} = clf(hParent{2});
      end
    end

    % find axes
    foundAxes = false(1, 4);
    hKids = get(hParent{1}, 'Children');
    hKids(~isvalid(hKids)) = [];
    if ~isempty(hKids)
      allTags = get(hKids, 'Tag');
      if ~iscell(allTags), allTags = {allTags}; end
      isForeign = cellfun(@isempty, regexp(allTags, 'kSpaceAndImage_.*_2d'));
      hAxesKids = hKids(~isForeign);
      delete(hKids(isForeign));

      % check if we have all necessary axes
      axesTags = {'kSpaceAndImage_kspace_amp_2d', 'kSpaceAndImage_image_amp_2d', ...
        'kSpaceAndImage_kspace_phase_2d', 'kSpaceAndImage_image_phase_2d'};
      actualAxesTags = get(hAxesKids, 'Tag');
      hAxes = zeros(1, numel(axesTags));
      for iAxes = 1:numel(foundAxes)
        matchedAxes = strcmp(axesTags{iAxes}, actualAxesTags);
        if sum(matchedAxes) == 1
          foundAxes(iAxes) = true;
          hAxes(iAxes) = hAxesKids(matchedAxes);
        end
      end
    end

    if all(foundAxes)
      % assign axes
      AQSlice.plotImagehAxes = num2cell(hAxes);
    else
      if ~isempty(hKids)
        delete(hKids(isvalid(hKids)));
      end
      % create all axes
      AQSlice.plotImagehAxes{1} = subplot(2,2,1, 'Parent', hParent{1}, ...
        'Tag', 'kSpaceAndImage_kspace_amp_2d', 'Visible', 'off');
      AQSlice.plotImagehAxes{2} = subplot(2,2,2, 'Parent', hParent{1}, ...
        'Tag', 'kSpaceAndImage_image_amp_2d', 'Visible', 'off');
      AQSlice.plotImagehAxes{3} = subplot(2,2,3, 'Parent', hParent{1}, ...
        'Tag', 'kSpaceAndImage_kspace_phase_2d', 'Visible', 'off');
      AQSlice.plotImagehAxes{4} = subplot(2,2,4, 'Parent', hParent{1}, ...
        'Tag', 'kSpaceAndImage_image_phase_2d', 'Visible', 'off');
    end

    figName{1} = '';
    iOs = 1;
    imageTypes = {};
    % set visibility for axes and position them correctly
    if (AQSlice.plotImage ~= 0) || (AQSlice.plotImageOs ~= 0)
      if AQSlice.plotImage ~= 0
        iOs = iOs+1;
        imageTypes{end+1} = 'plotImagehAxes';
        if AQSlice.plotkSpace ~= 0
          if AQSlice.plotPhase
            figName{1} = 'Magnitude and Phase of k-Space, Magnitude and Phase of Image';
            for iAxes = 1:4
              subplot(2,2,iAxes, AQSlice.plotImagehAxes{iAxes}, 'Visible', 'on');
            end
          else
            figName{1} = 'k-Space and Image';
            for iAxes = 1:2
              subplot(1,2,iAxes, AQSlice.plotImagehAxes{iAxes}, 'Visible', 'on');
            end
            for iAxes = 3:4
              cla(AQSlice.plotImagehAxes{iAxes});
              set(AQSlice.plotImagehAxes{iAxes}, 'Visible', 'off');
            end
          end
        else
          if AQSlice.plotPhase
            figName{1} = 'Magnitude and Phase of Image';
            visAxes = [2, 4];
            for iAxes = 1:2
              subplot(2,1,iAxes, AQSlice.plotImagehAxes{visAxes(iAxes)}, 'Visible', 'on');
            end
            for iAxes = [1, 3]
              cla(AQSlice.plotImagehAxes{iAxes});
              set(AQSlice.plotImagehAxes{iAxes}, 'Visible', 'off');
            end
          else
            figName{1} = 'Image';
            subplot(1,1,1, AQSlice.plotImagehAxes{2}, 'Visible', 'on');
            for iAxes = [1, 3:4]
              cla(AQSlice.plotImagehAxes{iAxes});
              set(AQSlice.plotImagehAxes{iAxes}, 'Visible', 'off');
            end
          end
        end
      end
      % FIXME: Implement re-using of graphics objects in second figure
      if AQSlice.plotImageOs ~= 0
        imageTypes{end+1} = 'plotImageOshAxes';
        if AQSlice.plotkSpace ~= 0
          if AQSlice.plotPhase
            figName{iOs} = 'Magnitude and Phase of k-Space, Magnitude and Phase of Over-Sampled Image';
            AQSlice.plotImageOshAxes{1} = subplot(2,2,1, 'Parent', hParent{iOs});
            AQSlice.plotImageOshAxes{2} = subplot(2,2,2, 'Parent', hParent{iOs});
            AQSlice.plotImageOshAxes{3} = subplot(2,2,3, 'Parent', hParent{iOs});
            AQSlice.plotImageOshAxes{4} = subplot(2,2,4, 'Parent', hParent{iOs});
          else
            figName{iOs} = 'k-Space and Over-Sampled Image';
            AQSlice.plotImageOshAxes{1} = subplot(1,2,1, 'Parent', hParent{iOs});
            AQSlice.plotImageOshAxes{2} = subplot(1,2,2, 'Parent', hParent{iOs});
          end
        else
          if AQSlice.plotPhase
            figName{iOs} = 'Magnitude and Phase of Over-Sampled Image';
            AQSlice.plotImageOshAxes{2} = subplot(2,1,1, 'Parent', hParent{iOs});
            AQSlice.plotImageOshAxes{4} = subplot(2,1,2, 'Parent', hParent{iOs});
          else
            figName{iOs} = 'Over-sampled Image';
            AQSlice.plotImageOshAxes{2} = axes('Parent', hParent{iOs});
          end
        end
      else
        iOs = 1;
      end
    elseif AQSlice.plotkSpace ~= 0
      imageTypes{end+1} = 'plotImagehAxes';
      if AQSlice.plotPhase
        figName{1} = 'Magnitude and Phase of k-Space';
        visAxes = [1, 3];
        for iAxes = 1:2
          subplot(2,1,iAxes, AQSlice.plotImagehAxes{visAxes(iAxes)}, 'Visible', 'on');
        end
        for iAxes = [2, 4]
          cla(AQSlice.plotImagehAxes{iAxes});
          set(AQSlice.plotImagehAxes{iAxes}, 'Visible', 'off');
        end
      else
        figName{1} = 'k-Space';
        subplot(1,1,1, AQSlice.plotImagehAxes{1}, 'Visible', 'on');
        for iAxes = 2:4
          cla(AQSlice.plotImagehAxes{iAxes});
          set(AQSlice.plotImagehAxes{iAxes}, 'Visible', 'off');
        end
      end
    end
    for jOs = 1:iOs
      if ~ishghandle(hParent{jOs}, 'uipanel') && ~isempty(figName{jOs})
        set(hParent{jOs}, 'Name', figName{jOs});
      end
    end

    if AQSlice.plotkSpace ~= 0
      for imageType = imageTypes
        hKids = get(AQSlice.(imageType{1}){1}, 'Children');
        hImg = findobj(hKids, 'flat', 'Tag', 'kSpaceAndImage_2d_image');
        cData = 20*log10((abs(squeeze(data.kSpaceOs(:,:,:,:,AQSlice.iSlice)) * ...
          AQSlice.VoxelVolume / AQSlice.AreaCoil * data.Amplitude2Uin(AQSlice.UseAQWindow(1),1)/2^0.5)./1e-6));
        if sum(AQSlice.nPhase==1) == 2
          % read-phase encoded
          cData = permute(cData, [2,1]);
        else
          % phase-phase encoded
        end
        if isempty(hImg)
          imagesc(cData, ...
            'Parent', AQSlice.(imageType{1}){1}, 'Tag', 'kSpaceAndImage_2d_image');
        else
          set(hImg, 'CData', cData);
          set(AQSlice.(imageType{1}){1}, 'CLimMode', 'auto');
          axis(AQSlice.(imageType{1}){1}, 'tight');
        end
        set(AQSlice.(imageType{1}){1}, 'YDir', 'normal');
        title(AQSlice.(imageType{1}){1}, ['k-Space in dB' char(181) 'V']);
        k_clim = [-40,0] + [0,1].*get(AQSlice.(imageType{1}){1}, 'CLim');
        if k_clim(1) < k_clim(2), set(AQSlice.(imageType{1}){1}, 'CLim', k_clim); end
        if sum(AQSlice.nPhase==1) == 2
          % read-phase encoded
          t = find(AQSlice.nPhase>1, 1, 'first');
          xlabel(AQSlice.(imageType{1}){1}, [AQSlice.ReadCartesianAxis{1} 'read in samples']);
          ylabel(AQSlice.(imageType{1}){1}, [AQSlice.PhaseCartesianAxis{t} 'phase(' num2str(t(1)) ') in samples']);
          aspect=[AQSlice.nRead/AQSlice.nPhase(find(AQSlice.nPhase>1, 1, 'first')), ...
            AQSlice.sizeRead/AQSlice.sizePhase(find(AQSlice.nPhase>1, 1, 'first')), 1];
        else
          t = find(AQSlice.nPhase>1, 2, 'first');
          xlabel(AQSlice.(imageType{1}){1}, [AQSlice.PhaseCartesianAxis{t(2)} 'phase(' num2str(t(2)) ') in samples']);
          ylabel(AQSlice.(imageType{1}){1}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in samples']);
          aspect = [AQSlice.nPhase(t(2))    / AQSlice.nPhase(t(1)), ...
                    AQSlice.sizePhase(t(2)) / AQSlice.sizePhase(t(1)), 1];
        end
        if any([isnan(aspect), isinf(aspect), aspect==0])
          set(AQSlice.(imageType{1}){1}, 'PlotBoxAspectRatioMode', 'auto');
        else
          set(AQSlice.(imageType{1}){1}, 'PlotBoxAspectRatio', aspect);
        end
      end
    end
    if AQSlice.plotPhase && (AQSlice.plotkSpace ~= 0)
      for imageType = imageTypes
        hKids = get(AQSlice.(imageType{1}){3}, 'Children');
        hImg = findobj(hKids, 'flat', 'Tag', 'kSpaceAndImage_2d_image');
        cData = angle(squeeze(data.kSpaceOs(:,:,:,:,AQSlice.iSlice)));
        if sum(AQSlice.nPhase==1) == 2
          % read-phase encoded
          cData = permute(cData, [2, 1]);
        else
          % phase-phase encoded
        end
        if isempty(hImg)
          imagesc(cData, ...
            'Parent', AQSlice.(imageType{1}){3}, 'Tag', 'kSpaceAndImage_2d_image');
        else
          set(hImg, 'CData', cData);
          axis(AQSlice.(imageType{1}){3}, 'tight');
        end
        set(AQSlice.(imageType{1}){3}, 'CLim', [-pi,pi]);
        set(AQSlice.(imageType{1}){3}, 'YDir', 'normal');
        colormap(AQSlice.(imageType{1}){3}, gray);
        if sum([AQSlice.nPhase(1)>1, AQSlice.nPhase(2)>1, AQSlice.nPhase(3)>1])==1
          aspect = [AQSlice.nRead    / AQSlice.nPhase(find(AQSlice.nPhase>1,1,'first')), ...
                    AQSlice.sizeRead / AQSlice.sizePhase(find(AQSlice.nPhase>1,1,'first')), 1];
        else
          aspect=[AQSlice.nPhase(t(2))    /AQSlice.nPhase(t(1)), ...
                  AQSlice.sizePhase(t(2)) /AQSlice.sizePhase(t(1)), 1];
        end
        if any([isnan(aspect), isinf(aspect), aspect==0])
          set(AQSlice.(imageType{1}){3}, 'PlotBoxAspectRatioMode', 'auto');
        else
          set(AQSlice.(imageType{1}){3}, 'PlotBoxAspectRatio', aspect);
        end
        linkaxes([AQSlice.(imageType{1}){1},AQSlice.(imageType{1}){3}],'xy');
      end
    end
    if AQSlice.plotImage ~= 0
      data = plot_Image(AQSlice.plotImagehAxes{2}, data, AQSlice, 'kSpaceAndImage_2d_image');
      if AQSlice.plotPhase
        hKids = get(AQSlice.plotImagehAxes{4}, 'Children');
        hImg = findobj(hKids, 'flat', 'Tag', 'kSpaceAndImage_2d_image');
        if sum(AQSlice.nPhase==1) == 2
          % read-phase encoded
          aspect = [AQSlice.nRead    / AQSlice.nPhase(find(AQSlice.nPhase>1, 1, 'first')), ...
                    AQSlice.sizeRead / AQSlice.sizePhase(find(AQSlice.nPhase>1, 1, 'first')), 1];
          xData = data.Ticks(1).ReadZ;
          if ~isinf(AQSlice.sizeRead)
            xData = xData / AQSlice.LengthUnitScale;
          end
          t = find(AQSlice.nPhase>1, 1, 'first');
          yData = data.Ticks(t).PhaseZ;
          if ~isinf(AQSlice.sizePhase(t))
            yData = yData / AQSlice.LengthUnitScale;
          end
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
        end
        if sum(AQSlice.nPhase==1) == 2
          cData = angle(permute(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice)),[2,1]));
        else
          cData = angle(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice)));
        end
        if isempty(hImg)
          imagesc(xData, yData, cData, ...
            'Parent', AQSlice.plotImagehAxes{4}, 'Tag', 'kSpaceAndImage_2d_image');
        else
          set(hImg, 'XData', xData, 'YData', yData, 'CData', cData);
          axis(AQSlice.plotImagehAxes{4}, 'tight');
        end
        set(AQSlice.plotImagehAxes{4}, 'CLim', [-pi,pi])
        set(AQSlice.plotImagehAxes{4}, 'YDir', 'normal');
        % colorbar
        colormap(AQSlice.plotImagehAxes{4}, gray);
        if sum(AQSlice.nPhase==1) == 2
          if isinf(AQSlice.sizeRead)
            lengthUnit = 'Hz';
          else
            lengthUnit = AQSlice.LengthUnit;
          end
          xlabel(AQSlice.plotImagehAxes{2}, [AQSlice.ReadCartesianAxis{1} 'read in ' lengthUnit]);
          if isinf(AQSlice.sizePhase(t))
            % FIXME: Correct unit?
            lengthUnit = 'px';
          else
            lengthUnit = AQSlice.LengthUnit;
          end
          ylabel(AQSlice.plotImagehAxes{2}, [AQSlice.PhaseCartesianAxis{t} 'phase(' num2str(t) ') in ' lengthUnit]);
        else
          if isinf(AQSlice.sizePhase(t(2)))
            % FIXME: Correct unit?
            lengthUnit = 'px';
          else
            lengthUnit = AQSlice.LengthUnit;
          end
          xlabel(AQSlice.plotImagehAxes{2}, [AQSlice.PhaseCartesianAxis{t(2)} 'phase(' num2str(t(2)) ') in ' lengthUnit]);
          if isinf(AQSlice.sizePhase(t(2)))
            % FIXME: Correct unit?
            lengthUnit = 'px';
          else
            lengthUnit = AQSlice.LengthUnit;
          end
          ylabel(AQSlice.plotImagehAxes{2}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' lengthUnit]);
        end
        linkaxes([AQSlice.plotImagehAxes{2}, AQSlice.plotImagehAxes{4}], 'xy');
        if any([isnan(aspect), isinf(aspect), aspect==0])
          set(AQSlice.plotImagehAxes{4}, 'DataAspectRatioMode', 'auto');
        else
          set(AQSlice.plotImagehAxes{4}, 'DataAspectRatio', [1 1 1]);
        end
      end
    end

    % position axes correctly again. imagesc in one axes might affect
    % the position of neighboring axes.
    if AQSlice.plotImage ~= 0
      if AQSlice.plotkSpace ~= 0
        if AQSlice.plotPhase
          for iAxes = 1:4
            subplot(2,2,iAxes, AQSlice.plotImagehAxes{iAxes}, 'Visible', 'on');
          end
        else
          for iAxes = 1:2
            subplot(1,2,iAxes, AQSlice.plotImagehAxes{iAxes}, 'Visible', 'on');
          end
        end
      else
        if AQSlice.plotPhase
          visAxes = [2, 4];
          for iAxes = 1:2
            subplot(2,1,iAxes, AQSlice.plotImagehAxes{visAxes(iAxes)}, 'Visible', 'on');
          end
        else
          subplot(1,1,1, AQSlice.plotImagehAxes{2}, 'Visible', 'on');
        end
      end
    elseif AQSlice.plotkSpace ~= 0
      if AQSlice.plotPhase
        visAxes = [1, 3];
        for iAxes = 1:2
          subplot(2,1,iAxes, AQSlice.plotImagehAxes{visAxes(iAxes)}, 'Visible', 'on');
        end
      else
        subplot(1,1,1, AQSlice.plotImagehAxes{1}, 'Visible', 'on');
      end
    end

    if AQSlice.plotImageOs ~= 0
      hKids = get(AQSlice.plotImageOshAxes{2}, 'Children');
      hImg = findobj(hKids, 'flat', 'Tag', 'kSpaceAndImage_2d_image');
      if sum(AQSlice.nPhase==1) == 2
        % read-phase encoded
        t = find(AQSlice.nPhase>1, 1, 'first');
        aspect = [AQSlice.nRead    / AQSlice.nPhase(t), ...
                  AQSlice.sizeRead / AQSlice.sizePhase(t) 1];

        xData = data.Ticks(1).ReadOsZ;
        if ~isinf(AQSlice.sizeRead)
          xData = xData / AQSlice.LengthUnitScale;
        end
        yData = data.Ticks(t).PhaseOsZ;
        if ~isinf(AQSlice.sizePhase(t))
          yData = yData / AQSlice.LengthUnitScale;
        end
        cData = abs(permute(squeeze(data.ImageOsZ(:,:,:,:,AQSlice.iSlice)), [2,1]) / ...
          AQSlice.AmplitudeUnitScale);
        if isempty(hImg)
          imagesc(xData, yData, cData, ...
            'Parent', AQSlice.plotImageOshAxes{2}, 'Tag', 'kSpaceAndImage_2d_image');
        else
          set(hImg, 'XData', xData, 'YData', yData, 'CData', cData);
          axis(AQSlice.plotImageOshAxes{2}, 'tight');
        end
        set(AQSlice.plotImageOshAxes{2}, 'YDir', 'normal');
        colormap(AQSlice.plotImageOshAxes{2}, gray);
        title(AQSlice.plotImageOshAxes{2}, ['Image Amplitude in ' AQSlice.AmplitudeUnit]);
        if isinf(AQSlice.sizeRead)
          lengthUnit = 'Hz';
        else
          lengthUnit = AQSlice.LengthUnit;
        end
        xlabel(AQSlice.plotImageOshAxes{2}, [AQSlice.ReadCartesianAxis{1} 'read in ' lengthUnit]);
        if isinf(AQSlice.sizePhase(t))
          % FIXME: Correct unit?
          lengthUnit = 'px';
        else
          lengthUnit = AQSlice.LengthUnit;
        end
        ylabel(AQSlice.plotImageOshAxes{2}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' lengthUnit]);
      else
        % phase-phase encoded
        t = find(AQSlice.nPhase>1, 2, 'first');
        aspect = [AQSlice.nPhase(t(1))    / AQSlice.nPhase(t(2)), ...
                  AQSlice.sizePhase(t(1)) / AQSlice.sizePhase(t(2)) 1];

        xData = data.Ticks(t(2)).PhaseOsZ;
        if ~isinf(AQSlice.sizePhase(t(2)))
          xData = xData / AQSlice.LengthUnitScale;
        end
        yData = data.Ticks(t(1)).PhaseOsZ;
        if ~isinf(AQSlice.sizePhase(t(1)))
          yData = yData / AQSlice.LengthUnitScale;
        end
        cData = abs(squeeze(data.ImageOsZ) / AQSlice.AmplitudeUnitScale);
        if isempty(hImg)
          imagesc(xData, yData, cData, ...
            'Parent', AQSlice.plotImageOshAxes{2}, 'Tag', 'kSpaceAndImage_2d_image');
        else
          set(hImg, 'XData', xData, 'YData', yData, 'CData', cData);
          axis(AQSlice.plotImageOshAxes{2}, 'tight');
        end
        set(AQSlice.plotImageOshAxes{2}, 'YDir', 'normal');
        colormap(AQSlice.plotImageOshAxes{2}, gray);
        title(AQSlice.plotImageOshAxes{2}, ['Image Amplitude in ' AQSlice.AmplitudeUnit]);
        if isinf(AQSlice.sizePhase(t(2)))
          % FIXME: Correct unit?
          lengthUnit = 'px';
        else
          lengthUnit = AQSlice.LengthUnit;
        end
        xlabel(AQSlice.plotImageOshAxes{2}, [AQSlice.PhaseCartesianAxis{t(2)} 'phase(' num2str(t(2)) ') in ' lengthUnit]);
        if isinf(AQSlice.sizePhase(t(1)))
          % FIXME: Correct unit?
          lengthUnit = 'px';
        else
          lengthUnit = AQSlice.LengthUnit;
        end
        ylabel(AQSlice.plotImageOshAxes{2}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' lengthUnit]);
      end
      if any([isnan(aspect), isinf(aspect), aspect==0])
        set(AQSlice.plotImageOshAxes{2}, 'DataAspectRatioMode', 'auto');
      else
        set(AQSlice.plotImageOshAxes{2}, 'DataAspectRatio', [1 1 1]);
      end
      if AQSlice.plotPhase
        hKids = get(AQSlice.plotImageOshAxes{4}, 'Children');
        hImg = findobj(hKids, 'flat', 'Tag', 'kSpaceAndImage_2d_image');
        if sum(AQSlice.nPhase==1) == 2
          cData = angle(permute(squeeze(data.ImageOsZ(:,:,:,:,AQSlice.iSlice)), [2,1]));
        else
          cData = angle(squeeze(data.ImageOsZ(:,:,:,:,AQSlice.iSlice)));
        end
        if isempty(hImg)
          imagesc(xData, yData, cData, ...
            'Parent', AQSlice.plotImageOshAxes{4}, 'Tag', 'kSpaceAndImage_2d_image');
        else
          set(hImg, 'XData', xData, 'YData', yData, 'CData', cData);
          axis(AQSlice.plotImageOshAxes{4}, 'tight');
        end
        set(AQSlice.plotImageOshAxes{4}, 'CLim', [-pi,pi])
        set(AQSlice.plotImageOshAxes{4}, 'YDir', 'normal');
        % colorbar
        colormap(AQSlice.plotImageOshAxes{4}, gray)
        if sum(AQSlice.nPhase==1) == 2
          aspect = [AQSlice.nRead    / AQSlice.nPhase(find(AQSlice.nPhase>1, 1, 'first')), ...
                    AQSlice.sizeRead / AQSlice.sizePhase(find(AQSlice.nPhase>1, 1, 'first')) 1];
        else
          t = find(AQSlice.nPhase>1, 2, 'first');
          aspect = [AQSlice.nPhase(t(1))    / AQSlice.nPhase(t(2)), ...
                    AQSlice.sizePhase(t(1)) / AQSlice.sizePhase(t(2)) 1];
        end
        if any([isnan(aspect), isinf(aspect), aspect==0])
          set(AQSlice.plotImageOshAxes{4}, 'DataAspectRatioMode', 'auto');
        else
          set(AQSlice.plotImageOshAxes{4}, 'DataAspectRatio', [1 1 1]);
        end
        linkaxes([AQSlice.plotImageOshAxes{2}, AQSlice.plotImageOshAxes{4}], 'xy');
      end
    end

    % re-set tags
    % FIXME: Re-code to not lose them in the first place (when calling imagesc)
    set(AQSlice.plotImagehAxes{1}, 'Tag', 'kSpaceAndImage_kspace_amp_2d');
    set(AQSlice.plotImagehAxes{2}, 'Tag', 'kSpaceAndImage_image_amp_2d');
    set(AQSlice.plotImagehAxes{3}, 'Tag', 'kSpaceAndImage_kspace_phase_2d');
    set(AQSlice.plotImagehAxes{4}, 'Tag', 'kSpaceAndImage_image_phase_2d');


    % if (AQSlice.plotB0ppm ~= 0) || (AQSlice.plotB0Hz ~= 0)
    % find Region of Interest
    if isemptyfield(data, 'RoI') || ...
        (~(isscalar(data.RoI) && data.RoI == 1) && ...
        ~all(size(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice))) == size(data.RoI(:,:,min(end,AQSlice.iSlice)))))
      if sum([AQSlice.nPhase(1)>1, AQSlice.nPhase(2)>1, AQSlice.nPhase(3)>1]) == 1
        ZeroFillFactorPermuted = data.ZeroFillFactor([1, (find(AQSlice.nPhase>1)+1)]);
      else
        ZeroFillFactorPermuted = data.ZeroFillFactor([false, AQSlice.nPhase>1]);
      end
      data = smoothRoI(data, AQSlice, squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice)), ....
        2*ZeroFillFactorPermuted+1);
      data.RoI = double(data.RoI);  % There is no logical NaN.
      data.RoI(~data.RoI) = NaN;
    end
    % end

    if any([AQSlice.plotB0ppm, AQSlice.plotB0Hz, AQSlice.plotB0PpmPhase, ...
            AQSlice.plotB0Gradient, isinf(AQSlice.sizeRead), AQSlice.plotFft1_data] ~= 0)
      data.fCenter = data.f_fft1_data(floor(AQSlice.nRead*AQSlice.ReadOS/2)+1,AQSlice.UseAQWindow(1),AQSlice.UsetRep(1));
    end

    if AQSlice.plotB0ppm ~= 0
      hParent = AQSlice.plotB0ppm;
      if ~ishghandle(hParent, 'uipanel')
        if (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
          hParent = figure(hParent);
        end
        hParent = clf(hParent);
      else
        delete(get(hParent, 'Children'));
      end
      AQSlice.plotB0ppmhAxes{1} = axes('Parent', hParent);

      if sum(AQSlice.nPhase==1) == 2
        t=find(AQSlice.nPhase>1,1,'first');
        imagesc(data.Ticks(1).ReadZ/AQSlice.LengthUnitScale, ...
                data.Ticks(t).PhaseZ/AQSlice.LengthUnitScale, ...
                permute(data.RoI,[2,1]).*permute(squeeze(unwrap3Dmiddle(-angle(data.ImageZ(:,:,:,:,AQSlice.iSlice)))/pi/2/AQSlice.tEcho/data.fCenter*1e6),[2,1]), ...
                'Parent', AQSlice.plotB0ppmhAxes{1});
        set(AQSlice.plotB0ppmhAxes{1}, 'YDir', 'normal');
        cbB0ppm = colorbar('peer', AQSlice.plotB0ppmhAxes{1});
        ylabel(cbB0ppm, 'ppm');
        colormap(AQSlice.plotB0ppmhAxes{1}, jet);
        xlabel(AQSlice.plotB0ppmhAxes{1}, [AQSlice.ReadCartesianAxis{1} 'read in ' AQSlice.LengthUnit]);
        ylabel(AQSlice.plotB0ppmhAxes{1}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit]);
      else
        t = find(AQSlice.nPhase>1, 2, 'first');
        imagesc(data.Ticks(t(2)).PhaseZ/AQSlice.LengthUnitScale, ...
                data.Ticks(t(1)).PhaseZ/AQSlice.LengthUnitScale, ...
                data.RoI.*squeeze(unwrap3Dmiddle(-angle(data.ImageZ(:,:,:,:,AQSlice.iSlice)))/pi/2/AQSlice.tEcho/data.fCenter*1e6), ...
                'Parent', AQSlice.plotB0ppmhAxes{1});
        set(AQSlice.plotB0ppmhAxes{1}, 'YDir', 'normal');
        cbB0ppm = colorbar('peer', AQSlice.plotB0ppmhAxes{1});
        ylabel(cbB0ppm, 'ppm');
        colormap(AQSlice.plotB0ppmhAxes{1}, jet);
        xlabel(AQSlice.plotB0ppmhAxes{1}, [AQSlice.PhaseCartesianAxis{t(2)} 'phase(' num2str(t(2)) ') in ' AQSlice.LengthUnit]);
        ylabel(AQSlice.plotB0ppmhAxes{1}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit]);
      end
      set(AQSlice.plotB0ppmhAxes{1}, 'DataAspectRatio', [1 1 1])
      title(AQSlice.plotB0ppmhAxes{1}, ...
        ['Offset to = ' num2str(data.fCenter*1e-6, '%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3, '%6.6f') ' mT, in ppm']);
    end

    if AQSlice.plotB0Hz ~= 0
      hParent = AQSlice.plotB0Hz;
      if ~ishghandle(hParent, 'uipanel')
        if (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
          hParent = figure(hParent);
        end
        hParent = clf(hParent);
      else
        delete(get(hParent, 'Children'));
      end
      AQSlice.plotB0HzhAxes{1} = axes('Parent', hParent);

      if sum(AQSlice.nPhase==1) == 2
        t=find(AQSlice.nPhase>1,1,'first');
        imagesc(data.Ticks(1).ReadZ/AQSlice.LengthUnitScale, ...
                data.Ticks(t).PhaseZ/AQSlice.LengthUnitScale, ...
                permute(data.RoI,[2,1]).*permute(squeeze(unwrap3Dmiddle(-angle(data.ImageZ(:,:,:,:,AQSlice.iSlice)))/pi/2/AQSlice.tEcho),[2,1]), ...
                'Parent', AQSlice.plotB0HzhAxes{1});
        set(AQSlice.plotB0HzhAxes{1}, 'YDir', 'normal');
        cbB0Hz = colorbar('peer', AQSlice.plotB0HzhAxes{1});
        ylabel(cbB0Hz, 'Hz');
        colormap(AQSlice.plotB0HzhAxes{1}, jet);
        xlabel(AQSlice.plotB0HzhAxes{1}, [AQSlice.ReadCartesianAxis{1} 'read in ' AQSlice.LengthUnit]);
        ylabel(AQSlice.plotB0HzhAxes{1}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit]);
      else
        t = find(AQSlice.nPhase>1, 2, 'first');
        imagesc(data.Ticks(t(2)).PhaseZ/AQSlice.LengthUnitScale, ...
                data.Ticks(t(1)).PhaseZ/AQSlice.LengthUnitScale, ...
                data.RoI.*(squeeze(unwrap3Dmiddle(-angle(data.ImageZ(:,:,:,:,AQSlice.iSlice)))/pi/2/AQSlice.tEcho)), ...
                'Parent', AQSlice.plotB0HzhAxes{1});
        set(AQSlice.plotB0HzhAxes{1}, 'YDir', 'normal');
        cbB0Hz = colorbar('peer', AQSlice.plotB0HzhAxes{1});
        ylabel(cbB0Hz, 'Hz');
        colormap(AQSlice.plotB0HzhAxes{1}, jet);
        xlabel(AQSlice.plotB0HzhAxes{1}, [AQSlice.PhaseCartesianAxis{t(2)} 'phase(' num2str(t(2)) ') in ' AQSlice.LengthUnit]);
        ylabel(AQSlice.plotB0HzhAxes{1}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit]);
      end
      set(AQSlice.plotB0HzhAxes{1}, 'DataAspectRatio', [1 1 1])
      title(AQSlice.plotB0HzhAxes{1}, ...
        ['Offset to = ' num2str(data.fCenter*1e-6, '%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3, '%6.6f') ' mT, in Hz']);
    end

    if AQSlice.plotB0PpmPhase ~= 0 && isfield(data, 'ImageZFrequency')
      hParent = AQSlice.plotB0PpmPhase;
      if ~ishghandle(hParent, 'uipanel')
        if (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
          hParent = figure(hParent);
        end
        hParent = clf(hParent);
      else
        delete(get(hParent, 'Children'));
      end
      AQSlice.plotB0ppmhAxes{2} = axes('Parent', hParent);

      t = find(AQSlice.nPhase>1, 2, 'first');
      imagesc(data.Ticks(t(2)).PhaseZ/AQSlice.LengthUnitScale, ...
              data.Ticks(t(1)).PhaseZ/AQSlice.LengthUnitScale, ...
              data.RoI.*(squeeze(data.ImageZFrequency(:,:,:,:,AQSlice.iSlice)))./data.fCenter.*1e6, ...
              'Parent', AQSlice.plotB0ppmhAxes{2});
      set(AQSlice.plotB0ppmhAxes{2}, 'YDir', 'normal');
      cbB0ppm = colorbar('peer', AQSlice.plotB0ppmhAxes{2});
      ylabel(cbB0ppm, 'ppm');
      colormap(AQSlice.plotB0ppmhAxes{2}, jet);
      xlabel(AQSlice.plotB0ppmhAxes{2}, [AQSlice.PhaseCartesianAxis{t(2)} 'phase(' num2str(t(2)) ') in ' AQSlice.LengthUnit]);
      ylabel(AQSlice.plotB0ppmhAxes{2}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit]);

      set(AQSlice.plotB0ppmhAxes{2}, 'DataAspectRatio', [1 1 1])
      title(AQSlice.plotB0ppmhAxes{2}, ...
        ['Offset to = ' num2str(data.fCenter*1e-6, '%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3, '%6.6f') ' mT, in ppm']);
    end

    if AQSlice.plotB0HzPhase ~= 0 && isfield(data, 'ImageZFrequency')
      hParent = AQSlice.plotB0HzPhase;
      if ~ishghandle(hParent, 'uipanel')
        if (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
          hParent = figure(hParent);
        end
        hParent = clf(hParent);
      else
        delete(get(hParent, 'Children'));
      end
      AQSlice.plotB0HzhAxes{3} = axes('Parent', hParent);

      t = find(AQSlice.nPhase>1, 2, 'first');
      imagesc(data.Ticks(t(2)).PhaseZ/AQSlice.LengthUnitScale, ...
              data.Ticks(t(1)).PhaseZ/AQSlice.LengthUnitScale, ...
              data.RoI.*(squeeze(data.ImageZFrequency(:,:,:,:,AQSlice.iSlice))), ...
              'Parent', AQSlice.plotB0HzhAxes{3});
      set(AQSlice.plotB0HzhAxes{3}, 'YDir', 'normal');
      cbB0ppm = colorbar('peer', AQSlice.plotB0HzhAxes{3});
      ylabel(cbB0ppm, 'ppm');
      colormap(AQSlice.plotB0HzhAxes{3}, jet);
      xlabel(AQSlice.plotB0HzhAxes{3}, [AQSlice.PhaseCartesianAxis{t(2)} 'phase(' num2str(t(2)) ') in ' AQSlice.LengthUnit]);
      ylabel(AQSlice.plotB0HzhAxes{3}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit]);

      set(AQSlice.plotB0HzhAxes{3}, 'DataAspectRatio', [1 1 1])
      title(AQSlice.plotB0HzhAxes{3}, ...
        ['Offset to = ' num2str(data.fCenter*1e-6, '%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3, '%6.6f') ' mT, in Hz']);
    end

    if AQSlice.plotB0Gradient ~= 0
      hParent = AQSlice.plotB0Gradient;
      if ~ishghandle(hParent, 'uipanel')
        if (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
          hParent = figure(hParent);
        end
        hParent = clf(hParent);
      else
        delete(get(hParent, 'Children'));
      end
      AQSlice.plotB0GradienthAxes{1} = axes('Parent', hParent);

      B0Gradient_xyz=cell(2,1);
      if sum(AQSlice.nPhase==1) == 2
        t = find(AQSlice.nPhase>1, 1, 'first');
        B0=permute(data.RoI,[2,1]).*permute(squeeze(unwrap3Dmiddle(-angle(data.ImageZ(:,:,:,:,AQSlice.iSlice)))/pi/2/AQSlice.tEcho*2*pi/AQSlice.Gamma),[2,1]);
        [B0Gradient_xyz{[2,1]}] = gradient(B0, diff(data.Ticks(1).ReadZ(1:2)), diff(data.Ticks(t).PhaseZ(1:2)));
        B0Gradient_mean(2) = mean(reshape(B0Gradient_xyz{2}, [], 1), 'omitnan');
        B0Gradient_mean(1) = mean(reshape(B0Gradient_xyz{1}, [], 1), 'omitnan');
        B0Gradient_xyz = cat(4, B0Gradient_xyz{:});
        B0Gradient = sqrt(sum(B0Gradient_xyz.^2, 4));

        xData = data.Ticks(1).ReadZ;
        if ~isinf(AQSlice.sizeRead)
          xData = xData / AQSlice.LengthUnitScale;
        end
        yData = data.Ticks(t).PhaseZ;
        if ~isinf(AQSlice.sizePhase(t))
          yData = yData / AQSlice.LengthUnitScale;
        end
        imagesc(xData, yData, B0Gradient*1e3, ...
                'Parent', AQSlice.plotB0GradienthAxes{1});
        set(AQSlice.plotB0GradienthAxes{1}, 'YDir', 'normal');
        cbB0Gradient = colorbar('peer', AQSlice.plotB0GradienthAxes{1});
        ylabel(cbB0Gradient, 'mT/m');
        colormap(AQSlice.plotB0GradienthAxes{1}, jet);
        if isinf(AQSlice.sizeRead)
          lengthUnit = 'Hz';
        else
          lengthUnit = AQSlice.LengthUnit;
        end
        xlabel(AQSlice.plotB0GradienthAxes{1}, ...
          [AQSlice.ReadCartesianAxis{1} 'read in ' lengthUnit ' (mean ' num2str(B0Gradient_mean(2)*1e3,'%10.3f') ' mT/m)']);
        if isinf(AQSlice.sizePhase(t))
          % FIXME: Correct unit?
          lengthUnit = 'px';
        else
          lengthUnit = AQSlice.LengthUnit;
        end
        ylabel(AQSlice.plotB0GradienthAxes{1}, ...
          [AQSlice.PhaseCartesianAxis{t} 'phase(' num2str(t) ') in ' lengthUnit ' (mean ' num2str(B0Gradient_mean(1)*1e3,'%10.3f') ' mT/m)']);
      else
        t = find(AQSlice.nPhase>1, 2, 'first');
        B0 = data.RoI ...
          .* squeeze(unwrap3Dmiddle(-angle(data.ImageZ(:,:,:,:,AQSlice.iSlice)))) ...
          / AQSlice.tEcho / AQSlice.Gamma;
        [B0Gradient_xyz{[2,1]}] = gradient(B0, diff(data.Ticks(t(2)).Phase(1:2)), diff(data.Ticks(t(1)).PhaseZ(1:2)));
        B0Gradient_mean(2) = mean(reshape(B0Gradient_xyz{2}, [], 1), 'omitnan');
        B0Gradient_mean(1) = mean(reshape(B0Gradient_xyz{1}, [], 1), 'omitnan');
        B0Gradient_xyz = cat(4, B0Gradient_xyz{:});
        B0Gradient = sqrt(sum(B0Gradient_xyz.^2, 4));

        imagesc(data.Ticks(t(2)).PhaseZ/AQSlice.LengthUnitScale, ...
                data.Ticks(t(1)).PhaseZ/AQSlice.LengthUnitScale, ...
                B0Gradient*1e3, ...
                'Parent', AQSlice.plotB0GradienthAxes{1});
        set(AQSlice.plotB0GradienthAxes{1}, 'YDir', 'normal');
        cbB0Gradient = colorbar('peer', AQSlice.plotB0GradienthAxes{1});
        ylabel(cbB0Gradient, 'mT/m');
        colormap(AQSlice.plotB0GradienthAxes{1}, jet);
        xlabel(AQSlice.plotB0GradienthAxes{1}, ...
          [AQSlice.PhaseCartesianAxis{t(2)} 'phase(' num2str(t(2)) ') in ' AQSlice.LengthUnit ' (mean ' num2str(B0Gradient_mean(2)*1e3,'%10.3f') ' mT/m)']);
        ylabel(AQSlice.plotB0GradienthAxes{1}, ...
          [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit ' (mean ' num2str(B0Gradient_mean(1)*1e3,'%10.3f') ' mT/m)']);
      end
      set(AQSlice.plotB0GradienthAxes{1}, 'DataAspectRatio', [1 1 1])
      title(AQSlice.plotB0GradienthAxes{1}, ...
        ['Gradient at = ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT, in mT/m']);
    end
  end
elseif any(szImage(1:4)>1) || AQSlice.plotImageOs
  %% 1d image
  if any([AQSlice.plotB0ppm, AQSlice.plotB0Hz, AQSlice.plotB0HzPhase, ...
        AQSlice.plotB0PpmPhase, isinf(AQSlice.sizeRead), AQSlice.plotFft1_data] ~= 0)
    data.fCenter = data.f_fft1_data(floor(AQSlice.nRead*AQSlice.ReadOS/2)+1,AQSlice.UseAQWindow(1),AQSlice.UsetRep(1));
  end

  if (AQSlice.plotImage ~= 0) || (AQSlice.plotImageOs ~= 0) || (AQSlice.plotkSpace ~= 0)
    hParent{1} = AQSlice.plotImageHandle;

    % clean up parent
    if ishghandle(hParent{1}, 'uipanel')
      delete(get(AQSlice.plotImageHandle, 'Children'));
    else
      if AQSlice.raiseFigures || ~ishghandle(hParent{1}, 'figure')
        hParent{1} = figure(hParent{1});
      end
      hParent{1} = clf(hParent{1});
    end

    if (AQSlice.plotImage ~= 0) && (AQSlice.plotImageOs ~= 0)
      % open second figure
      hParent{2} = AQSlice.plotImageOs;
      if ~ishghandle(hParent{2}, 'uipanel')
        if AQSlice.raiseFigures || ~ishghandle(hParent(2), 'figure'), hParent{2} = figure(hParent{2}); end
        hParent{2} = clf(hParent{2});
      end
    end
    figName{1} = '';
    iOs = 1;
    imageTypes = {};
    if (AQSlice.plotImage ~= 0) || (AQSlice.plotImageOs ~= 0)
      if AQSlice.plotImage ~= 0
        iOs = iOs+1;
        imageTypes{end+1} = 'plotImagehAxes';
        if AQSlice.plotkSpace ~= 0
          if AQSlice.plotPhase
            figName{1} = 'Magnitude and Phase of k-Space, Magnitude and Phase of Image';
            AQSlice.plotImagehAxes{1} = subplot(2,2,1, 'Parent', hParent{1});
            AQSlice.plotImagehAxes{2} = subplot(2,2,2, 'Parent', hParent{1});
            AQSlice.plotImagehAxes{3} = subplot(2,2,3, 'Parent', hParent{1});
            AQSlice.plotImagehAxes{4} = subplot(2,2,4, 'Parent', hParent{1});
          else
            figName{1} = 'k-Space and Image';
            AQSlice.plotImagehAxes{1} = subplot(1,2,1, 'Parent', hParent{1});
            AQSlice.plotImagehAxes{2} = subplot(1,2,2, 'Parent', hParent{1});
          end
        else
          if AQSlice.plotPhase
            figName{1} = 'Magnitude and Phase of Image';
            AQSlice.plotImagehAxes{2} = subplot(2,1,1, 'Parent', hParent{1});
            AQSlice.plotImagehAxes{4} = subplot(2,1,2, 'Parent', hParent{1});
          else
            figName{1} = 'Image';
            AQSlice.plotImagehAxes{2} = axes('Parent', hParent{1});
          end
        end
      end
      if AQSlice.plotImageOs ~= 0
        imageTypes{end+1} = 'plotImageOshAxes';
        if AQSlice.plotkSpace ~= 0
          if AQSlice.plotPhase
            figName{iOs} = 'Magnitude and Phase of k-Space, Magnitude and Phase of Over-Sampled Image';
            AQSlice.plotImageOshAxes{1} = subplot(2,2,1, 'Parent', hParent{iOs});
            AQSlice.plotImageOshAxes{2} = subplot(2,2,2, 'Parent', hParent{iOs});
            AQSlice.plotImageOshAxes{3} = subplot(2,2,3, 'Parent', hParent{iOs});
            AQSlice.plotImageOshAxes{4} = subplot(2,2,4, 'Parent', hParent{iOs});
          else
            figName{iOs} = 'k-Space and Over-Sampled Image';
            AQSlice.plotImageOshAxes{1} = subplot(1,2,1, 'Parent', hParent{iOs});
            AQSlice.plotImageOshAxes{2} = subplot(1,2,2, 'Parent', hParent{iOs});
          end
        else
          if AQSlice.plotPhase
            figName{iOs} = 'Magnitude and Phase of Over-Sampled Image';
            AQSlice.plotImageOshAxes{2} = subplot(2,1,1, 'Parent', hParent{iOs});
            AQSlice.plotImageOshAxes{4} = subplot(2,1,2, 'Parent', hParent{iOs});
          else
            figName{iOs} = 'Over-Sampled Image';
            AQSlice.plotImageOshAxes{2} = axes('Parent', hParent{iOs});
          end
        end
      else
        iOs = 1;
      end
    elseif AQSlice.plotkSpace ~= 0
      imageTypes{end+1} = 'plotImagehAxes';
      if AQSlice.plotPhase
        figName{1} =  'Magnitude and Phase of k-Space';
        AQSlice.plotImagehAxes{1} = subplot(2,1,1, 'Parent', hParent{1});
        AQSlice.plotImagehAxes{3} = subplot(2,1,2, 'Parent', hParent{1});
      else
        figName{1} = 'k-Space';
        AQSlice.plotImagehAxes{1} = axes('Parent', hParent{1});
      end
    end
    for jOs = 1:iOs
      if ~ishghandle(hParent{jOs}, 'uipanel') && ~isempty(figName{jOs})
        set(hParent{jOs}, 'Name', figName{jOs});
      end
    end
    if AQSlice.plotkSpace ~= 0
      for imageType = imageTypes
        if sum(AQSlice.nPhase==1) == 3
          plot(AQSlice.(imageType{1}){1}, data.kTicks(1).ReadOs, ...
            20*log10(abs(squeeze(data.kSpaceOs(:,:,:,:,AQSlice.iSlice)) * ...
                         AQSlice.VoxelVolume / AQSlice.AreaCoil * data.Amplitude2Uin(AQSlice.UseAQWindow(1),1)/2^0.5*1e6)));
          if isinf(AQSlice.sizeRead)
            xlabel(AQSlice.(imageType{1}){1}, [AQSlice.ReadCartesianAxis{1} 'read points']);
          else
            xlabel(AQSlice.(imageType{1}){1}, [AQSlice.ReadCartesianAxis{1} 'read spatial frequency in 1/m']);
          end
        else
          t=find(AQSlice.nPhase>1,1,'first');
          plot(AQSlice.(imageType{1}){1}, data.kTicks(t).PhaseOs, ...
            20*log10(abs(squeeze(data.kSpaceOs(:,:,:,:,AQSlice.iSlice)) * ...
                         AQSlice.VoxelVolume / AQSlice.AreaCoil * data.Amplitude2Uin(AQSlice.UseAQWindow(1),1)/2^0.5*1e6)));
          if isinf(AQSlice.sizePhase(t))
            xlabel(AQSlice.(imageType{1}){1}, [AQSlice.PhaseCartesianAxis{t(1)} 'read points']);
          else
            xlabel(AQSlice.(imageType{1}){1}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') spatial frequency in 1/m']);
          end
        end
        title(AQSlice.(imageType{1}){1}, ['k-Space in dB' char(181) 'V']);
        ylabel(AQSlice.(imageType{1}){1}, ['Amplitude in dB' char(181) 'V']);
        grid(AQSlice.(imageType{1}){1}, 'on');
        set(AQSlice.(imageType{1}){1}, 'XMinorGrid', 'on');
        set(AQSlice.(imageType{1}){1}, 'YMinorGrid', 'on');
        k_ylim = get(AQSlice.(imageType{1}){1}, 'YLim');
        % FIXME: Optional value for lower limit (-40dB)?
        if k_ylim(1) < -40
          % don't show anything below -40 dB if upper limit is above -40dB
          k_ylim(1) = -40;
          if k_ylim(1) < k_ylim(2), set(AQSlice.(imageType{1}){1}, 'YLim', k_ylim); end
        elseif k_ylim(2) > -40
          % add a "line" with one singular point to force extending the y axis
          % to -40 dB
          k_xlim = get(AQSlice.(imageType{1}){1}, 'XLim');
          line(AQSlice.(imageType{1}){1}, mean(k_xlim), -40);
        end
        if AQSlice.plotPhase
          if sum(AQSlice.nPhase==1) == 3
            plot(AQSlice.(imageType{1}){3}, data.kTicks(1).ReadOs, angle(squeeze(data.kSpaceOs(:,:,:,:,AQSlice.iSlice))));
          else
            plot(AQSlice.(imageType{1}){3}, data.kTicks(find(AQSlice.nPhase>1,1,'first')).PhaseOs, angle(squeeze(data.kSpaceOs(:,:,:,:,AQSlice.iSlice))));
          end
          ylabel(AQSlice.(imageType{1}){3}, 'Phase in rad');
          ylim(AQSlice.(imageType{1}){3}, [-pi, pi]);
          grid(AQSlice.(imageType{1}){3}, 'on');
          set(AQSlice.(imageType{1}){3}, 'XMinorGrid', 'on');
          set(AQSlice.(imageType{1}){3}, 'YMinorGrid', 'on');
          linkaxes([AQSlice.(imageType{1}){3}, AQSlice.(imageType{1}){1}], 'x');
        end
      end
    end
    if AQSlice.plotImage ~= 0
      if sum(AQSlice.nPhase==1) == 3
        if ~isinf(AQSlice.sizeRead)
          plot(AQSlice.plotImagehAxes{2}, data.Ticks(1).ReadZ/AQSlice.LengthUnitScale, ...
               abs(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice))/AQSlice.AmplitudeUnitScale));
          xlabel(AQSlice.plotImagehAxes{2}, [AQSlice.ReadCartesianAxis{1} 'read in ' AQSlice.LengthUnit]);
          ylabel(AQSlice.plotImagehAxes{2}, ['Image Amplitude in ' AQSlice.AmplitudeUnit]);
        else
          plot(AQSlice.plotImagehAxes{2}, data.Ticks(1).ReadZ, ...
               abs(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice))/AQSlice.AmplitudeUnitScale));
          xlabel(AQSlice.plotImagehAxes{2}, ['offset to ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz in Hz']);
          ylabel(AQSlice.plotImagehAxes{2}, ['Image Amplitude in ' AQSlice.AmplitudeUnit]);
        end
      else
        t=find(AQSlice.nPhase>1,1,'first');
        plot(AQSlice.plotImagehAxes{2}, data.Ticks(t).PhaseZ/AQSlice.LengthUnitScale, ...
             abs(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice))/AQSlice.AmplitudeUnitScale));
        xlabel(AQSlice.plotImagehAxes{2}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit]);
        ylabel(AQSlice.plotImagehAxes{2}, ['Image Amplitude in ' AQSlice.AmplitudeUnit]);
      end
      grid(AQSlice.plotImagehAxes{2}, 'on');
      set(AQSlice.plotImagehAxes{2}, 'XMinorGrid', 'on');
      set(AQSlice.plotImagehAxes{2}, 'YMinorGrid', 'on');
      ylim(AQSlice.plotImagehAxes{2},ylim(AQSlice.plotImagehAxes{2}).*[0,1]);
      if AQSlice.plotPhase
        if sum(AQSlice.nPhase==1) == 3
          if ~isinf(AQSlice.sizeRead)
            plot(AQSlice.plotImagehAxes{4}, data.Ticks(1).ReadZ/AQSlice.LengthUnitScale, ...
              angle(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice))));
          else
            plot(AQSlice.plotImagehAxes{4}, data.Ticks(1).ReadZ, ...
              angle(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice))));
          end
        else
          plot(AQSlice.plotImagehAxes{4}, data.Ticks(find(AQSlice.nPhase>1,1,'first')).PhaseZ/AQSlice.LengthUnitScale, ...
            angle(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice))));
        end
        ylim(AQSlice.plotImagehAxes{4}, [-pi, pi]);
        ylabel(AQSlice.plotImagehAxes{4}, 'Phase in rad');
        grid(AQSlice.plotImagehAxes{4}, 'on');
        set(AQSlice.plotImagehAxes{4}, 'XMinorGrid', 'on');
        set(AQSlice.plotImagehAxes{4}, 'YMinorGrid', 'on');
        linkaxes([AQSlice.plotImagehAxes{4}, AQSlice.plotImagehAxes{2}], 'x');
      end
    end

    % find Region of Interest
    if isemptyfield(data, 'RoI') || ...
        (~(isscalar(data.RoI) && data.RoI == 1) && ...
        ~all(size(squeeze(data.ImageZ(:,:,:,:,AQSlice.iSlice))) == size(data.RoI)))
      if sum(AQSlice.nPhase==1) == 3
        ZeroFillFactorPermuted = data.ZeroFillFactor(1);
      else
        ZeroFillFactorPermuted = data.ZeroFillFactor(find(AQSlice.nPhase>1)+1);
      end
      % FIXME: Is the image always along the first dimension?
      data = smoothRoI(data, AQSlice, reshape(data.ImageZ(:,:,:,:,AQSlice.iSlice), [], 1), ...
        [2*ZeroFillFactorPermuted+1, 1]);
      data.RoI = double(data.RoI);  % There is no logical NaN.
      data.RoI(~data.RoI) = NaN;
    end

    if AQSlice.plotB0Hz ~= 0
      hParent = AQSlice.plotB0Hz;
      if ~ishghandle(hParent, 'uipanel')
        if (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
          hParent = figure(hParent);
        end
        hParent = clf(hParent);
      else
        delete(get(hParent, 'Children'));
      end
      AQSlice.plotB0HzhAxes{1} = axes('Parent', hParent);

     if sum(AQSlice.nPhase==1) == 3
        plot(AQSlice.plotB0HzhAxes{1}, data.Ticks(1).ReadZ/AQSlice.LengthUnitScale, ...
             data.RoI.*unwrap3Dmiddle(-angle(reshape(data.ImageZ(:,:,:,:,AQSlice.iSlice),[],1)))/pi/2/AQSlice.tEcho);
        xlabel(AQSlice.plotB0HzhAxes{1}, [AQSlice.ReadCartesianAxis{1} 'read in ' AQSlice.LengthUnit]);
        ylabel(AQSlice.plotB0HzhAxes{1}, 'Offset frequency in Hz');
      else
        t = find(AQSlice.nPhase>1, 1, 'first');
        plot(AQSlice.plotB0HzhAxes{1}, data.Ticks(t).PhaseZ/AQSlice.LengthUnitScale, ...
             data.RoI.*unwrap3Dmiddle(-angle(reshape(data.ImageZ(:,:,:,:,AQSlice.iSlice),[],1)))/pi/2/AQSlice.tEcho);
        xlabel(AQSlice.plotB0HzhAxes{1}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit]);
        ylabel(AQSlice.plotB0HzhAxes{1}, 'Offset frequency in Hz');
      end
      grid(AQSlice.plotB0HzhAxes{1}, 'on');
      set(AQSlice.plotB0HzhAxes{1}, 'XMinorGrid', 'on');
      set(AQSlice.plotB0HzhAxes{1}, 'YMinorGrid', 'on');
      title(AQSlice.plotB0HzhAxes{1}, ['Offset to = ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT']);
    end

    if AQSlice.plotB0ppm ~= 0
      hParent = AQSlice.plotB0ppm;
      if ~ishghandle(hParent, 'uipanel')
        if (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
          hParent = figure(hParent);
        end
        hParent = clf(hParent);
      else
        delete(get(hParent, 'Children'));
      end
      AQSlice.plotB0ppmhAxes{1} = axes('Parent', hParent);
      if sum(AQSlice.nPhase==1) == 3
        plot(AQSlice.plotB0ppmhAxes{1}, data.Ticks(1).ReadZ/AQSlice.LengthUnitScale, ...
             data.RoI.*unwrap3Dmiddle(-angle(reshape(data.ImageZ(:,:,:,:,AQSlice.iSlice),[],1)))/pi/2/AQSlice.tEcho/data.fCenter*1e6);
        xlabel(AQSlice.plotB0ppmhAxes{1}, [AQSlice.ReadCartesianAxis{1} 'read in ' AQSlice.LengthUnit]);
        ylabel(AQSlice.plotB0ppmhAxes{1}, 'Offset in ppm');
      else
        t = find(AQSlice.nPhase>1, 1, 'first');
        plot(AQSlice.plotB0ppmhAxes{1}, data.Ticks(t).PhaseZ/AQSlice.LengthUnitScale, ...
             data.RoI.*unwrap3Dmiddle(-angle(reshape(data.ImageZ(:,:,:,:,AQSlice.iSlice),[],1)))/pi/2/AQSlice.tEcho./data.fCenter.*1e6);
        xlabel(AQSlice.plotB0ppmhAxes{1}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit]);
        ylabel(AQSlice.plotB0ppmhAxes{1}, 'Offset in ppm');
      end
      grid(AQSlice.plotB0ppmhAxes{1}, 'on');
      set(AQSlice.plotB0ppmhAxes{1}, 'XMinorGrid', 'on');
      set(AQSlice.plotB0ppmhAxes{1}, 'YMinorGrid', 'on');
      title(AQSlice.plotB0ppmhAxes{1}, ['Offset to = ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT']);
    end

    if (AQSlice.plotB0HzPhase ~= 0) && isfield(data, 'ImageZFrequency')
      hParent = AQSlice.plotB0HzPhase;
      if ~ishghandle(hParent, 'uipanel')
        if (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
          hParent = figure(hParent);
        end
        hParent = clf(hParent);
      else
        delete(get(hParent, 'Children'));
      end
      AQSlice.plotB0HzhAxes{3} = axes('Parent', hParent);

      t = find(AQSlice.nPhase>1, 1, 'first');
      plot(AQSlice.plotB0HzhAxes{3}, data.Ticks(t).PhaseZ/AQSlice.LengthUnitScale, ...
           data.RoI.*(squeeze(data.ImageZFrequency(:,:,:,:,AQSlice.iSlice))));
      xlabel(AQSlice.plotB0HzhAxes{3}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit]);
      ylabel(AQSlice.plotB0HzhAxes{3}, 'Offset frequency in Hz');
      grid(AQSlice.plotB0HzhAxes{3}, 'on');
      set(AQSlice.plotB0HzhAxes{3}, 'XMinorGrid', 'on');
      set(AQSlice.plotB0HzhAxes{3}, 'YMinorGrid', 'on');
      title(AQSlice.plotB0HzhAxes{3}, ['Offset to = ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT']);
    end

    if (AQSlice.plotB0PpmPhase ~= 0) && isfield(data, 'ImageZFrequency')
      hParent = AQSlice.plotB0PpmPhase;
      if ~ishghandle(hParent, 'uipanel')
        if (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
          hParent = figure(hParent);
        end
        hParent = clf(hParent);
      else
        delete(get(hParent, 'Children'));
      end
      AQSlice.plotB0ppmhAxes{2} = axes('Parent', hParent);

      t = find(AQSlice.nPhase>1, 1, 'first');
      plot(AQSlice.plotB0ppmhAxes{2}, data.Ticks(t).PhaseZ/AQSlice.LengthUnitScale, ...
           data.RoI.*(squeeze(data.ImageZFrequency(:,:,:,:,AQSlice.iSlice)))./data.fCenter.*1e6);
      xlabel(AQSlice.plotB0ppmhAxes{2}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit]);
      ylabel(AQSlice.plotB0ppmhAxes{2}, 'Offset in ppm');
      grid(AQSlice.plotB0ppmhAxes{2}, 'on');
      set(AQSlice.plotB0ppmhAxes{2}, 'XMinorGrid', 'on');
      set(AQSlice.plotB0ppmhAxes{2}, 'YMinorGrid', 'on');
      title(AQSlice.plotB0ppmhAxes{2}, ['Offset to = ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz, ' num2str(data.fCenter*2*pi/AQSlice.Gamma*1e3,'%6.6f') ' mT']);
    end

    if AQSlice.plotImageOs ~= 0
      if sum(AQSlice.nPhase==1) == 3
        if isinf(AQSlice.sizeRead)
          lengthTicks = data.Ticks(1).ReadOsZ;
          xTickLabel = ['offset to ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz in Hz'];
        else
          lengthTicks = data.Ticks(1).ReadOsZ/AQSlice.LengthUnitScale;
          xTickLabel = [AQSlice.ReadCartesianAxis{1} 'read in ' AQSlice.LengthUnit];
        end
        plot(AQSlice.plotImageOshAxes{2}, lengthTicks, ...
             abs(squeeze(data.ImageOsZ(:,:,:,:,AQSlice.iSlice))/AQSlice.AmplitudeUnitScale));
        xlabel(AQSlice.plotImageOshAxes{2}, xTickLabel);
        ylabel(AQSlice.plotImageOshAxes{2}, ['Image Amplitude in ' AQSlice.AmplitudeUnit]);
      else
        t = find(AQSlice.nPhase>1, 1, 'first');
        lengthTicks = data.Ticks(t).PhaseOsZ / AQSlice.LengthUnitScale;
        plot(AQSlice.plotImageOshAxes{2}, lengthTicks, ...
             abs(squeeze(data.ImageOsZ(:,:,:,:,AQSlice.iSlice))/AQSlice.AmplitudeUnitScale));
        xlabel(AQSlice.plotImageOshAxes{2}, [AQSlice.PhaseCartesianAxis{t(1)} 'phase(' num2str(t(1)) ') in ' AQSlice.LengthUnit]);
        ylabel(AQSlice.plotImageOshAxes{2}, ['Image Amplitude in ' AQSlice.AmplitudeUnit]);
      end
      grid(AQSlice.plotImageOshAxes{2}, 'on');
      set(AQSlice.plotImageOshAxes{2}, 'XMinorGrid', 'on');
      set(AQSlice.plotImageOshAxes{2}, 'YMinorGrid', 'on');
      ylim(AQSlice.plotImageOshAxes{2},ylim(AQSlice.plotImageOshAxes{2}).*[0,1]);
      if AQSlice.plotPhase
        plot(AQSlice.plotImageOshAxes{4}, lengthTicks, ...
          angle(squeeze(data.ImageOsZ(:,:,:,:,AQSlice.iSlice))));
        ylim(AQSlice.plotImageOshAxes{4}, [-pi, pi]);
        ylabel(AQSlice.plotImageOshAxes{4}, 'Phase in rad');
        grid(AQSlice.plotImageOshAxes{4}, 'on');
        set(AQSlice.plotImageOshAxes{4}, 'XMinorGrid', 'on');
        set(AQSlice.plotImageOshAxes{4}, 'YMinorGrid', 'on');
        linkaxes([AQSlice.plotImageOshAxes{4}, AQSlice.plotImageOshAxes{2}], 'x');
      end
    end
  end

  if AQSlice.plotData ~= 0 || AQSlice.plotFft1_data ~= 0
    hParent = max(AQSlice.plotData, AQSlice.plotFft1_data);
    if ~ishghandle(hParent, 'uipanel')
      if (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
        hParent = figure(hParent);
      end
      hParent = clf(hParent);
    else
      delete(get(hParent, 'Children'));
    end
    if AQSlice.plotData ~= 0 && AQSlice.plotFft1_data ~= 0
      if AQSlice.plotPhase
        figName{1} = 'Magnitude and phase of data and fft1_data';
        AQSlice.plotDatahAxes{1} = subplot(2,2,1, 'Parent', hParent);
        AQSlice.plotDatahAxes{2} = subplot(2,2,2, 'Parent', hParent);
        AQSlice.plotDatahAxes{3} = subplot(2,2,3, 'Parent', hParent);
        AQSlice.plotDatahAxes{4} = subplot(2,2,4, 'Parent', hParent);
      else
        figName{1} = 'Magnitude of data and fft1_data';
        AQSlice.plotDatahAxes{1} = subplot(1,2,1, 'Parent', hParent);
        AQSlice.plotDatahAxes{2} = subplot(1,2,2, 'Parent', hParent);
      end
    elseif AQSlice.plotData ~= 0
      if AQSlice.plotPhase
        figName{1} = 'Magnitude and phase of data';
        AQSlice.plotDatahAxes{1} = subplot(2,1,1, 'Parent', hParent);
        AQSlice.plotDatahAxes{3} = subplot(2,1,2, 'Parent', hParent);
      else
        figName{1} = 'Magnitude of data';
        AQSlice.plotDatahAxes{1} = axes('Parent', hParent);
      end
    elseif AQSlice.plotFft1_data ~= 0
      if AQSlice.plotPhase
        figName{1} = 'Magnitude and phase of fft1_data';
        AQSlice.plotDatahAxes{2} = subplot(2,1,1, 'Parent', hParent);
        AQSlice.plotDatahAxes{4} = subplot(2,1,2, 'Parent', hParent);
      else
        figName{1} = 'Magnitude of fft1_data';
        AQSlice.plotDatahAxes{2} = axes('Parent', hParent);
      end
    end
    if ~ishghandle(hParent, 'uipanel') && ~isempty(figName{1})
      set(hParent, 'Name', figName{1});
    end

    if AQSlice.plotData ~= 0
      d=20*log10(abs(squeeze(data.kSpaceOsRaw(:,:,:,:,AQSlice.iSlice)) ...
        * data.Amplitude2Uin(AQSlice.UseAQWindow(1),1)/2^0.5*1e6));
      if isinf(AQSlice.sizeRead)
        hLinesCurr = plot(AQSlice.plotDatahAxes{1}, data.kTicks(1).ReadOs*diff(data.time_of_tRep([1,2],AQSlice.UseAQWindow(1),AQSlice.UsetRep(1)))*1e3, d);
        xlabel(AQSlice.plotDatahAxes{1}, ['Time offset to ' num2str(data.tImageZ*1e3,'%10.3f') ' ms, in ms']);
      else
        hLinesCurr = plot(AQSlice.plotDatahAxes{1}, data.kTicks(1).ReadOs, d);
        xlabel(AQSlice.plotDatahAxes{1}, [AQSlice.ReadCartesianAxis{1} 'read spatial frequency in 1/m']);
      end
      set(hLinesCurr, {'Color'}, mat2cell(parula(size(d, 2)), ones(1, size(d, 2)), 3));
      title(AQSlice.plotDatahAxes{1}, ['k-Space in dB' char(181) 'V']);
      ylabel(AQSlice.plotDatahAxes{1}, ['Amplitude in dB' char(181) 'V']);
      grid(AQSlice.plotDatahAxes{1}, 'on');
      set(AQSlice.plotDatahAxes{1}, 'XMinorGrid', 'on');
      set(AQSlice.plotDatahAxes{1}, 'YMinorGrid', 'on');
      k_ylim = [-40,0] + [0,1].*get(AQSlice.plotDatahAxes{1}, 'YLim');
      if k_ylim(1) < k_ylim(2)
        set(AQSlice.plotDatahAxes{1}, 'YLim', k_ylim);
      end
      xlim(AQSlice.plotDatahAxes{1}, [-Inf, Inf]);
      if AQSlice.plotPhase
      if isinf(AQSlice.sizeRead)
         hLinesCurr = plot(AQSlice.plotDatahAxes{3}, data.kTicks(1).ReadOs*diff(data.time_of_tRep([1,2],AQSlice.UseAQWindow(1),AQSlice.UsetRep(1)))*1e3, ...
           angle(squeeze(data.kSpaceOsRaw(:,:,:,:,AQSlice.iSlice))));
      else
         hLinesCurr = plot(AQSlice.plotDatahAxes{3}, data.kTicks(1).ReadOs, angle(squeeze(data.kSpaceOsRaw(:,:,:,:,AQSlice.iSlice))));
      end
        ylabel(AQSlice.plotDatahAxes{3}, 'Phase in rad');
        ylim(AQSlice.plotDatahAxes{3}, [-pi, pi]);
        grid(AQSlice.plotDatahAxes{3}, 'on');
        set(AQSlice.plotDatahAxes{3}, 'XMinorGrid', 'on');
        set(AQSlice.plotDatahAxes{3}, 'YMinorGrid', 'on');
        set(hLinesCurr, {'Color'}, mat2cell(parula(size(d, 2)), ones(1, size(d, 2)), 3));
        linkaxes([AQSlice.plotDatahAxes{1}, AQSlice.plotDatahAxes{3}], 'x');
      end

    end
    if AQSlice.plotFft1_data ~= 0
      d = abs(squeeze(data.fft1_dataCut(:,:,:,:,AQSlice.iSlice)) ...
              * data.Amplitude2Uin(AQSlice.UseAQWindow(1),1) / 2^0.5 * 1e6);
      if isinf(AQSlice.sizeRead)
        hLinesCurr=plot(AQSlice.plotDatahAxes{2}, squeeze(data.f_fft1_dataCut(:,:,:,:,AQSlice.iSlice))-data.fCenter, d);
        ylabel(AQSlice.plotDatahAxes{2}, ['Image Amplitude in ' char(181) 'V']);
        xlabel(AQSlice.plotDatahAxes{2}, [AQSlice.ReadCartesianAxis{1} 'read offset to ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz, in Hz']);
      else
        hLinesCurr=plot(AQSlice.plotDatahAxes{2}, data.Ticks(1).ReadOs/AQSlice.LengthUnitScale, d);
        ylabel(AQSlice.plotDatahAxes{2}, ['Image Amplitude in ' char(181) 'V']);
        xlabel(AQSlice.plotDatahAxes{2}, [AQSlice.ReadCartesianAxis{1} 'read in ' AQSlice.LengthUnit]);
      end
      set(hLinesCurr, {'Color'}, mat2cell(parula(size(d, 2)), ones(1, size(d, 2)), 3));
      grid(AQSlice.plotDatahAxes{2}, 'on');
      set(AQSlice.plotDatahAxes{2}, 'XMinorGrid', 'on');
      set(AQSlice.plotDatahAxes{2}, 'YMinorGrid', 'on');
      xlim(AQSlice.plotDatahAxes{2}, [-Inf, Inf]);
      ylim(AQSlice.plotDatahAxes{2}, ylim(AQSlice.plotDatahAxes{2}) .* [0, 1]);
      if AQSlice.plotPhase
        if isinf(AQSlice.sizeRead)
          hLinesCurr=plot(AQSlice.plotDatahAxes{4}, squeeze(data.f_fft1_dataCut(:,:,:,:,AQSlice.iSlice))-data.fCenter, ...
               angle(squeeze(data.fft1_dataCut(:,:,:,:,AQSlice.iSlice)) ...
               * data.Amplitude2Uin(AQSlice.UseAQWindow(1),1)/2^0.5*1e6));
          ylabel(AQSlice.plotDatahAxes{2}, ['Image Amplitude in ' char(181) 'V']);
          xlabel(AQSlice.plotDatahAxes{2}, [AQSlice.ReadCartesianAxis{1} 'read offset to ' num2str(data.fCenter*1e-6,'%10.6f') ' MHz, in Hz']);
        else
          hLinesCurr=plot(AQSlice.plotDatahAxes{4}, data.Ticks(1).ReadOs/AQSlice.LengthUnitScale, ...
               angle(squeeze(data.fft1_dataCut(:,:,:,:,AQSlice.iSlice)) ...
               * data.Amplitude2Uin(AQSlice.UseAQWindow(1),1)/2^0.5*1e6));
          ylabel(AQSlice.plotDatahAxes{4}, 'Phase in rad');
          xlabel(AQSlice.plotDatahAxes{4}, [AQSlice.ReadCartesianAxis{1} 'read in ' AQSlice.LengthUnit]);
        end
        set(hLinesCurr, {'Color'}, mat2cell(parula(size(d, 2)), ones(1, size(d, 2)), 3));
        ylim(AQSlice.plotDatahAxes{4}, [-pi, pi]);
        ylabel(AQSlice.plotDatahAxes{4}, 'Phase in rad');
        grid(AQSlice.plotDatahAxes{4}, 'on');
        set(AQSlice.plotDatahAxes{4}, 'XMinorGrid', 'on');
        set(AQSlice.plotDatahAxes{4}, 'YMinorGrid', 'on');
        linkaxes([AQSlice.plotDatahAxes{2}, AQSlice.plotDatahAxes{4}], 'x');
      end
    end
  end
end

drawnow();

end


function AQSlice = setPhaseFigure(AQSlice, absField, phaseField, figureDefault)
%% Function that sets a default graphics value for the phase figure

if isemptyfield(AQSlice, phaseField)
  if AQSlice.plotPhase && ~isempty(AQSlice.(absField)) && AQSlice.(absField) ~= 0
    if ishghandle(AQSlice.(absField), 'figure') || ...
        (isa(AQSlice.(absField), 'double') && mod(AQSlice.(absField), 1) == 0)
      AQSlice.(phaseField) = double(AQSlice.(absField)) + 1;
    else
      AQSlice.(phaseField) = figureDefault;
    end
  else
    AQSlice.(phaseField) = 0;
  end
end

end


function hsl = plot3Ddata(hParent, data, AQSlice, imageData, isImage, readField, phaseField)
%% Plot 3d data using sliceomatic and set up the correct descriptions

if ~ishghandle(hParent, 'uipanel') && ...
    (AQSlice.raiseFigures || ~ishghandle(hParent, 'figure'))
  hParent = figure(hParent);
end

hsl = hParent;
if isappdata(hParent, 'sliceomatic')
  d = getappdata(hParent, 'sliceomatic');
  if isvalid(d.axmain_obj)
    hsl = d.axmain_obj;
  end
end

autoAspectRatio = false;
if isImage
  tickType = 'Ticks';
  if sum(AQSlice.nPhase(:)>1) == 3
    % CSI image (FIXME: This does not handle 4d images completely)
    ticks{3} = data.(tickType)(3).(phaseField);
    if isinf(AQSlice.sizePhase(3))
      autoAspectRatio = true;
      lengthUnit = 'Hz';
    else
      lengthUnit = AQSlice.LengthUnit;
      ticks{3} = ticks{3} / AQSlice.LengthUnitScale;
    end
    xLabelStr = [AQSlice.PhaseCartesianAxis{3} 'phase(3) in ' lengthUnit];
    xTitleStr = AQSlice.PhaseCartesianAxis{3};
  else
    ticks{3} = data.(tickType)(1).(readField);
    if isinf(AQSlice.sizeRead)
      autoAspectRatio = true;
      lengthUnit = 'Hz';
    else
      lengthUnit = AQSlice.LengthUnit;
      ticks{3} = ticks{3} / AQSlice.LengthUnitScale;
    end
    xLabelStr = [AQSlice.ReadCartesianAxis{1} 'read in ' lengthUnit];
    xTitleStr = AQSlice.ReadCartesianAxis{1};
  end
  ticks{1} = data.(tickType)(1).(phaseField);
  if isinf(AQSlice.sizePhase(1))
    autoAspectRatio = true;
    lengthUnit = 'Hz';
  else
    lengthUnit = AQSlice.LengthUnit;
    ticks{1} = ticks{1} / AQSlice.LengthUnitScale;
  end
  yLabelStr = [AQSlice.PhaseCartesianAxis{1} 'phase(1) in ' lengthUnit];
  ticks{2} = data.(tickType)(2).(phaseField);
  if isinf(AQSlice.sizePhase(2))
    autoAspectRatio = true;
    lengthUnit = 'Hz';
  else
    lengthUnit = AQSlice.LengthUnit;
    ticks{2} = ticks{2} / AQSlice.LengthUnitScale;
  end
  zLabelStr = [AQSlice.PhaseCartesianAxis{2} 'phase(2) in ' lengthUnit];
else  % k-space
  tickType = 'kTicks';
  if sum(AQSlice.nPhase(:)>1) == 3
    % CSI image (FIXME: This does not handle 4d images completely)
    ticks{3} = data.(tickType)(3).(phaseField);
    if isinf(AQSlice.sizePhase(3))
      autoAspectRatio = true;
      lengthUnit = 'samples';
    else
      lengthUnit = '1/m';
    end
    xLabelStr = [AQSlice.PhaseCartesianAxis{3}, 'phase(3) spatial frequency in ', lengthUnit];
    xTitleStr = AQSlice.PhaseCartesianAxis{3};
  else
    ticks{3} = data.(tickType)(1).(readField);
    if isinf(AQSlice.sizeRead)
      autoAspectRatio = true;
      lengthUnit = 'samples';
    else
      lengthUnit = '1/m';
    end
    xLabelStr = [AQSlice.ReadCartesianAxis{1}, 'read spatial frequency in ', lengthUnit];
    xTitleStr = AQSlice.ReadCartesianAxis{1};
  end
  ticks{1} = data.(tickType)(1).(phaseField);
  if isinf(AQSlice.sizePhase(1))
    autoAspectRatio = true;
    lengthUnit = 'samples';
  else
    lengthUnit = '1/m';
  end
  yLabelStr = [AQSlice.PhaseCartesianAxis{1}, 'phase(1) spatial frequency in ', lengthUnit];
  ticks{2} = data.(tickType)(2).(phaseField);
  if isinf(AQSlice.sizePhase(2))
    autoAspectRatio = true;
    lengthUnit = 'samples';
  else
    lengthUnit = '1/m';
  end
  zLabelStr = [AQSlice.PhaseCartesianAxis{2} 'phase(2) spatial frequency in ' lengthUnit];
end

hsl = sliceomatic(hsl, imageData, ...
  ticks{3}, ticks{1}, ticks{2}, AQSlice.sliceomaticProps);
light('Parent', hsl.hAxes);
xlabel(hsl.hAxes, xLabelStr);
ylabel(hsl.hAxes, yLabelStr);
zlabel(hsl.hAxes, zLabelStr);
if autoAspectRatio
  % FIXME: Is there a better way to have the correct aspect ratio between axes
  % that are encoded but a "free" aspect ratio for axes that are not encoded?
  set(hsl.hAxes, 'DataAspectRatioMode', 'auto');
  set(hsl.hAxes, 'PlotboxAspectRatioMode', 'auto');
end

title(hsl.GetSliderX(), xTitleStr);
title(hsl.GetSliderY(), AQSlice.PhaseCartesianAxis{1});
title(hsl.GetSliderZ(), AQSlice.PhaseCartesianAxis{2});

labelStrs = {xTitleStr, AQSlice.PhaseCartesianAxis{1}, AQSlice.PhaseCartesianAxis{2}};
if ~any(cellfun(@isempty, labelStrs))
  %% override sliceomatic's data cursor update function
  hFigure = ancestor(hParent, 'figure');
  DataCursorUpdateFcns = getappdata(hFigure, 'DataCursorUpdateFcns');
  % keep the original data cursor update function around as a work horse
  sliceomaticDatacursorUpdateFcn = DataCursorUpdateFcns.sliceomatic;
  DataCursorUpdateFcns.sliceomatic = @(pointDataTip, eventData) ...
    kSpaceAndImage_3d_DataCursorFcn(pointDataTip, eventData, ...
    sliceomaticDatacursorUpdateFcn, labelStrs);
  setappdata(hFigure, 'DataCursorUpdateFcns', DataCursorUpdateFcns);
  hdcm = datacursormode(hFigure);
  set(hdcm, 'UpdateFcn', @DataCursorUpdateFcnHandler);

  set(hParent, 'DeleteFcn', @kSpaceAndImage_DeleteFcn);

  %% Hijack sliceomatic's "orientation" popup menu and install our own callback
  d = getappdata(hsl.hParent, 'sliceomatic');
  set(d.orientationdropdown, 'Callback', @(h,e) rotate3dView(h, e, labelStrs, hsl.hAxes));
end

end


function outputTxt = kSpaceAndImage_3d_DataCursorFcn(pointDataTip, eventData, ...
  sliceomaticDatacursorUpdateFcn, labelStrs)
%% Data cursor update function that assigns the correct axes labels

% let sliceomatic's data cursor update function get the correct data
outputTxt = sliceomaticDatacursorUpdateFcn(pointDataTip, eventData);

if isempty(outputTxt)
  return;
end

% parse output of form 'X: %.3f\nY: %.3f\nZ: %.3f\nV: %.3f'
output = sscanf(outputTxt, '%*c: %f\n');

% sort X-Y-Z
[~, idxSort] = sort(cellfun(@(x) x(numel(x)-1), labelStrs, ...
  'UniformOutput', false));
labelStrs = labelStrs(idxSort);

% create output string
labelStrs{4} = 'V';
data = [labelStrs; num2cell(reshape(output([idxSort, 4]), 1, []))];
outputTxt = sprintf('%s: %.3f\n', data{:,1:3});
outputTxt = [outputTxt, sprintf('%s: %.3g\n', data{:,4})];

end


function kSpaceAndImage_DeleteFcn(hParent, eventData)
%% Executes on deletion of parent (axes).
%
%       kSpaceAndImage_DeleteFcn(hParent, eventData)
%
% This function un-registers the DataCursorFcn when the (sliceomatic) axes is
% deleted.

% FIXME: handle multiple sliceomatic in one figure
hFigure = ancestor(hParent, 'figure');
DataCursorUpdateFcns = getappdata(hFigure, 'DataCursorUpdateFcns');
if isfield(DataCursorUpdateFcns, 'kSpaceAndImage')
  DataCursorUpdateFcns = rmfield(DataCursorUpdateFcns, 'kSpaceAndImage');
  setappdata(hFigure, 'DataCursorUpdateFcns', DataCursorUpdateFcns);
end

end


function rotate3dView(hObject, eventData, labelStrs, hAxes)
%% Callback for the Orientation popupmenu
str = get(hObject, 'String');
str = str{get(hObject, 'Value')};

if strcmp(str, 'default')
  set(hAxes, 'View', [-37.5, 30]);
  return;
end

axesChar = cellfun(@(x) x(numel(x)-1), labelStrs, ...
  'UniformOutput', false);
k = cellfun(@(x) ~isempty(strfind(str, x)), axesChar);

if k(1) && k(2)
  set(hAxes, 'View', [0, 90]);
elseif k(1) && k(3)
  set(hAxes, 'View', [0, 0]);
elseif k(2) && k(3)
  set(hAxes, 'View', [90, 0]);
end

end
