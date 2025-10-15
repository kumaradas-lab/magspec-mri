function SeqLoop = Image3D2VolView(SeqLoop, FileName, ZeroFillFactor, ZeroFillWindowSize)
%% Store 3D data to VolView file
%
%   SeqLoop = Image3D2VolView(SeqLoop, FileName, ZeroFillFactor, ZeroFillWindowSize)
%   SeqLoop = Image3D2VolView(Image3D, FileName, ZeroFillFactor, ZeroFillWindowSize)
%
% This function converts the measurement data in "SeqLoop" or the 3D image
% "Image3D" to a format readable by VolView and stores it in the file
% "FileName". Additionally, a .vvi file is created with the same file name that
% can be used to load the data in VolView.
%
% INPUT:
%   SeqLoop:            Structure with results of a measurement.
%   Image3D:            3D matrix with image data.
%   FileName:           Optional. File name for the VolView binary data.
%                       (Default: 'output/VolView3D_yyyymmdd_HHMMSSFFF')
%   ZeroFillFactor:     Optional. Scalar or 1x3 zero fill factor for padding of the
%                       k-space before the FFT. (Default: 1)
%   ZeroFillWindowSize: Optional. To reduce the "jitter" in the image, k-space
%                       frequencies further away from the center frequency can
%                       be damped. A low value for "ZeroFillWindowSize" leads to
%                       higher damping, a "ZeroFillWindowSize" of Inf means no
%                       damping. A good value in most circumstances is 1.4
%                       (approx. sqrt(2)).
%                       (Default: 1.4).
%
% OUTPUT:
%   SeqLoop:            The same as input "SeqLoop". Or: A structure with
%                       "Image3D" and the values that where written to the
%                       VolView files.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% check input
if ~isstruct(SeqLoop)
  SeqLoop.data.Image = SeqLoop;
  SeqLoop.data.ImageOs = SeqLoop.data.Image;
  SeqLoop.AQSlice(1).ReadOS = 1;
  SeqLoop.AQSlice(1).PhaseOS = [1, 1, 1];
  SeqLoop.AQSlice(1).sizeRead = 1;
  SeqLoop.AQSlice(1).nRead = size(SeqLoop.data.Image, 1);
  SeqLoop.AQSlice(1).sizePhase = [1, 1, 1];
  SeqLoop.AQSlice(1).nPhase = [size(SeqLoop.data.Image,2), size(SeqLoop.data.Image,3), 1];
end

if isemptyfield(SeqLoop.AQSlice(1), 'VoxelVolume'),         SeqLoop.AQSlice(1).VoxelVolume = 1;         end
if isemptyfield(SeqLoop.AQSlice(1), 'AreaCoil'),            SeqLoop.AQSlice(1).AreaCoil = 1;            end
if isemptyfield(SeqLoop.AQSlice(1), 'AmplitudeUnit'),       SeqLoop.AQSlice(1).AmplitudeUnit = 'au';    end
if isemptyfield(SeqLoop.AQSlice(1), 'AmplitudeUnitScale'),  SeqLoop.AQSlice(1).AmplitudeUnitScale = 1;  end


if nargin<2 || isempty(FileName)
  FileName = ['output/VolView3D_ ' datestr(now, 'yyyymmdd_HHMMSSFFF')];
end

if nargin < 3
  if isemptyfield(SeqLoop.AQSlice(1), 'ZeroFillFactor'), SeqLoop.AQSlice(1).ZeroFillFactor = 1; end
  ZeroFillFactor = SeqLoop.AQSlice(1).ZeroFillFactor;
end
ZeroFillFactor=reshape(ZeroFillFactor(1:min(end,3)),1,[]);
if nargin < 4
  if isemptyfield(SeqLoop.AQSlice(1), 'ZeroFillWindowSize'), SeqLoop.AQSlice(1).ZeroFillWindowSize = 1.4; end
  ZeroFillWindowSize = SeqLoop.AQSlice(1).ZeroFillWindowSize;
end
ZeroFillWindowSize=reshape(ZeroFillWindowSize(1:min(end,3)),1,[]);

%% collect data
SeqLoop.AQSlice(1).resolutionRead = SeqLoop.AQSlice(1).sizeRead / SeqLoop.AQSlice(1).nRead;
SeqLoop.AQSlice(1).resolutionPhase = SeqLoop.AQSlice(1).sizePhase ./ SeqLoop.AQSlice(1).nPhase;

sImage = [SeqLoop.AQSlice(1).nRead, SeqLoop.AQSlice(1).nPhase(1), SeqLoop.AQSlice(1).nPhase(2)];
sImageZ = sImage .* ZeroFillFactor;
sImageZ(sImage==1) = 1;

lengthScaling = 1e-3;
lengthScalingUnit = 'm';  % '', 'm', Greek character mu, 'n' (note: consider using char(181) for micro)
centerOfImage = SeqLoop.AQSlice(1).Center2OriginImage(1:3) / lengthScaling;

resolutionImage = [SeqLoop.AQSlice(1).resolutionRead(1), ...
  SeqLoop.AQSlice(1).resolutionPhase(1:2)] ./ lengthScaling ./ ZeroFillFactor;
resolutionImage(isinf(resolutionImage)) = 1;

amplitudeScaling = SeqLoop.AQSlice(1).AmplitudeUnitScale;
amplitudeScalingUnit = SeqLoop.AQSlice(1).AmplitudeUnit;
if isfield(SeqLoop.data, 'ImageZ') && nargin <= 2
  % take ImageZ as it currently is
  ImageZ = SeqLoop.data.ImageZ;
elseif isfield(SeqLoop.data, 'ImageOs')
  % calculate ImageZ with new zero filling values
  szImageOs = arrayfun(@(n) size(SeqLoop.data.ImageOs, n), 1:4);
  szImageOsZ = szImageOs .* ZeroFillFactor;
  szImageOsZ(szImageOs==1) = 1;

  ImageOsZ = zeroFill_image(...
    SeqLoop.data.ImageOs / SeqLoop.AQSlice(1).VoxelVolume * ...
    SeqLoop.AQSlice(1).AreaCoil / SeqLoop.AQSlice(1).AmplitudeUnitScale, ...
    szImageOsZ, ZeroFillWindowSize);

  % cut oversampled data (ReadOS, PhaseOS(1) and PhaseOS(2) )
  MySizeOS = [ ...
    SeqLoop.AQSlice(1).nRead*SeqLoop.AQSlice(1).ReadOS, ...
    SeqLoop.AQSlice(1).nPhase(1)*SeqLoop.AQSlice(1).PhaseOS(1), ...
    SeqLoop.AQSlice(1).nPhase(2)*SeqLoop.AQSlice(1).PhaseOS(2) ] .* ZeroFillFactor;
  MySize = [ ...
    SeqLoop.AQSlice(1).nRead, ...
    SeqLoop.AQSlice(1).nPhase(1), ...
    SeqLoop.AQSlice(1).nPhase(2) ] .* ZeroFillFactor;
  ImageZ = ImageOsZ(floor(MySizeOS(1)/2)+  (1-floor(MySize(1)/2):ceil(MySize(1)/2)),...
                    floor(MySizeOS(2)/2)+  (1-floor(MySize(2)/2):ceil(MySize(2)/2)),...
                    floor(MySizeOS(3)/2)+  (1-floor(MySize(3)/2):ceil(MySize(3)/2)),:,:);
elseif nargin > 2
  error('PD:Image3D2VolView:ZeroFilling', 'Zero filling not supported on reduced data set.')
else
  error('PD:Image3D2VolView:ImageData', 'Data must contain either "ImageOs" or "ImageZ"');
end

nImages = size(ImageZ, 5);
for iImage = 1:size(ImageZ, 5)
  if nImages > 1
    imageSuffix = sprintf('_%02d', iImage);
  else
    imageSuffix = '';
  end

  %% write file with VolView Information
  fid = fopen([FileName, imageSuffix, '.vvi'], 'w');

  fprintf(fid, '<KWOpenFileProperties Version="1.5"\n');
  fprintf(fid, '                      ClassName="vtkKWOpenFileProperties"\n');
  fprintf(fid, '                      Spacing="%16.16e %16.16e %16.16e"\n', resolutionImage(:));
  fprintf(fid, '                      Origin="%16.16e %16.16e %16.16e"\n', centerOfImage(:));
  fprintf(fid, '                      DistanceUnits="%sm"\n', lengthScalingUnit);
  fprintf(fid, '                      ScalarUnits0="%s"\n', amplitudeScalingUnit);
  fprintf(fid, '                      ScalarType="11"\n');
  fprintf(fid, '                      WholeExtent="0 %u 0 %u 0 %u"\n', sImageZ-1);
  fprintf(fid, '                      NumberOfScalarComponents="1"\n');
  fprintf(fid, '                      IndependentComponents="1"\n');
  fprintf(fid, '                      FileOrientation="4 2 0"\n');
  fprintf(fid, '                      BigEndianFlag="0"\n');
  fprintf(fid, '                      FileDimensionality="3"\n');
  fprintf(fid, '                      Scope="2"/>\n');

  fclose(fid);

  %% write file with data
  fid = fopen([FileName, imageSuffix], 'w', 'l');
  fwrite(fid, abs(ImageZ(:,:,:,:,iImage))./amplitudeScaling, 'double');
  fclose(fid);
end

end
