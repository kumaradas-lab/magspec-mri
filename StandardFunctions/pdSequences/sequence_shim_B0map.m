function [SeqOut, mySave] = sequence_shim_B0map(HW, Seq, varargin)
%% Determine magnet shim by measuring the B0 map
%
%   [SeqOut, mySave] = sequence_shim_B0map(HW, Seq, AQ, TX, Grad, mySave)
%
% The function sequence_Flash is used for the actual measurement of the B0 map.
% See the documentation of that function for further settings.
% If successful, this function adds a new line to the file in HW.MagnetShimPath.
%
% INPUT:
%
%   HW
%           HW object or structure.
%
%   Seq
%           Structure with settings for the measurement. See sequence_Flash for
%           more details. Additionally, the following fields are used. Default
%           values are used if they are omitted or empty:
%
%     plotB0map
%             Display (the part of) the B0 map that was used to determine the
%             new magnet shim values. (Default: 0)
%
%     Shim
%             Structure with settings specific to this function. Default values
%             are used if the following fields are omitted or empty:
%
%       y_fov_max
%               Maximum absolute value of the y-coordinate in the magnet
%               coordinate system that is included in the B0 map for shimming.
%               (Default: 6.5e-3)
%
%       r_fov_max
%               Maximum radius (in x-z-direction of the magnet coordinate
%               system) that is included in the B0 map for shimming.
%               (Default: 3e-3)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2020-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% default input

if isemptyfield(Seq, 'Loops'), Seq.Loops = 0; end  % number of loop averages 1...
if isemptyfield(Seq, 'LoopsBreak'), Seq.LoopsBreak = 2.5; end  % Pause between two loop averages in seconds ([]= fast as possible)

if isemptyfield(Seq, 'T1'), Seq.T1 = 100e-3; end  % T1 of sample; excitation angle is acos(exp(-Seq.tRep/Seq.T1))/pi*180
if isemptyfield(Seq, 'tEcho'), Seq.tEcho = 3e-3; end  % echo time in seconds e.g. 4e-3
if isemptyfield(Seq, 'RepetitionTime'), Seq.RepetitionTime = 22e-3; end  % repetition time in seconds (default is Seq.tEcho*2)

% B0 map settings
if isemptyfield(Seq, 'CorrectB0Read'), Seq.CorrectB0Read = struct(); end
if isemptyfield(Seq.CorrectB0Read, 'Use'), Seq.CorrectB0Read.Use = true; end  % correct read offset with B0 map
if isemptyfield(Seq.CorrectB0Read, 'Get'), Seq.CorrectB0Read.Get = true; end  % get data for B0 correction
if isemptyfield(Seq.CorrectB0Read, 'Plot'), Seq.CorrectB0Read.Plot = true; end  % plot results for both acquired images
if isemptyfield(Seq.CorrectB0Read, 'tEchoIncr'), Seq.CorrectB0Read.tEchoIncr = 0.5e-3; end  % echo time increment for second measurement in s
if isemptyfield(Seq.CorrectB0Read, 'MinRelAmp'), Seq.CorrectB0Read.MinRelAmp = 0.1; end
if isemptyfield(Seq.CorrectB0Read, 'MaxRelAmpDiff'), Seq.CorrectB0Read.MaxRelAmpDiff = 2; end
if isemptyfield(Seq.CorrectB0Read, 'MaxFreqOffset'), Seq.CorrectB0Read.MaxFreqOffset = 2000; end
if isemptyfield(Seq.CorrectB0Read, 'RoIExtension'), Seq.CorrectB0Read.RoIExtension = 0; end
if isemptyfield(Seq.CorrectB0Read, 'ZeroFillWindowSize'), Seq.CorrectB0Read.ZeroFillWindowSize = 1; end
if isemptyfield(Seq, 'plotB0map'), Seq.plotB0map = 0; end  % plot B0 map

% Shim settings
if isemptyfield(Seq, 'Shim'), Seq.Shim = struct(); end
if isemptyfield(Seq.Shim, 'y_fov_max'), Seq.Shim.y_fov_max = 6.5e-3; end  % maximum absolute y in FOV used for shim (for solenoid: 5e-3, for loop-gap: 6.5e-3)
if isemptyfield(Seq.Shim, 'r_fov_max'), Seq.Shim.r_fov_max = 3e-3; end  % maximum radius of FOV for shim

% Pixels and size
if isemptyfield(Seq, 'AQSlice'), Seq.AQSlice = struct(); end
if isemptyfield(Seq.AQSlice(1), 'nRead'), Seq.AQSlice(1).nRead = 25; end  % number of pixels in read direction
if isemptyfield(Seq.AQSlice(1), 'nPhase'), Seq.AQSlice(1).nPhase = [23, 24]; end  % number of pixels in phase(1:2) direction
if isemptyfield(Seq.AQSlice(1), 'HzPerPixMin'), Seq.AQSlice(1).HzPerPixMin = 1000; end  % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ window, 0: longest possible)
if isemptyfield(Seq.AQSlice(1), 'sizeRead'), Seq.AQSlice(1).sizeRead = 14e-3; end  % size in read direction in meter
if isemptyfield(Seq.AQSlice(1), 'sizePhase'), Seq.AQSlice(1).sizePhase = ones(1,2)*14e-3; end  % size in phase(1:2) direction in meter
if isemptyfield(Seq.AQSlice(1), 'PhaseOS'), Seq.AQSlice(1).PhaseOS = [2, 2]; end  % oversampling factor in phase(1:2)

% Orientation in space
if isemptyfield(Seq.AQSlice(1), 'alfa'), Seq.AQSlice(1).alfa = 0.0*pi; end  % 1st rotation around x axis in RAD
if isemptyfield(Seq.AQSlice(1), 'phi'), Seq.AQSlice(1).phi  = 0.0*pi; end  % 2nd rotation around y axis in RAD
if isemptyfield(Seq.AQSlice(1), 'theta'), Seq.AQSlice(1).theta= 0.0*pi; end  % 3rd rotation around z axis in RAD

if isemptyfield(Seq.AQSlice(1), 'excitationPulse'), Seq.AQSlice(1).excitationPulse = @Pulse_RaisedCos; end  % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
if isemptyfield(Seq.AQSlice(1), 'SpoilLengthFactor'), Seq.AQSlice(1).SpoilLengthFactor = 3; end

% Plot
if isemptyfield(Seq, 'plotSeqAQ'), Seq.plotSeqAQ = 1:3; end  % plot sequence on real timeline, plots RF, AQ and Grad (1==x, 2==y, 3==z, 0 no gradient)
if isemptyfield(Seq, 'LoopPlot'), Seq.LoopPlot = 1; end  % plot result at every loop
if isemptyfield(Seq, 'LoopPlotAverages'), Seq.LoopPlotAverages = 1; end  % plot average at every loop
if isemptyfield(Seq, 'LoopSeqPlot'), Seq.LoopSeqPlot = 1; end
if isemptyfield(Seq.AQSlice(1), 'plotkSpace'), Seq.AQSlice(1).plotkSpace = 0; end  % plot k-space
if isemptyfield(Seq.AQSlice(1), 'plotImage'), Seq.AQSlice(1).plotImage = 1; end  % plot image
if isemptyfield(Seq.AQSlice(1), 'plotPhase'), Seq.AQSlice(1).plotPhase = 0; end  % plot phase of k-space or image
if isemptyfield(Seq.AQSlice(1), 'plotB0ppm'), Seq.AQSlice(1).plotB0ppm = 0; end  % plot B0 ppm (only 3D)
if isemptyfield(Seq.AQSlice(1), 'plotB0Hz'), Seq.AQSlice(1).plotB0Hz = 0; end  % plot B0 Hz
if isemptyfield(Seq.AQSlice(1), 'ZeroFillWindowSize'), Seq.AQSlice(1).ZeroFillWindowSize = 1.4; end  % zero fill window size (high k-space values are damped by a cos^2 law)
if isemptyfield(Seq.AQSlice(1), 'ZeroFillFactor'), Seq.AQSlice(1).ZeroFillFactor = 2; end  % zero fill resolution factor


%% actual measurement
disp('measuring B0 map...');
% FIXME: Optionally, select FLASH or Spin Echo measurement
[SeqOut, mySave] = sequence_Flash(HW, Seq, varargin{:});


%% fit shim values
[dataB0Shim, AQSlice] = get_kSpaceAndImageTicks(SeqOut(1).dataB0, SeqOut(1).dataB0.AQSlice);

[rB0, p1B0, p2B0] = ndgrid(dataB0Shim.Ticks(1).ReadZ, dataB0Shim.Ticks(1).PhaseZ, dataB0Shim.Ticks(2).PhaseZ);
% Rotate to magnet coordinate system
% For alpha=phi=theta=0 -> read : z , phase(1) : x , phase(2) : y
[RzB0, RxB0, RyB0] = get_aptDegRotationMatrix(dataB0Shim.AQSlice.theta/pi*180, ...
  dataB0Shim.AQSlice.alfa/pi*180, dataB0Shim.AQSlice.phi/pi*180);
R = (RxB0' * RyB0' * RzB0')';
CoorB0 = R * [rB0(:).'; p1B0(:).'; p2B0(:).'];
xMagnet = reshape(CoorB0(2,:), size(dataB0Shim.ImageZ));
yMagnet = reshape(CoorB0(3,:), size(dataB0Shim.ImageZ));
zMagnet = reshape(CoorB0(1,:), size(dataB0Shim.ImageZ));
roi = ~isnan(dataB0Shim.ImageZ) & ...
  (abs(yMagnet) <= Seq.Shim.y_fov_max) & ...
  (sqrt(xMagnet.^2+zMagnet.^2) <= Seq.Shim.r_fov_max);
C = [ones(sum(roi(:)),1), xMagnet(roi), yMagnet(roi), zMagnet(roi)];
d = dataB0Shim.ImageZ(roi);
opt = C \ d;
disp('Deviations in B0 map:')
fprintf('B0 Offset: %8.2f %cT\n', opt(1)*1e6, 181);
fprintf('Shim X:    %8.2f %cT/m\n', opt(2)*1e6, 181);
fprintf('Shim Y:    %8.2f %cT/m\n', opt(3)*1e6, 181);
fprintf('Shim Z:    %8.2f %cT/m\n', opt(4)*1e6, 181);

% new magnet shim values
magnetShim = HW.Grad(SeqOut.AQSlice(1).iDevice).AmpOffset;
magnetShim(1:3) = magnetShim(1:3) - opt(2:4).';

shimElemStr = sprintf('%6.9f, ', magnetShim(1:3));
if numel(HW.Grad) > 1
  % FIXME: Currently shimming with channels on one single device only is supported
  newCalLine = sprintf(['HW.Grad(%d).AmpOffset([1,2,3]) = [%s]; ', ...
    ' %% %s by B0 map using sequence_Flash, x y z in T/m and B0 in T\n'], ...
    SeqOut.AQSlice(1).iDevice, shimElemStr(1:end-2), ...
    datestr(now, 'yyyy-mm-ddTHH:MM:SS'));
else
  newCalLine = sprintf(['HW.MagnetShim([1,2,3]) = [%s]; ', ...
    ' %% %s by B0 map using sequence_Flash, x y z in T/m and B0 in T\n'], ...
    shimElemStr(1:end-2), ...
    datestr(now, 'yyyy-mm-ddTHH:MM:SS'));
end

% FIXME: Add consistency checks.
if ~exist(fileparts(HW.MagnetShimPath), 'dir')
  mkdir(fileparts(HW.MagnetShimPath));
end
fid = fopen(HW.MagnetShimPath, 'a+');
if fid < 0
  warning('PD:sequence_shim_B0map:InvalidFile', ...
    'Could not open file "%s" for writing.\nThe best-fit shim values for the selected region are:', ...
    HW.MagnetShimPath);
else
  fwrite(fid, newCalLine);
  fclose(fid);
  disp('A new line was added to the following file: ');
  disp(HW.MagnetShimPath);
  fprintf('\n');
end
disp(newCalLine);


if (ishghandle(SeqOut.plotB0map) || SeqOut.plotB0map > 0)
  if SeqOut.plotB0map == 1 || ... % (1) might be used as logical (do not overwrite Figure 1)
      ~ishghandle(SeqOut.plotB0map) && ... % no HG handle and ...
      (abs(round(SeqOut.plotB0map)-SeqOut.plotB0map) > eps) % no integer
    SeqOut.plotB0map = 171;
  end
  % clear figure
  if ishghandle(SeqOut.plotB0map, 'figure')
    clf(SeqOut.plotB0map);
  elseif ~ishghandle(SeqOut.plotB0map) % FIXME: What are valid parents here? COMPLETE list!
    % only create a new window if Seq.plot is not a valid parent
    figure(SeqOut.plotB0map);
  end
end

if ishghandle(SeqOut.plotB0map, 'figure')
  % permuteOrder = SeqLoop(1).data(1).PermuteOrder;
  permuteOrder = [2, 1, 3];
  B0map = dataB0Shim.ImageZ;
  B0map(~roi) = NaN;
  hsl = sliceomatic(SeqOut.plotB0map, permute(B0map, permuteOrder), ...
    dataB0Shim.Ticks(1).ReadZ, dataB0Shim.Ticks(1).PhaseZ, dataB0Shim.Ticks(2).PhaseZ);
  % set(hsl, 'CLim', [-1e-6, 1e-6]);
  title(hsl.hAxes, 'B0 map (read corrected) in T');
  xlabel(hsl.hAxes, [AQSlice.ReadCartesianAxis{1}, 'in m']);
  ylabel(hsl.hAxes, [AQSlice.PhaseCartesianAxis{1}, 'in m']);
  zlabel(hsl.hAxes, [AQSlice.PhaseCartesianAxis{2}, 'in m']);
  title(hsl.GetSliderX(), AQSlice.ReadCartesianAxis{1});
  title(hsl.GetSliderY(), AQSlice.PhaseCartesianAxis{1});
  title(hsl.GetSliderZ(), AQSlice.PhaseCartesianAxis{2});
end

end
