function [HW, mySave, SliceSelectOut, SeqOut] = Find_Shim(HW, mySave, minTime, doPlot, Seq, SliceSelect)
%% Find magnet shim values and store in file
%
%   [HW, mySave, SliceSelectOut, SeqOut] = Find_Shim(HW, mySave, minTime, doPlot, Seq, SliceSelect)
%
% This function searches the magnet shim and saves it in HW as well as in the
% LoadMagnet file in the User folder.
% It uses fminsearch to maximize the mean signal of the selected spin echo.
%
% After the magnet shim is determined, the execution is paused for
% HW.FindFrequencyPause seconds.
%
%
% INPUT:
%
% All input parameters but HW are optional. If they are omitted or empty, the
% default values as set in the structure HW.FindShim are used.
%
%   HW        HW structure or object
%
%   mySave    mySave structure (necessary for minTime)
%
%   minTime   Minimum time in seconds since the last time Find_Shim was executed
%             before new shim values are searched again. (default: 1000)
%
%   doPlot    Plot sequence and data (bool, default: false)
%
%   Seq       structure with the following fields (default values are used when
%             omitted or empty):
%     iterations      approximate maximum number of iterations for the optimizer
%                     (default: HW.FindShim.iterations)
%     tFID            acquisition time in seconds for the FID (default: [])
%     tEcho           acquisition time in seconds for the echo
%                     (default: HW.FindShim.tEcho)
%     T1              estimated T1 time in seconds of the used sample (oil is
%                     recommended, default: HW.FindShim.T1)
%     RepetitionTime  rate in seconds for each optimization step
%                     (default: 3*Seq.T1)
%     ShimStart       1xn vector that contains the initial shim values where the
%                     search is started (in T/m; default: HW.FindShim.ShimStart)
%     ShimStep        1xn vector with initial step widths for optimization
%                     (in T/m; default: HW.FindShim.ShimStep)
%     AQEcho          relative part of the time tEcho to be acquired between the
%                     inversion pulses [0...1[ (0 => only one sample)
%                     (default: HW.FindShim.AQEcho)
%     nEchos          which echo is used for optimization (0: FID)
%                     (default: HW.FindShim.nEchos)
%     fSample         sample frequency for acquisition in Hz
%                     (default: HW.FindShim.fSample)
%     InvertPulse     function_handle for the pulse shape of the inversion pulse
%                     (default: HW.FindShim.InvertPulse)
%     SlicePulse      function_handle for the pulse shape of the slice selective
%                     excitation pulse (default: HW.FindShim.SlicePulse)
%     use_Find_Frequency_FID
%                     use Find_Frequency_FID before every iteration step
%                     (boolean, default: HW.FindShim.use_Find_Frequency_FID)
%     use_Find_Frequency_Sweep
%                     use Find_Frequency_Sweep before every iteration step
%                     (boolean, default: HW.FindShim.use_Find_Frequency_Sweep)
%     use_nEchos_Frequency
%                     use frequency of selected echo (or FID) for next iteration
%                     step (boolean). The frequency is determined by the maximum
%                     amplitude in the frequency spectrum.
%                     (default: HW.FindShim.use_nEchos_Frequency)
%     plotRaiseWindow
%                     raise and focus figure window in plot_data_1D for every
%                     iteration step (boolean, default: false)
%     plotProgress    show figure with progress of optimization
%     useFminsearch   actually do the optimization (boolean, default: true).
%                     If false, the measurements are performed without
%                     optimizing the shim values. (default: true)
%     ShimMatrixOff   Switch the ShimMatrix off during optimization.
%                     (default: true)
%     useSliceSelect  use values in SliceSelect (see below) (boolean)
%                     (default: HW.FindShim.useSliceSelect)
%     FIDWindowFunction
%                     function_handle with weighing function for optimization of
%                     the acquired signal (default: @RectWin)
%     tauMin          If the resulting T2* after shimming is lower than this
%                     value in seconds, a warning is emitted and the new shim
%                     values aren't used automatically (default: 1e-3).
%
%   SliceSelect
%             structure with the following fields. See the documentation for
%             "sequence_EchoStandard" for further details:
%     alfa            first rotation about x-axis
%     phi             seconds rotation about (new) y-axis
%     theta           third rotation about (new) z-axis
%     CenterRot       center of rotation in m; [0 0 0] is the center of the FOV.
%     nRead           not used in this sequence
%     nPhase          not used in this sequence
%     sizeRead        not used in this sequence
%     sizePhase       not used in this sequence
%     thickness       thickness of slice in meters
%
%
% OUTPUT:
%
%   HW        HW structure or object
%
%   mySave    mySave structure
%
%   sliceSelectOut
%             structure with actually used values (see input SliceSelect)
%
%   SeqOut    structure with sequence parameters and results
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
if nargin < 1  % HW
  error('PD:Find_Shim:noHW', 'The first input argument "HW" is mandatory.');
end
if nargin<2, mySave = []; end
if nargin<3 || isempty(minTime), minTime = 1000; end  % Minimal time in seconds since the last Find_Shim.
if nargin<4 || isempty(doPlot), doPlot = 0; end  % Plot sequence and data
if nargin<5, Seq = []; end
if nargin<6, SliceSelect = []; end

if isemptyfield(Seq, 'iterations'), Seq.iterations = HW.FindShim.iterations; end  % number of iterations
if isemptyfield(Seq, 'tFID'), Seq.tFID = []; end  % FID  time
if isemptyfield(Seq, 'tEcho'), Seq.tEcho = HW.FindShim.tEcho; end  % echo time
if isemptyfield(Seq, 'T1'), Seq.T1 = HW.FindShim.T1; end  % estimated T1 time of the used sample (oil is recommened)
if isemptyfield(Seq, 'RepetitionTime'), Seq.RepetitionTime = 3*Seq.T1; end

% relative part of the time tEcho to be acquired between the inversion pulses [0...1[ (0 => only one sample)
if isemptyfield(Seq, 'AQEcho'), Seq.AQEcho = HW.FindShim.AQEcho; end
if isemptyfield(Seq, 'nEchos'), Seq.nEchos = HW.FindShim.nEchos; end  % 0,1,2,3....
if isemptyfield(Seq, 'fSample'), Seq.fSample = HW.FindShim.fSample; end  % HW.RX.fSample = 125e6 Hz
if isemptyfield(Seq, 'InvertPulse'), Seq.InvertPulse = HW.FindShim.InvertPulse; end  % 180 degrees pulse function handle
if isemptyfield(Seq, 'SlicePulse'), Seq.SlicePulse = HW.FindShim.SlicePulse; end  % 90 degrees pulse function handle
% find frequency between shim steps   use_Find_Frequency_FID
if isemptyfield(Seq, 'use_Find_Frequency_FID')
  Seq.use_Find_Frequency_FID = HW.FindShim.use_Find_Frequency_FID;
end
% find frequency between shim steps   use_Find_Frequency_Sweep
if isemptyfield(Seq, 'use_Find_Frequency_Sweep')
  Seq.use_Find_Frequency_Sweep = HW.FindShim.use_Find_Frequency_Sweep;
end
% use frequency of n-th echo for next step
if isemptyfield(Seq, 'use_nEchos_Frequency')
  Seq.use_nEchos_Frequency = HW.FindShim.use_nEchos_Frequency;
end
if isemptyfield(Seq, 'plotRaiseWindow'), Seq.plotRaiseWindow = false; end  % raise and focus figure window in plot_data_1D
if isemptyfield(Seq, 'plotProgress'), Seq.plotProgress = true; end  % plot shim and mean(abs(RX)) during optimization
if isemptyfield(Seq, 'useFminsearch'), Seq.useFminsearch = true; end  % use fminsearch
if isemptyfield(Seq, 'nPreLoop'), Seq.nPreLoop = 3; end  % number of loops before fminsearch (T1)
if isemptyfield(Seq, 'ShimMatrixOff'), Seq.ShimMatrixOff = true; end  % set ShimMatrix to zero during linear shimming
if isemptyfield(Seq, 'useSliceSelect'), Seq.useSliceSelect = HW.FindShim.useSliceSelect; end  % use slice gradient
% filter window @RaisedCosine, @LorentzToGauss or @SineBel
if isemptyfield(Seq, 'FIDWindowFunction'), Seq.FIDWindowFunction = @RectWin; end
if isemptyfield(Seq, 'tauMin'),  Seq.tauMin = 1e-3;  end

if ~Seq.nEchos
  if ~isempty(Seq.tFID) && isemptyfield(Seq, 'ShimSequence');
    Seq.tEcho = Seq.tFID*2;  % echo time (if Seq.AQFID = 2 -> FID AQ Time)
    Seq.AQFID = 1;  % relative part of (echo time)/2 for data acquisition at FID ]0...1[
    Seq.AQEcho = 0.5;  % relative part of echo time for data acquistion at echos [0...1[;  0 meaning 1 sample
  end
end

if isemptyfield(Seq, 'LoadSystemAtEnd'), Seq.LoadSystemAtEnd = ~isa(HW, 'PD.HW'); end

if isemptyfield(SliceSelect, 'alfa'), SliceSelect.alfa = HW.FindShim.SliceSelect.alfa; end   % first rotation around x axis
if isemptyfield(SliceSelect, 'phi'), SliceSelect.phi = HW.FindShim.SliceSelect.phi; end   % second rotation around (new) y axis
if isemptyfield(SliceSelect, 'theta'), SliceSelect.theta = HW.FindShim.SliceSelect.theta; end   % third rotation around (new) z axis
if isemptyfield(SliceSelect, 'CenterRot'), SliceSelect.CenterRot = HW.FindShim.SliceSelect.CenterRot; end
if isemptyfield(SliceSelect, 'nRead'), SliceSelect.nRead = HW.FindShim.SliceSelect.nRead; end
if isemptyfield(SliceSelect, 'nPhase'), SliceSelect.nPhase = HW.FindShim.SliceSelect.nPhase; end
if isemptyfield(SliceSelect, 'sizeRead'), SliceSelect.sizeRead = HW.FindShim.SliceSelect.sizeRead; end
if isemptyfield(SliceSelect, 'sizePhase'), SliceSelect.sizePhase = HW.FindShim.SliceSelect.sizePhase; end
if isemptyfield(SliceSelect, 'thickness'), SliceSelect.thickness = HW.FindShim.SliceSelect.thickness; end
if isempty(SliceSelect.thickness) && Seq.useSliceSelect
  SliceSelect.thickness = 0.008;
end

if ~(Seq.useSliceSelect)
  SliceSelect.thickness = Inf;
end

% FIXME: Optionally use all available gradient channels for shimming
if isemptyfield(SliceSelect, 'iDevice'), SliceSelect.iDevice = 1; end


%% initialization
% first time clear mySave
if isempty(mySave);
  mySave.HW = HW;
  mySave.lastTime_Shim = 0;
end
if isemptyfield(mySave, 'lastTime_Shim'), mySave.lastTime_Shim = 0; end

% set ShimMatrix currents to 0
if Seq.ShimMatrixOff && evalin('base', 'exist(''ShimMatrix'', ''var'')')
  disp('Setting ShimMatrix currents to 0.');
  ShimMatrix = evalin('base', 'ShimMatrix');
  ShimMatrix.SetCurrents(zeros(48, 1));
end

% Search again only if time since last shim search has exceeded "minTime"
if (now*24*3600 - mySave.lastTime_Shim) <= minTime
  return;
end


%% actual optimization
disp('Searching shim...');

% define parameters

% Seq
Seq.plot = doPlot;  % Plot data
% length of the 1st TX pulse (excitation)
Seq.p90 = HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/HW.TX(SliceSelect.iDevice).AmpDef/2;
% length of the other TX pulses (refocus)
Seq.p180 = HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/HW.TX(SliceSelect.iDevice).AmpDef;

Seq.SliceSelect = SliceSelect;
SliceSelectOut = SliceSelect;

% run measurement
magnetShim = sequence_shim(HW, Seq);

% B0 shift off
if (isstruct(HW.Grad(SliceSelect.iDevice)) && ~isemptyfield(HW.Grad(SliceSelect.iDevice), 'ShimGradients')) || ...
    (isa(HW.Grad(SliceSelect.iDevice), 'PD.Grad') && ~isempty(HW.Grad(SliceSelect.iDevice).ShimGradients))
  magnetShim(HW.Grad(SliceSelect.iDevice).ShimGradients==0) = 0;
else
  magnetShim([0, 0, 0, ones(1,length(magnetShim)-3)]) = 0;
end

if isemptyfield(Seq, 'ShimSequence')
  ShimSequenceStr = '';
else
  ShimSequenceStr = sprintf(' using "%s"', func2str(Seq.ShimSequence));
end

[SeqOut, data_1D] = opt_maxtomeanEcho('return');
% FIXME: Support multi-channel?
iAQ = find([SeqOut.AQ(:).Device] == SeqOut.SliceSelect.iDevice, 1, 'first');
% FIXME: This always assesses the last acquisition window. Is that ok?
iNaN = find(isnan(data_1D(iAQ).data));
if numel(iNaN) == 1, iNaN = [0, iNaN]; end

% get T2*
% FIXME: This assumes that the center of excitation is at t=0. Is this always
% the case?
dt = data_1D(iAQ).data((iNaN(end-1)+1):(iNaN(end)-1));
tm = data_1D(iAQ).time_all((iNaN(end-1)+1):(iNaN(end)-1));
if SeqOut.nEchos ~= 0
  % not an FID
  dt = dt(tm>=SeqOut.tEcho*SeqOut.nEchos);
  tm = tm(tm>=SeqOut.tEcho*SeqOut.nEchos) - SeqOut.tEcho*SeqOut.nEchos;
end
fitExpSettings.CorrectFrequencyOffset = 1;
fitExpSettings.CorrectPhaseOffset = 1;
fitExpSettings.SingleExp = 1;
fitExpSettings.DoubleExp = 0;
T = fit_exp(dt(:), tm(:), fitExpSettings);
SeqOut.tau = T.tau;

shimGradIdxStr = sprintf('%d,', find(HW.Grad(SliceSelect.iDevice).ShimGradients~=0));
shimElemStr = sprintf('%6.9f, ', magnetShim(HW.Grad(SliceSelect.iDevice).ShimGradients~=0));
if numel(HW.Grad) > 1
  % FIXME: Currently shimming with channels on one single device only is supported
  newCalLine = sprintf(['HW.Grad(%d).AmpOffset([%s]) = [%s]; ', ...
    ' %% %s by FindShim%s (T2* = %.1f ms), x y z in T/m and B0 in T\n'], ...
    SliceSelect.iDevice, shimGradIdxStr(1:end-1), shimElemStr(1:end-2), ...
    datestr(now, 'yyyy-mm-ddTHH:MM:SS'), ShimSequenceStr, SeqOut.tau*1e3);
else
  newCalLine = sprintf(['HW.MagnetShim([%s]) = [%s]; ', ...
    ' %% %s by FindShim%s (T2* = %.1f ms), x y z in T/m and B0 in T\n'], ...
    shimGradIdxStr(1:end-1), shimElemStr(1:end-2), ...
    datestr(now, 'yyyy-mm-ddTHH:MM:SS'), ShimSequenceStr, SeqOut.tau*1e3);
end

useNewShim = true;
% check whether at end of shimming process, there was any significant signal
[~, phaseDiffWeightedStd] = get_MeanPhaseDiffWeighted(data_1D(iAQ).data((iNaN(end-1)+1):(iNaN(end)-1)));
fOffsetStd = phaseDiffWeightedStd * Seq.fSample/2/pi;
% FIXME: This uses the same threshold like Find_Frequency_Sweep. Do we want to
% have a separate different threshold for Find_Shim?
if fOffsetStd > HW.FindFrequencySweep.fOffsetFIDsStdMaxValue
  warning('PD:Find_Shim:fOffsetFIDStdMaxValue', ...
    'The standard deviation of the signal frequency is too high (%.3f Hz)', ...
    fOffsetStd);
  useNewShim = false;
end

if SeqOut.tau < Seq.tauMin
  warning('PD:Find_Shim:tauMinExceeded', ...
    'The T2* of the signal after shimming is too low (%.1f ms)', ...
    SeqOut.tau*1e3);
  useNewShim = false;
end

if ~useNewShim
  warning('PD:Find_Shim:discardingShim', ...
    ['Automatic shimming failed. Old shim values are still in use.\n', ...
    'Please check the sample.\nPlease check HW.fLarmor.\n', ...
    'If you want to use the new shim values anyway, add the following line to your MagnetShimCal.m file:\n', ...
    strrep(newCalLine, '%', '%%')]);
  return;
end


HW.MagnetShim = magnetShim;


%% write new shim values to LoadMagnet file
if ~isempty(HW.MagnetShimPath) && (isemptyfield(mySave, 'DummySerial') || mySave.DummySerial <= 0)
  fid = fopen(HW.MagnetShimPath, 'a+');
  fwrite(fid, newCalLine);
  fclose(fid);
  disp('A new line was added to the following file: ');
  disp(HW.MagnetShimPath);
  fprintf('\n');
  disp(newCalLine);
  fprintf('\n');
else
  if Seq.use_Find_Frequency_FID == 1
    mySave = [];
    [HW, mySave] = Find_Frequency_FID(HW, mySave, 0, 1, 0.2e-3);
  end
  disp('Please add the following line to your MagnetShimCal.m file.' );
  disp(newCalLine);
  Seq.plot = true;  % Plot data
  [~] = sequence_EchoStandard(HW, Seq);
end

mySave.lastTime_Shim = now*24*3600;
if HW.FindFrequencyPause > 0
  disp(['Waiting ' num2str(HW.FindFrequencyPause) ' seconds after finding shim']);
  sleep(HW.FindFrequencyPause);
end


%% Reload System with new parameters
if Seq.LoadSystemAtEnd
  LoadSystem;
end

end
