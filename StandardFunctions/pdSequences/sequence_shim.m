function shimopt = sequence_shim(HW, Seq)
%% Search the magnet shim using fminsearch
%
%       shimopt = sequence_shim(HW, Seq)
%
% The actual measurement is done with the function "opt_maxtomeanEcho". By
% default, that function calls "sequence_EchoStandard" (or the function set in
% Seq.ShimSequence). See also the documentation of the used function for
% additional parameters that can be used as fields of Seq.
%
%
% INPUT:
%
%   HW
%         HW structure or object.
%
%   Seq
%         structure with the following fields (default values are used when
%         omitted or empty):
%
%     T1
%           approximate T1 relaxation time of the used sample in seconds.
%           (Default: HW.FindShim.T1)
%
%     RepetitionTime
%           Repetition time (between iteration steps) in seconds.
%           (Default: 3*Seq.T1)
%
%     nEchos
%           number of the echo that is used for the optimization (0: FID).
%           (Default: HW.FindShim.nEchos)
%
%     ShimStart
%           Start the optimization process with these values (see
%           HW.MagnetShim). This can be either a row vector with values for each
%           gradient channel, or a scalar in which case it applies to all
%           gradient channels.
%           (Default: HW.FindShim.ShimStart or if that is empty HW.MagnetShim)
%
%     ShimStep
%           Initial step size used in optimization. This can be either a row
%           vector with values for each gradient channel, or a scalar in which
%           case it applies to all gradient channels.
%           (Default: If Seq.ShimStart is non-zero, [0.001, 0.001, 0.001, 0.0000002].
%                     Otherwise, [0.005, 0.005, 0.005, 0.000001])
%
%     iterations
%           Approximate maximum number of iterations for fminsearch.
%           (Default: HW.FindShim.iterations)
%
%     nPreLoop
%           Number of loops before the optimization process with fminsearch
%           starts. This can be used to approximate a steady-state of the
%           magnetization. (Default: 3)
%
%     useFminsearch
%           If false, the measurement is repeated without actually optimizing
%           with fminsearch. (Default: true)
%
%     FIDWindowFunction
%           Function handle for the filter window for the acquisition that is
%           used in the optimization process. Possible filter window functions
%           include: @RectWin, @RaisedCosine, @LorentzToGauss or @SineBel
%           (Default: @RectWin)
%
%     use_nEchos_Frequency
%           Track the Larmor frequency during the optimization and use the
%           frequency that was detected in each respectively previous step.
%           (Default: HW.FindShim.use_nEchos_Frequency)
%
%     use_Find_Frequency_FID
%           Use the function "Find_Frequency_FID" before the optimization
%           process to detect the Larmor frequency.
%           (Default: HW.FindShim.use_Find_Frequency_FID)
%
%     use_Find_Frequency_Sweep
%           Use the function "Find_Frequency_Sweep" before the optimization
%           process to detect the Larmor frequency.
%           (Default: HW.FindShim.use_Find_Frequency_Sweep)
%
%     verbose
%           Show information about each shimming step in command window.
%           (Default: true).
%
%     plot
%           Plot pulse sequence and measured data of each step.
%           (Default: false)
%
%     plotRaiseWindow
%           Raise the window with the measured data at each step.
%           (Default: false)
%
%     plotProgress
%           Plot graph showing progress of optimization.
%           (Default: true)
%
%     useSliceSelect
%           Use slice selective pulse sequence for the measurement. Whether this
%           is actually used depends on the use shim sequence (see
%           Seq.ShimSequence). The settings in Seq.SliceSelect are used.
%           (Default: false)
%
%     SliceSelect
%           Structure with settings for slice selection (see
%           Seq.useSliceSelect). Default values are used when omitted or empty:
%       alfa
%             First rotation around magnet x-axis.
%             (Default: HW.FindShim.SliceSelect.alfa)
%       phi
%             Secound rotation around rotated y'-axis.
%             (Default: HW.FindShim.SliceSelect.phi)
%       theta
%             Third rotation around rotated z''-axis.
%             (Default: HW.FindShim.SliceSelect.theta)
%       thickness
%             Thickness of slice in meter.
%             (Default: HW.FindShim.SliceSelect.thickness)
%
%
% OUTPUT:
%
%   shimopt
%         Magnet shim values with "best" measurement results
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

%% check input
if nargin < 2, error('PD:sequence_shim:NoInput', 'Function must be called with two inputs (HW, Seq).'); end

%% default parameters
if isemptyfield(Seq, 'T1'), Seq.T1 = HW.FindShim.T1; end  % estimated T1 time of the used sample in seconds (oil is recommened)
if isemptyfield(Seq, 'RepetitionTime'), Seq.RepetitionTime = 3*Seq.T1; end % repetition time in seconds
if isemptyfield(Seq, 'LoopsRepetitionTime'), Seq.LoopsRepetitionTime = Seq.RepetitionTime; end % repetition time in seconds
if isemptyfield(Seq, 'nEchos'), Seq.nEchos = HW.FindShim.nEchos; end  % number of the echo that is used for the optimization (0: FID)

if isemptyfield(Seq, 'iterations'), Seq.iterations = HW.FindShim.iterations; end  % (maximum) number of iterations
if isemptyfield(Seq, 'nPreLoop'), Seq.nPreLoop = 3; end  % number of loops before fminsearch (T1)
if isemptyfield(Seq, 'useFminsearch'), Seq.useFminsearch = true; end  % use fminsearch (dummy run otherwise)
if isemptyfield(Seq, 'FIDWindowFunction'), Seq.FIDWindowFunction = @RectWin; end  % filter window @RectWin, @RaisedCosine, @LorentzToGauss or @SineBel
if isemptyfield(Seq, 'use_nEchos_Frequency'), Seq.use_nEchos_Frequency = HW.FindShim.use_nEchos_Frequency; end   % use frequency of n-th echo for next step
if isemptyfield(Seq, 'use_Find_Frequency_FID'), Seq.use_Find_Frequency_FID = HW.FindShim.use_Find_Frequency_FID; end   % find frequency between shim steps   use_Find_Frequency_FID
if isemptyfield(Seq, 'use_Find_Frequency_Sweep'), Seq.use_Find_Frequency_Sweep = HW.FindShim.use_Find_Frequency_Sweep; end   % find frequency between shim steps   use_Find_Frequency_Sweep
if isemptyfield(Seq, 'verbose'),  Seq.verbose = true;  end

if isemptyfield(Seq, 'plot'), Seq.plot = 0; end  % Plot sequence and data
if isemptyfield(Seq, 'plotRaiseWindow'), Seq.plotRaiseWindow = false; end  % raise and focus figure window in plot_data_1D
if isemptyfield(Seq, 'plotProgress'), Seq.plotProgress = true; end  % plot shim and mean(abs(RX)) during optimization
if ~isfield(Seq, 'guiInterruptHandle'), Seq.guiInterruptHandle = []; end

if isemptyfield(Seq,'useSliceSelect'), Seq.useSliceSelect = HW.FindShim.useSliceSelect; end % use slice gradient
if ~isfield(Seq, 'SliceSelect'), Seq.SliceSelect = struct(); end
if isemptyfield(Seq.SliceSelect, 'alfa'), Seq.SliceSelect.alfa = HW.FindShim.SliceSelect.alfa; end   % around x axis
if isemptyfield(Seq.SliceSelect, 'phi'), Seq.SliceSelect.phi = HW.FindShim.SliceSelect.phi; end   % around y axis
if isemptyfield(Seq.SliceSelect, 'theta'), Seq.SliceSelect.theta = HW.FindShim.SliceSelect.theta; end   % around z axis
if isemptyfield(Seq.SliceSelect, 'CenterRot'), Seq.SliceSelect.CenterRot = HW.FindShim.SliceSelect.CenterRot; end
if isemptyfield(Seq.SliceSelect, 'nRead'), Seq.SliceSelect.nRead = HW.FindShim.SliceSelect.nRead; end
if isemptyfield(Seq.SliceSelect, 'nPhase'), Seq.SliceSelect.nPhase = HW.FindShim.SliceSelect.nPhase; end
if isemptyfield(Seq.SliceSelect, 'sizeRead'), Seq.SliceSelect.sizeRead = HW.FindShim.SliceSelect.sizeRead; end
if isemptyfield(Seq.SliceSelect, 'sizePhase'), Seq.SliceSelect.sizePhase = HW.FindShim.SliceSelect.sizePhase; end
if isemptyfield(Seq.SliceSelect, 'thickness'), Seq.SliceSelect.thickness = HW.FindShim.SliceSelect.thickness; end
if isempty(Seq.SliceSelect.thickness)
  if Seq.useSliceSelect
    Seq.SliceSelect.thickness = 0.008;
  end
end
if ~Seq.useSliceSelect
  Seq.SliceSelect.thickness = 1e12;
end

% FIXME: Optionally use all available gradient channels for shimming
if isemptyfield(Seq.SliceSelect, 'iDevice'), Seq.SliceSelect.iDevice = 1; end

if isemptyfield(Seq, 'ShimStart')  % start with these shim values
  if isempty(HW.FindShim.ShimStart)
    Seq.ShimStart = HW.Grad(Seq.SliceSelect.iDevice).AmpOffset;
  else
    deviceOffset = cumsum([0, HW.Grad(:).n]);
    Seq.ShimStart = HW.FindShim.ShimStart(deviceOffset(Seq.SliceSelect.iDevice) + (1:HW.Grad(Seq.SliceSelect.iDevice).n));
  end
end
if isscalar(Seq.ShimStart)
  Seq.ShimStart = ones(1, HW.Grad(Seq.SliceSelect.iDevice).n) * Seq.ShimStart;
end
if isemptyfield(Seq, 'ShimStep')  % stepwidth
  if isempty(HW.FindShim.ShimStep)
    if any(Seq.ShimStart ~= 0)
      % If ShimStart is set, use smaller step width (assuming we are close to
      % the correct shim values).
      Seq.ShimStep = [0.001, 0.001, 0.001, repelem(0.0000002, 1, numel(HW.Grad(Seq.SliceSelect.iDevice).B))];
    else
      Seq.ShimStep = [0.005, 0.005, 0.005, repelem(0.000001, 1, numel(HW.Grad(Seq.SliceSelect.iDevice).B))];
    end
  else
    deviceOffset = cumsum([0, HW.Grad(:).n]);
    Seq.ShimStep = HW.FindShim.ShimStep(deviceOffset(Seq.SliceSelect.iDevice) + (1:HW.Grad(Seq.SliceSelect.iDevice).n));
  end
end
if isscalar(Seq.ShimStep)
  Seq.ShimStep = ones(1, HW.Grad(Seq.SliceSelect.iDevice).n) * Seq.ShimStep;
end

%%
global fOffsetFIDw
% HW.MagnetShim=Seq.ShimStart;
if Seq.use_nEchos_Frequency == 1
  fOffsetFIDw = 0;
end

% don't reinitialize the mri device
if HW.Grad(Seq.SliceSelect.iDevice).HoldShim || (Seq.RepetitionTime <= 2)
  Seq.TimeToNextSequence = ...
    floor((2^HW.MMRT(Seq.SliceSelect.iDevice).CLTimeBitWidth-12) ...
          / HW.MMRT(Seq.SliceSelect.iDevice).fSystem) - 3;
  Seq.IgnoreTimingError = true;
  if isemptyfield(Seq, 'ShimSequence')
    [data, SeqOut] = sequence_EchoStandard(HW, Seq, Seq.SliceSelect);
    Seq.plotAllHandle = SeqOut.plotAllHandle;
  else
    SeqOut = Seq.ShimSequence(HW, Seq);
    data = SeqOut.data;
  end
  Seq.TimeFromLastSequence = Seq.LoopsRepetitionTime - sum(SeqOut.tRep);
  Seq.TimeToNextSequence = Seq.TimeFromLastSequence;
  Seq.CLTime = SeqOut.CLTime*0 + 10e-6;
  % Seq.tOffset = HW.Grad(Seq.SliceSelect.iDevice).tEC + HW.Grad(Seq.SliceSelect.iDevice).tRamp;
  Seq.Reinitialize = false;
  Seq.IgnoreTimingError = false;
  Seq.plotSeq = [];
  if Seq.use_nEchos_Frequency == 1
    iAQ = find([SeqOut.AQ(:).Device] == Seq.SliceSelect.iDevice, 1, 'first');  % FIXME: Support multi-channel?
    fOffsetFID = get_MeanPhaseDiffWeighted(data(iAQ).data(1:SeqOut.AQ(iAQ).nSamples(1,SeqOut.nEchos+1),1,SeqOut.nEchos+1)) * ...
      SeqOut.AQ(iAQ).fSample(1,SeqOut.nEchos+1) / 2/pi;
    HW.fLarmor = HW.fLarmor - double(fOffsetFID);
  end
else
  oldReInit = HW.ReInit;
  reInitProtect = onCleanup(@() setfield(HW, 'ReInit', oldReInit));  %#ok<SFLD>
  HW.ReInit = 1;
  Seq.Reinitialize = 1;
end


if (ishghandle(Seq.plot) || Seq.plot > 0) && ~Seq.plotRaiseWindow
  if Seq.plot == 1 || ... % (1) might be used as logical (do not overwrite Figure 1)
      ~ishghandle(Seq.plot) && ... % no HG handle and ...
      (abs(round(Seq.plot)-Seq.plot) > eps) % no integer
    Seq.plot = 100;
  end
  % clear figure
  if ishghandle(Seq.plot, 'figure')
    % focus the window here once instead of in every loop
    figure(Seq.plot);
  elseif ~ishghandle(Seq.plot) % FIXME: What are valid parents here? COMPLETE list!
    % only create a new window if Seq.plot is not a valid parent
    figure(Seq.plot);
  end
end

if ishghandle(Seq.plotProgress) || Seq.plotProgress > 0
  if Seq.plotProgress == 1 || ... % (1) might be used as logical (do not overwrite Figure 1)
     ~ishghandle(Seq.plotProgress) && ... % no HG handle and ...
      (abs(round(Seq.plotProgress)-Seq.plotProgress) > eps) % no integer
    Seq.plotProgress = 200;
  end
  keepProgressPos = false;
  % clear figure
  if ishghandle(Seq.plotProgress, 'figure')
    clf(Seq.plotProgress);
    % focus the window here once instead of in every loop
    figure(Seq.plotProgress);
    keepProgressPos = true;
  elseif ~ishghandle(Seq.plotProgress) % FIXME: What are valid parents here? COMPLETE list!
    % only create a new window if Seq.plotProgress is not a valid parent
    figure(Seq.plotProgress);
  end

  % Move figure to uncover "Timeline continuous"
  if ~keepProgressPos && (ishghandle(Seq.plot) || Seq.plot > 0)
    if ishghandle(Seq.plot, 'figure')
      % "Timeline continuous" figure already open
      tcPosition = get(Seq.plot, 'Position');
    elseif ishghandle(Seq.plot)
      % "Timeline continuous" is plotted as part of another window
      tcPosition = [];
    else
      % "Timeline continuous" not yet open
      % Open hidden figure to get default position
      hf = figure('Visible', 'off');
      tcPosition = get(hf, 'Position');
      close(hf)
    end
    if ~isempty(tcPosition)
      % Check if "Shim Progress" covers "Timeline continuous"
      spPosition = get(Seq.plotProgress, 'Position');
      coversX = sum(spPosition([1,3])) >= tcPosition(1) && spPosition(1) <= sum(tcPosition([1,3]));
      coversY = sum(spPosition([2,4])) >= tcPosition(2) && spPosition(2) <= sum(tcPosition([2,4]));
      if coversX && coversY
        % Get resolution of monitor with lower left corner of the figure
        monitorPos = get(0, 'MonitorPositions');
        onMonitor = all(bsxfun(@lt, tcPosition(1:2), monitorPos(:,1:2)+monitorPos(:,3:4)) & ...
          bsxfun(@ge, tcPosition(1:2), monitorPos(:,1:2)), 2);
        thisMonitorPos = monitorPos(onMonitor,:);
        % Only if there is a minimum number of free pixels on this screen
        minFreeWidth = 250;
        if tcPosition(1) > thisMonitorPos(1) + minFreeWidth || tcPosition(3) + minFreeWidth < thisMonitorPos(3)
          % check available pixels on left or right of "Timeline continuous"
          if tcPosition(1)-thisMonitorPos(1) < sum(thisMonitorPos([1,3])) - sum(tcPosition([1,3]))
            % move to right
            spPosition(1) = sum(tcPosition([1,3]))+1;
            if sum(spPosition([1,3])) > sum(thisMonitorPos([1,3]))
              % resize figure to fit on screen
              spPosition(3) = sum(thisMonitorPos([1,3])) - spPosition(1) - 1;
            end
          else
            % move to left
            spPosition(1) = tcPosition(1) - spPosition(3) - 1;
            if spPosition(1) < thisMonitorPos(1)
              % resize figure to fit on screen
              spPosition(1) = thisMonitorPos(1);
              spPosition(3) = tcPosition(1) - thisMonitorPos(1) - 1;
            end
          end
          set(Seq.plotProgress, 'Position', spPosition)
        end
      end
      if ~Seq.plotRaiseWindow
        figure(Seq.plot);
      end
    end
  end

  % clear persistent variables in function
  plot_data_store('clear');
end

opt_maxtomeanEcho('reset');
MagnetShimOld = HW.MagnetShim;
fLarmorOld = HW.fLarmor;
% optimization using fminsearch
if Seq.useFminsearch
  warning('on', 'PD:OpenMatlab:UseFindFrequency');
  Dac2Amp_xyzB = HW.Grad(Seq.SliceSelect.iDevice).Dac2Amp(HW.Grad(Seq.SliceSelect.iDevice).xyzB);

  uwp = [];
  if isa(HW, 'PD.HWClass')
    % Unwind protect (for user pressing Ctrl+C during operation)
    uwp = PD.UnwindProtectGuard(@() RestoreHWGuard(HW, fLarmorOld, MagnetShimOld));
  end

  usedShims = HW.Grad(Seq.SliceSelect.iDevice).ShimGradients~=0;
  if any(abs(Seq.ShimStep(usedShims)./Dac2Amp_xyzB(usedShims)) > 6400)
    warning('PD:FindShim:HighStep', 'The initial step size is very high (%d DAC)!', ...
      round(max(abs(Seq.ShimStep(usedShims)./Dac2Amp_xyzB(usedShims)))));
  end

  shimopt_offset = fminsearch(@(Shim) opt_maxtomeanEcho(HW, Seq, Shim), ...
      Seq.ShimStep(usedShims)./0.05, ...
      optimset('MaxFunEvals', Seq.iterations, 'Display', 'off', ...
               'TolX', min(abs(Dac2Amp_xyzB(usedShims)))*sqrt(3/2), ...
               'TolFun', inf));
  uwp.isUnwindProtect = false;

  shimopt = Seq.ShimStart - Seq.ShimStep./0.05;
  shimopt(usedShims) = shimopt(usedShims) + shimopt_offset;
  shimopt(~usedShims) = 0;
else
  for iLoop = 1:Seq.iterations
    opt_maxtomeanEcho(HW, Seq, Seq.ShimStep./0.05);
  end
  shimopt = Seq.ShimStart;
end

if Seq.use_nEchos_Frequency == 1
  clear global fOffsetFIDw
end
%--- end sequence_shim
end

function RestoreHWGuard(HW, fLarmor, MagnetShim)
HW.fLarmor = fLarmor;
HW.MagnetShim = MagnetShim;
end
