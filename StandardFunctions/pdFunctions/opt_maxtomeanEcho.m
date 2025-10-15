function [absmean, lastData_1D] = opt_maxtomeanEcho(HW, Seq, Shim)
%% Objective function for optimization in auto shim sequence
%
%   [absmean, data_1D] = opt_maxtomeanEcho(HW, Seq, Shim)
%
% This function acquires a Spin Echo (or FID) and returns the inverse of the
% average absolute amplitude as a measure for the optimization.
%
%
% INPUT:
%   HW
%       HW structure or object.
%       Alternatively, the first argument can be 'reset' which triggers the
%       pre-loops for the next call to this function. Or it can be 'return' in
%       which case, the function returns immediately after the output arguments
%       are set (see below).
%
%   Seq
%       Structure with sequence data which is passed to sequence_EchoStandard().
%
%   Shim
%       Free parameters for the optimization. Basically offsets to the start
%       value of the shimmming routine.
%
%
% OUTPUT:
%   absmean
%       Inverse of the average absolute amplitude of the selected Spin Echo(or
%       FID). If the first input argument is 'return', this is the last SeqOut.
%
%   data_1D
%       Structure with the 1d measurement data. If the first input argument is
%       'return', the data of the last measurement is returned.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%%
global fOffsetFIDw
persistent last_call firstLoop data_1D SeqOut


%% check input
if nargout == 2, absmean = []; lastData_1D = data_1D; end
if isempty(last_call), last_call = now; end
if ischar(HW),
  if strcmp(HW, 'reset')
    firstLoop = true;
  elseif ~strcmp(HW, 'return')
    warning('PD:opt_maxtomeanEcho:UnknownOption', 'Unknown option "%s".', HW);
  end
  absmean = SeqOut;
  return;
end


%% abort measurement
if ~isempty(Seq.guiInterruptHandle) && ishghandle(Seq.guiInterruptHandle) && ~get(Seq.guiInterruptHandle, 'Value')
  AbortMeasurement(HW);
  error('PD:opt_maxtomeanEcho:GuiAbort', 'Measurement aborted from GUI');
end


%% check heat development in coil groups due to shim
Dac2Amp_xyzB = HW.Grad(Seq.SliceSelect.iDevice).Dac2Amp(HW.Grad(Seq.SliceSelect.iDevice).xyzB);

% FIXME: Support shimming with all available gradient channels (of multiple
% devices)
usedShims = HW.Grad(Seq.SliceSelect.iDevice).ShimGradients~=0;

currentShim = Shim + Seq.ShimStart(usedShims) - Seq.ShimStep(usedShims)./0.05;

% check temperature
usedThermalGroup = HW.Grad(Seq.SliceSelect.iDevice).CoilThermalGroup(usedShims);
usedShimsIdx = find(usedShims);
for iCoilGroup = 1:max(usedThermalGroup)
  iGrad = usedShimsIdx(usedThermalGroup == iCoilGroup);
  currentCurrent = currentShim .* HW.Grad(Seq.SliceSelect.iDevice).Amp2LoadIin(iGrad);
  powerGroup = sum(HW.Grad(Seq.SliceSelect.iDevice).LoadRin(iGrad) .* currentCurrent.^2);
  if powerGroup >= HW.Grad(Seq.SliceSelect.iDevice).CoilPowerDissipation(iCoilGroup)/2
    % Power due to shim is >=50% of the power that can be dissipated.
    warning('PD:FindShim:HighShimStep', ...
      'Total power for shim in coil group %d is very high (%.3f W). Skipping step!', ...
      iCoilGroup, powerGroup);

    absmean = 1e12;
    plotProgress(HW, Seq, false);
    return;
  end
end


%% Spin Echo measurement
isHWstruct = ~isa(HW, 'PD.HWClass');

HW.Grad(Seq.SliceSelect.iDevice).AmpOffset(usedShims) = currentShim;
if isHWstruct
  HW.MagnetShim = [HW.Grad(:).AmpOffset];
end

if Seq.verbose
  shimAmpStr = sprintf(' % +8.3f,', HW.Grad(Seq.SliceSelect.iDevice).AmpOffset(usedShims)*1e3);
  shimDacStr = sprintf('% +7.1f, ', HW.Grad(Seq.SliceSelect.iDevice).AmpOffset(usedShims)./Dac2Amp_xyzB(usedShims));
  disp(['Shim', shimAmpStr(1:end-1), ' mT/m (', shimDacStr(1:end-2) ' DAC)']);
  % disp(HW.Grad(Seq.SliceSelect.iDevice).AmpOffset);
end

if Seq.use_nEchos_Frequency == 1
  if isempty(fOffsetFIDw), fOffsetFIDw = 0; end
  if ~isemptyfield(HW.FindFrequencySweep, 'shiftB0') ...
      && HW.FindFrequencySweep.shiftB0 ...
      && ~(isempty(HW.B0Target) || ~isfinite(HW.B0Target))
    newB0shift = HW.Grad(Seq.SliceSelect.iDevice).AmpOffset(4) + fOffsetFIDw / HW.GammaDef * 2*pi;
    newB0 = HW.B0Target;
    if abs(newB0shift) > ...
        HW.Grad(Seq.SliceSelect.iDevice).HoldShimNormMax(4) * HW.Grad(Seq.SliceSelect.iDevice).MaxAmp(4)
      desiredB0shift = newB0shift;
      newB0shift = sign(newB0shift) * HW.Grad(Seq.SliceSelect.iDevice).HoldShimNormMax(4) * HW.Grad(Seq.SliceSelect.iDevice).MaxAmp(4);
      newB0 = newB0 - (desiredB0shift - newB0shift);
      warning('PD:AutoFrequencyShift:B0ShiftExceeded', ...
        ['The required B0 shift (%.1f %cT) is larger than the maximum shift (%.1f %cT).\n', ...
        'Consider setting a different HW.B0Target. Adjusting HW.fLarmor to stay in range.'], ...
        desiredB0shift*1e6, 181, HW.Grad.HoldShimNormMax(4)*HW.Grad.MaxAmp(4)*1e6, 181);
    end
    HW.Grad(Seq.SliceSelect.iDevice).AmpOffset(4) = newB0shift;
    HW.B0 = newB0;
  else
    HW.fLarmor = HW.fLarmor - double(fOffsetFIDw);
  end
end
if Seq.use_Find_Frequency_FID == 1
  mySave = [];
  [HW, ~]  = Find_Frequency_FID( HW, mySave,0,1,0.2e-3);
  warning('PD:OpenMatlab:UseFindFrequency', ...
    '"Seq.use_Find_Frequency_FID" is deprecated. Use "Seq.use_nEchos_Frequency" instead.');
  warning('off', 'PD:OpenMatlab:UseFindFrequency');
end
if Seq.use_Find_Frequency_Sweep == 1
  mySave = [];
  [HW, ~] = Find_Frequency_Sweep(HW, mySave, 0, 100000, 0, HW.tFlip90Def, 20);
  warning('PD:OpenMatlab:UseFindFrequency', ...
    '"Seq.use_Find_Frequency_Sweep" is deprecated. Use "Seq.use_nEchos_Frequency" instead.');
  warning('off', 'PD:OpenMatlab:UseFindFrequency');
end
if Seq.Reinitialize
  pause(Seq.RepetitionTime - (now - last_call)*24*60*60);
  last_call = now;
end
if firstLoop
  % Pre-Loops for sample with T1 big compared to tRep
  SeqPreLoop = Seq;
  % SeqPreLoop.plot = 1;
  SeqPreLoop.IgnoreTimingError = 1;  % Don't emit warnings if timing fails in pre-loops
  plotProgress(HW, SeqPreLoop, true);  % create window before starting sequence
  if SeqPreLoop.nPreLoop > 0
    if isHWstruct
      fLarmor_orig = HW.fLarmor;
    end
    SeqPreLoop.Reinitialize = 1;
    for iPreLoop = 1:SeqPreLoop.nPreLoop
      if isemptyfield(Seq, 'ShimSequence')
        [dataPre, SeqOutPre] = sequence_EchoStandard(HW, SeqPreLoop, Seq.SliceSelect);
      else
        [SeqOutPre, ~] = Seq.ShimSequence(HW, SeqPreLoop);
        dataPre = SeqOutPre.data;
      end
      if Seq.use_nEchos_Frequency
        iAQ = find([SeqOutPre.AQ(:).Device] == Seq.SliceSelect.iDevice, 1, 'first');  % FIXME: Support multi-channel?
        if isHWstruct
          old_fOffsetFIDw = fOffsetFIDw;
        end
        fOffsetFIDw = get_MeanPhaseDiffWeighted(dataPre(iAQ).data(1:SeqOutPre.AQ(iAQ).nSamples(1,SeqOutPre.nEchos+1),1,SeqOutPre.nEchos+1)) * ...
          SeqOutPre.AQ(iAQ).fSample(1,SeqOutPre.nEchos+1) / 2/pi;
        if isHWstruct
          fOffsetFIDw = fOffsetFIDw + old_fOffsetFIDw;
          HW.fLarmor = fLarmor_orig;
        end
        if ~isemptyfield(HW.FindFrequencySweep, 'shiftB0') ...
            && HW.FindFrequencySweep.shiftB0 ...
            && ~(isempty(HW.B0Target) || ~isfinite(HW.B0Target))
          newB0shift = HW.Grad(Seq.SliceSelect.iDevice).AmpOffset(4) + fOffsetFIDw / HW.GammaDef * 2*pi;
          newB0 = HW.B0Target;
          if abs(newB0shift) > ...
              HW.Grad(Seq.SliceSelect.iDevice).HoldShimNormMax(4) * HW.Grad(Seq.SliceSelect.iDevice).MaxAmp(4)
            desiredB0shift = newB0shift;
            newB0shift = sign(newB0shift) * HW.Grad(Seq.SliceSelect.iDevice).HoldShimNormMax(4) * HW.Grad(Seq.SliceSelect.iDevice).MaxAmp(4);
            newB0 = newB0 - (desiredB0shift - newB0shift);
            warning('PD:AutoFrequencyShift:B0ShiftExceeded', ...
              ['The required B0 shift (%.1f %cT) is larger than the maximum shift (%.1f %cT).\n', ...
              'Consider setting a different HW.B0Target. Adjusting HW.fLarmor to stay in range.'], ...
              desiredB0shift*1e6, 181, HW.Grad.HoldShimNormMax(4)*HW.Grad.MaxAmp(4)*1e6, 181);
          end
          HW.Grad(Seq.SliceSelect.iDevice).AmpOffset(4) = newB0shift;
          HW.B0 = newB0;
        else
          HW.fLarmor = HW.fLarmor - double(fOffsetFIDw);
        end
      end
      SeqPreLoop.Reinitialize = Seq.Reinitialize;
    end
    clear SeqOutPre dataPre
  end
end
firstLoop = false;

if isemptyfield(Seq, 'ShimSequence')
  [data, SeqOut, data_1D] = sequence_EchoStandard(HW, Seq, Seq.SliceSelect);
else
  [SeqOut, ~] = Seq.ShimSequence(HW, Seq);
  data = SeqOut.data;
  [SeqOut, data_1D] = get_data_1D(SeqOut, data, 1);
  % plot_data_1D(HW, data_1D);
end

lastData_1D = data_1D;

iAQ = find([SeqOut.AQ(:).Device] == Seq.SliceSelect.iDevice, 1, 'first');  % FIXME: Support multi-channel?
if Seq.use_nEchos_Frequency
  if isHWstruct
    old_fOffsetFIDw = fOffsetFIDw;
  end
  fOffsetFIDw = get_MeanPhaseDiffWeighted(data(iAQ).data(1:SeqOut.AQ(iAQ).nSamples(1,SeqOut.nEchos+1),1,SeqOut.nEchos+1)) * ...
    SeqOut.AQ(iAQ).fSample(1,SeqOut.nEchos+1) / 2/pi;
  if isHWstruct
    fOffsetFIDw = fOffsetFIDw + old_fOffsetFIDw;
  end
end
echo = abs(data(iAQ).data(:,1,Seq.nEchos+1)) .* Seq.FIDWindowFunction(numel(data(iAQ).data(:,1,Seq.nEchos+1)));
absmean = 1/mean(echo);


%% plot graph with development of the optimization parameters and target
plotProgress(HW, Seq, false);

  function plotProgress(HW, Seq, discardPoints)
    persistent hg_props
    if ishghandle(Seq.plotProgress, 'figure') || ishghandle(Seq.plotProgress, 'uipanel') || ...
        (mod(Seq.plotProgress,1) == 0 && Seq.plotProgress > 0)
      if discardPoints
        plot_data{1,1} = NaN(2*sum(usedShims), 1);
        plot_data{2,1} = NaN;
      else
        % shim values
        OffsetinDAC_xyzB = HW.Grad(Seq.SliceSelect.iDevice).OffsetinDAC(HW.Grad(Seq.SliceSelect.iDevice).xyzB);
        plot_data{1,1} = [HW.Grad(Seq.SliceSelect.iDevice).AmpOffset(usedShims), ...
          (round(HW.Grad(Seq.SliceSelect.iDevice).AmpOffset(usedShims) ./ Dac2Amp_xyzB(usedShims) + ...
          OffsetinDAC_xyzB(usedShims)) - OffsetinDAC_xyzB(usedShims)) .* ...
          Dac2Amp_xyzB(usedShims)];
        % Mean RX B1
        if absmean == 1e12
          plot_data{2,1} = NaN;
        else
          plot_data{2,1} = 1/HW.RX(Seq.SliceSelect.iDevice).AmplitudeUnitScale/absmean;
        end
      end
      if isempty(hg_props)
        % shim values
        GradColors = HW.PlotSequence.GradColors;
        for iColors = 1:length(GradColors)
          if iscellstr(GradColors{iColors})
            GradColors{iColors} = str2rgb(GradColors{iColors});
          end
          if iscell(GradColors{iColors})
            GradColors{iColors} = cell2mat(reshape(GradColors{iColors}, [], 1));
          end
        end
        colorOrder = cellfun(@mean, GradColors, 'UniformOutput', false);
        colorOrder = cat(1, colorOrder{:});
        hg_props.axes{1}.props.ColorOrder = colorOrder(1:sum(usedShims),:);
        hg_props.axes{1}.props.LineStyleOrder = {'-', '--'};
        hg_props.axes{1}.props.Tag = 'shim_amp';
        hg_props.axes{1}.xlabel = 'Iteration';
        hg_props.axes{1}.ylabel = 'Shim in T/m';
        hg_props.axes{1}.legend.labels = HW.Grad(Seq.SliceSelect.iDevice).Name(usedShims);
        hg_props.axes{1}.legend.props.Box = 'off';
        hg_props.axes{1}.legend.props.Orientation = 'vertical';
        hg_props.axes{1}.legend.props.Location = 'northeast';

        % Mean RX B1
        hg_props.axes{2}.xlabel = 'Iteration';
        hg_props.axes{2}.ylabel = sprintf('Mean Abs %s in %s', ...
          HW.RX(Seq.SliceSelect.iDevice).AmplitudeName, HW.RX(Seq.SliceSelect.iDevice).AmplitudeUnit);
        hg_props.axes{2}.props.YLim = [0 Inf];
        hg_props.axes{2}.props.Tag = 'shim_mean_rx_amp';
        % all axes
        hg_props.allAxes.props.XLim = [0, Seq.iterations];
        hg_props.allAxes.props.Box = 'on';
        hg_props.allAxes.props.XMinorGrid = 'on';
        hg_props.allAxes.props.YMinorGrid = 'on';
        % figure
        if ishghandle(Seq.plotProgress) && isprop(Seq.plotProgress, 'Name')
          hg_props.figure.props.Name = 'Shim Progress';
        end
      end

      plot_data_store(Seq.plotProgress, plot_data, hg_props);
    end
  end

end
