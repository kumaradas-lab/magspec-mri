function [SeqOut, dataOut, data, data_1D] = sequence_MPS(HW, Seq)
%% MPS experiment
%
%   [SeqOut, dataOut, data, data_1D] = sequence_MPS(HW, Seq)
%
%
% INPUT:
%
%   HW
%         HW object or structure
%
%   Seq
%         structure with properties for the sequence. Amongst others, the
%         following fields can be used. If the properties aren't set or they are
%         empty, default values may apply:
%
%     GradMPS
%           Vector with the indices for the gradients to be used where x = 1,
%           y = 2, z = 3, B = 4. (Default: 4)
%
%     IinGrad
%           Amplitude of the current through the gradient coils in Ampere. If it
%           is a scalar, the set value applies to all gradients set in
%           Seq.GradMPS. It must have the same number of elements as Seq.GradMPS
%           otherwise. Seq.IinGrad and Seq.UinGrad are mutually exclusive. Use
%           this property in case the corresponding gradient channel is current
%           controlled (e.g. the DC-600 gradient amplifier). (Default in case
%           Seq.UinGrad is not set: 1)
%
%     UinGrad
%           Amplitude of the voltage at the gradient coils in Volts. If it is a
%           scalar, the set value applies to all gradients set in Seq.GradMPS.
%           It must have the same number of elements as Seq.GradMPS otherwise.
%           Seq.IinGrad and Seq.UinGrad are mutually exclusive. Use this
%           property in case the corresponding gradient channel is voltage
%           controlled (e.g. the "Grad" of the drive-l). (Default: not set)
%
%     fGrad
%           Frequency of the gradient output in Hertz. If it is a
%           scalar, the set value applies to all gradients set in Seq.GradMPS.
%           It must have the same number of elements as Seq.GradMPS otherwise.
%           (Default: 20e3)
%
%     nPeriodGrad
%           Number of period gradients. Limits may apply (e.g. 1-100 at 20 kHz;
%           1-50 at 10 kHz). If it is a scalar, the set value applies to all
%           gradients set in Seq.GradMPS. It must have the same number of
%           elements as Seq.GradMPS otherwise. (Default: 50)
%
%     tSampleGrad
%           Sample rate for generating the gradient pulses in seconds. If it is
%           a scalar, the set value applies to all gradients set in Seq.GradMPS.
%           It must have the same number of elements as Seq.GradMPS otherwise.
%           (Default: 6e-6)
%
%     phaseGrad
%           Phase of the cosinusoidal gradient pulse in radians. If it is a
%           scalar, the set value applies to all gradients set in Seq.GradMPS.
%           It must have the same number of elements as Seq.GradMPS otherwise.
%           (Default: 0)
%
%     phaseGradIncrement
%           Increment phaseGrad in each tRep (in radians). If it is a scalar,
%           the set value applies to all gradients set in Seq.GradMPS. It must
%           have the same number of elements as Seq.GradMPS otherwise.
%           (Default: 2*pi/Seq.nMeasurements)
%
%     nPhaseGrad
%           Number of (consecutive) phase increment steps. If it is a scalar,
%           the set value applies to all gradients set in Seq.GradMPS. It must
%           have the same number of elements as Seq.GradMPS otherwise.
%           (Default: 1)
%
%     RampGrad
%           Ramp gradient amplitude up at start and down at end of each
%           tRep. The ramp time is one period of the respective channel.
%           (Default: true)
%
%     fTX
%           rf frequency at TX output in Hertz. (Default: 10e-3)
%
%     ampTX
%           rf amplitude in Tesla. (Default: 33e-6)
%
%     phaseTX
%           Phase of rf pulse in radians. (Default: 0)
%
%     channelTX
%           Used TX channel. (Default: 2)
%
%     nPrepare
%           Number of tReps without acquisition (but with the MPS gradient
%           pulses as set) before the tReps with acquisition. (Default: 1)
%
%     nMeasurements
%           Number of tReps with acquisition. (Default: prod(Seq.nPhaseGrad) )
%
%     Gain
%           Acquisition gain (must be between HW.RX.GainMin and HW.RX.GainMax).
%           (Default: HW.RX.GainMax/10)
%
%     fSampleAQ
%           Sample frequency of the acquisition in Hertz. (Default: 1e6)
%
%     fAQ
%           Mixer frequency in Hertz (for down sampling of the acquired signal).
%           (Default: 0, i.e. no mixing)
%
%     AQmode
%           Selector for the mode in which the acquisition is paused at the tRep
%           border (command load time!). If 0, the gap is as small as possible.
%           If 1, the last period of the slowest MPS frequency is not acquired.
%           For values larger than one, additionally the *first* AQmode-1
%           periods of the slowest MPS frequency are not acquired.
%
%     AQDataMode
%           Select data mode (bit width) for transmission of received signal
%           between console and PC.
%           (Default: 3*(HW.MMRT.FPGA_Firmware >= 20221129) )
%
%     average
%           Number of averages (complete experiment including Seq.nPrepare).
%           (Default: 1)
%
%     averageBreak
%           Duration of the break between averages in seconds (must be > 1 ms).
%           (Default: 1)
%
%     plotData
%           Logical value to select if the acquired data should be plotted.
%           (Default: 1)
%
%
% OUTPUT:
%
%     SeqOut
%         Structure with the actually used settings.
%
%     dataOut
%         Structure with some evaluation results for the acquired data including
%         the following fields:
%
%       nAQ
%           Number of acquired samples.
%
%       tAQ
%           Total duration of acquired data in seconds.
%
%       time_all1D
%           Time stamps for all acquired samples in seconds (if they were
%           stitched next to next).
%
%       df
%           Frequency resolution in Hertz (inverse of dataOut.tAQ).
%
%       fstart
%           Lowest frequency in spectrum in Hertz.
%
%       fstop
%           Highest frequency in spectrum in Hertz.
%
%       data1D
%           Complex normalized amplitude of the acquired signal (full scale
%           corresponds to 1).
%
%       fft
%           Fast Fourier Transform of the acquired data.
%
%     data
%         Structure with the acquired measurement data (see "get_data").
%
%     data_1D
%         Structure with the acquires measurement data (see "get_data_1D").
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2022 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

%% input check
if nargin ~= 2, error('PD:sequence_MPS:nargin', 'Number of input arguments must be 2.'); end


iDevice = 1;  % FIXME: Support multiple MMRT devices

if isemptyfield(Seq, 'GradMPS'), Seq.GradMPS = 4; end  % index of gradient connected to MPS coil
if isemptyfield(Seq, 'IinGrad') && isemptyfield(Seq, 'UinGrad')
  % output current at gradient amplifier in Amperes
  Seq.IinGrad = 1*HW.Grad(iDevice).LoadIin2Amp(Seq.GradMPS);
end
if ~isemptyfield(Seq, 'IinGrad') && ~isemptyfield(Seq, 'UinGrad')
  error('PD:sequence_MPS:IinUinSet', 'Seq.IinGrad and Seq.UinGrad are mutually exclusive.');
end
if isemptyfield(Seq, 'fGrad'), Seq.fGrad = 20e3; end  % frequency of gradient output
if isemptyfield(Seq, 'nPeriodGrad'), Seq.nPeriodGrad = 50; end  % number of gradient periods (1-100 at 20 kHz; 1-50 at 10 kHz)
Seq.tEndGrad = max(Seq.nPeriodGrad./Seq.fGrad);  % end time of oscillations on gradient channels in seconds
if isemptyfield(Seq, 'tSampleGrad'), Seq.tSampleGrad = 6e-6; end  % sample rate for generating gradient pulses in seconds
if isemptyfield(Seq, 'AQmode')
  % 0 for smallest gap, 1 for one slow period at end,
  % >1 for AQmode-1 slow periods at start + one slow period at end
  Seq.AQmode = 3+1;
end
if isemptyfield(Seq, 'AQDataMode')
  % data mode (bit width) for transmission of received signal between console and PC
  Seq.AQDataMode = 3*(HW.MMRT(iDevice).FPGA_Firmware >= 20221129);
end

if isemptyfield(Seq, 'phaseGrad'), Seq.phaseGrad = 0; end  % phase of gradient oscillation in rad
if isemptyfield(Seq, 'nPhaseGrad'), Seq.nPhaseGrad = 1; end  % number of (consecutive) phase steps
if isemptyfield(Seq, 'RampGrad'), Seq.RampGrad = true; end  % ramp gradient amplitude at start and end of tRep

if ~isfield(Seq, 'fTX'), Seq.fTX = []; end  % rf frequency at TX in Hertz
if ~isfield(Seq, 'ampTX'), Seq.ampTX = []; end  % rf amplitude at TX in Tesla
if ~isfield(Seq, 'phaseTX'), Seq.phaseTX = []; end  % phase of rf pulse in radians
if ~isfield(Seq, 'channelTX'), Seq.channelTX = []; end  % channel for rf transmission
nTX = max([numel(Seq.fTX), numel(Seq.ampTX), numel(Seq.phaseTX), numel(Seq.phaseTX)]);
if nTX > 0
  if isemptyfield(Seq, 'fTX'), Seq.fTX = 10e-3; end
  if isemptyfield(Seq, 'ampTX'), Seq.ampTX = 33e-6; end
  if isemptyfield(Seq, 'phaseTX'), Seq.phaseTX = 0; end
  if isemptyfield(Seq, 'channelTX'), Seq.channelTX = 2; end

  if numel(Seq.fTX) < nTX, Seq.fTX = ones(1, nTX) * Seq.fTX(1); end
  if numel(Seq.ampTX) < nTX, Seq.ampTX = ones(1, nTX) * Seq.ampTX(1); end
  if numel(Seq.phaseTX) < nTX, Seq.phaseTX = ones(1, nTX) * Seq.phaseTX(1); end
  if numel(Seq.channelTX) < nTX
    error('PD:sequence_MPS:UnknownChannelTX', ...
      'Seq.channelTX must be defined if multiple rf channels are used.');
  end
end

if isemptyfield(Seq, 'Gain'), Seq.Gain = HW.RX(iDevice).GainMax/10; end  % AQ Gain (must be between HW.RX.GainMin and HW.RX.GainMax)
if isemptyfield(Seq, 'fSampleAQ'), Seq.fSampleAQ = 1e6; end  % AQ sample rate in Hz
if isemptyfield(Seq, 'fAQ'), Seq.fAQ = 0; end  % AQ mixing frequency in Hz

if isemptyfield(Seq, 'nPrepare'), Seq.nPrepare = 1; end  % number of tReps without AQ before tReps with AQ
if isemptyfield(Seq, 'nMeasurements'), Seq.nMeasurements = prod(Seq.nPhaseGrad); end  % number of tReps with AQ
if isemptyfield(Seq, 'phaseGradIncrement')
  % phase increment (between subsequent acquisitions) of gradient oscillation in rad
  Seq.phaseGradIncrement = 2*pi ./ Seq.nPhaseGrad;
end

if any(Seq.phaseGradIncrement ~= 0)
  Seq.CLTime = 100e-6;
else
  Seq.CLTime = 0.8e-6;
end
if isemptyfield(Seq, 'tRep')
  % set one tRep
  Seq.tRep = Seq.tEndGrad + (max(HW.Grad(iDevice).TimeDelay(Seq.GradMPS)) + Seq.CLTime)*(Seq.AQmode~=0);
  if HW.RX(iDevice).ClampCoil.Enable
    % extent for coil blank signal if necessary
    Seq.tRep = max(Seq.tRep, ...
      Seq.tEndGrad + HW.RX(iDevice).ClampCoil.tPostset + 0.1e-3);
  end
end

if isemptyfield(Seq, 'average'), Seq.average = 1; end  % number of averages
if isemptyfield(Seq, 'averageBreak'), Seq.averageBreak = []; end  % break between averages in seconds (must be > 1 ms)
if isemptyfield(Seq, 'plotData'), Seq.plotData = 1; end  % plot acquired data
if isemptyfield(Seq, 'plotSeq'), Seq.plotSeq = Seq.GradMPS; end  % Plot the sequence, 1:4 => all Gradient outputs
if isemptyfield(Seq, 'plotSeqTR'), Seq.plotSeqTR = Seq.GradMPS; end  % Plot the sequence TR, 1:4 => all Gradient outputs

% expand gradient properties to all used gradient channels if needed
nGrads = numel(Seq.GradMPS);
if isfield(Seq, 'UinGrad') && (numel(Seq.UinGrad) < nGrads), Seq.UinGrad = Seq.UinGrad(1) * ones(size(Seq.GradMPS)); end
if isfield(Seq, 'IinGrad') && (numel(Seq.IinGrad) < nGrads), Seq.IinGrad = Seq.IinGrad(1) * ones(size(Seq.GradMPS)); end
if numel(Seq.fGrad) < nGrads, Seq.fGrad = Seq.fGrad(1) * ones(size(Seq.GradMPS)); end
if numel(Seq.nPeriodGrad) < nGrads, Seq.nPeriodGrad = Seq.nPeriodGrad(1) * ones(size(Seq.GradMPS)); end
if numel(Seq.tSampleGrad) < nGrads, Seq.tSampleGrad = Seq.tSampleGrad(1) * ones(size(Seq.GradMPS)); end
if numel(Seq.phaseGrad) < nGrads, Seq.phaseGrad = Seq.phaseGrad(1) * ones(size(Seq.GradMPS)); end
if numel(Seq.phaseGradIncrement) < nGrads, Seq.phaseGradIncrement = Seq.phaseGradIncrement(1) * ones(size(Seq.GradMPS)); end


%%
Seq.tRep = Seq.tRep * ones(1, Seq.nMeasurements+Seq.nPrepare); %  tRep

AQ.fSample = HW.RX(iDevice).fSample./round(HW.RX(iDevice).fSample./Seq.fSampleAQ);

% round Seq.fGrad such that it matches an integer number of AQ samples
Seq.fGrad = AQ.fSample./round(AQ.fSample./Seq.fGrad);

AQ.Start = [nan(1, Seq.nPrepare), ...
  ones(1, Seq.nMeasurements)*1./AQ.fSample*3+(HW.RX(iDevice).nSampleRXLatenz+2)./HW.RX(iDevice).fSample];
switch Seq.AQmode
  case 0
    AQ.nSamples = floor((Seq.tEndGrad - 100e-9 ...
      - max(max(HW.Grad(iDevice).TimeDelay(Seq.GradMPS))+Seq.CLTime+ AQ.Start(end-1), ...  % end of AQ window to CL time
            (HW.RX(iDevice).nSamplesExtraCIC)/AQ.fSample+HW.RX(iDevice).nSamplesRXExtraCIC./HW.RX(iDevice).fSample))*AQ.fSample);  % time between AQ windows (CIC)
  case 1
    AQ.nSamples = floor((Seq.tEndGrad*AQ.fSample ...  % samples in one (total) tRep
      - ceil( (Seq.CLTime+AQ.Start(end-1)+100e-9+HW.RX(iDevice).nSamplesExtra/AQ.fSample )*AQ.fSample ...  % minimum pause
              / max(AQ.fSample./Seq.fGrad) ) ...  %  samples per Grad period
        * max(AQ.fSample./Seq.fGrad)));  % samples per Grad period
  otherwise
    AQ.nSamples  = round((Seq.tEndGrad*AQ.fSample -    Seq.AQmode * max(AQ.fSample./Seq.fGrad)));
    AQ.Start = [nan(1, Seq.nPrepare), ...
      ones(1, Seq.nMeasurements) * (floor((Seq.tEndGrad*AQ.fSample(1) - AQ.nSamples(1)) - max(AQ.fSample./Seq.fGrad)) - 0.5)./AQ.fSample];
end

% acquisition must span an integer number of the slowest fGrad
% FIXME: This doesn't work correctly if the higher gradient frequencies aren't
% (integer) harmonics of the lowest gradient frequency.
maxSamplesPerPeriod = round(max(AQ.fSample./Seq.fGrad));
AQ.nSamples = floor(AQ.nSamples/maxSamplesPerPeriod) * maxSamplesPerPeriod;

Seq.MissingSamples = Seq.tEndGrad*AQ.fSample(1)-AQ.nSamples(1);
% AQ.Frequency = 125e6/6000;
AQ.Frequency = Seq.fAQ;
AQ.Phase = 0;
AQ.Gain = Seq.Gain;
AQ.Repeat = 0;
% reset DDS phase
% Important for rf pulse phase: The DDS (Direct Digital Synthesis) reference
% clock will be started at the beginning of each tRep with the frequency of the
% first rf pulse. (Make sure that that frequency is 0 by prepending a short rf
% pulse with 0 amplitude and 0 frequency.)
AQ.ResetPhases = 1;
AQ.DataMode = Seq.AQDataMode;


oldPaEnable = HW.Grad(iDevice).PaEnable;
guard = onCleanup(@() ResetPaEnable(HW, oldPaEnable, iDevice));
HW.Grad(iDevice).PaEnable = 1;

% rf transmission
[TX(1:max(1,nTX)).Start] = deal(NaN);
for iTX = 1:nTX
  TX(iTX).Channel = Seq.channelTX(iTX);
  TX(iTX).Start = repmat([-16e-9; 0] + Seq.startTX(iTX), 1, numel(Seq.tRep));
  TX(iTX).Amplitude = repmat([0; Seq.ampTX(iTX)], 1, numel(Seq.tRep));
  TX(iTX).Frequency = repmat([0; Seq.fTX(iTX)], 1, numel(Seq.tRep));
  TX(iTX).Phase = repmat([0; Seq.phaseTX(iTX)], 1, numel(Seq.tRep));
  TX(iTX).Duration = repmat([16e-9; Seq.tEndGrad-100e-9], 1, numel(Seq.tRep));
end

% differential output voltage at the Grad output (static offset)
%               |
%               v
Grad(1).Shim = 0.0 ./ HW.Grad(iDevice).Amp2LoadUin(HW.Grad(iDevice).x);  % The HW.Grad.x indexing is used because HW.Grad.Amp2PaUout is DAC channel indexed and Grad(1).Shim is [x y z B0] = [1 2 3 4] indexed.
Grad(2).Shim = 0.0 ./ HW.Grad(iDevice).Amp2LoadUin(HW.Grad(iDevice).y);  % So if you don't use HW.Grad.x indexing, you get in trouble if HW.Grad.x ~= 1 or HW.Grad.y ~= 2 or HW.Grad.z ~= 3 or HW.Grad.B ~= 4
Grad(3).Shim = 0.0 ./ HW.Grad(iDevice).Amp2LoadUin(HW.Grad(iDevice).z);
Grad(4).Shim = 0.0 ./ HW.Grad(iDevice).Amp2LoadUin(HW.Grad(iDevice).B(1));

% phaseInc = cellfun(@(inc, n) inc * (-Seq.nPrepare:1:n-1), num2cell(Seq.phaseGradIncrement), num2cell(Seq.nPhaseGrad), 'UniformOutput', false);
phaseInc = cellfun(@(inc, n) inc * (0:1:n-1), ...
  num2cell(Seq.phaseGradIncrement), num2cell(Seq.nPhaseGrad), 'UniformOutput', false);
[phaseIncAll{1:numel(Seq.GradMPS(:))}] = ndgrid(phaseInc{:});
% differential output voltage at the Grad output (dynamic)
% Generate sine wave on selected channels
for iGradMPS = 1:numel(Seq.GradMPS(:))
  Gn = Seq.GradMPS(iGradMPS);
  tGrid = 0:Seq.tSampleGrad(iGradMPS):(Seq.tEndGrad-Seq.CLTime(1));
  if (tGrid(end) - Seq.tEndGrad) == 0, tGrid(end) = []; end
  % Grad(Gn).Time = [linspace(0, Seq.tEndGrad, round(Seq.tEndGrad/Seq.tSampleGrad+1)).' + HW.Grad.TimeDelay ];
  phaseInc = [(-Seq.nPrepare:1:-1)*Seq.phaseGradIncrement(iGradMPS), ...
    reshape(phaseIncAll{iGradMPS}, 1, [])];
  Grad(Gn).Time = (tGrid.') * ones(1, size(Seq.tRep,2));
  Grad(Gn).Amp = sin(bsxfun(@plus, 2*pi*(Grad(Gn).Time(:,1))*Seq.fGrad(iGradMPS), ...
    phaseInc+Seq.phaseGrad(iGradMPS)));
  if isfield(Seq, 'UinGrad')
    Grad(Gn).Amp = Seq.UinGrad(iGradMPS) * Grad(Gn).Amp;
  else
    Grad(Gn).Amp = Seq.IinGrad(iGradMPS) * Grad(Gn).Amp;
  end
  % Grad(Gn).Repeat = [0, ones(1, Seq.nMeasurements+Seq.nPrepare-1), 0];

  if Seq.RampGrad
    % ramp up amplitude at start of tRep
    isFirstCycle = (Grad(Gn).Time(:,1) - Grad(Gn).Time(1,1)) <= 1/Seq.fGrad(iGradMPS);
    Grad(Gn).Amp(isFirstCycle,:) = Grad(Gn).Amp(isFirstCycle,:) .* ...
      sin(2*pi * Seq.fGrad(iGradMPS)/4 * bsxfun(@minus, Grad(Gn).Time(isFirstCycle,:), Grad(Gn).Time(1,:))).^2;
    % ramp down amplitude at end of tRep
    isLastCycle = (Grad(Gn).Time(end,1) - Grad(Gn).Time(:,1)) <= 1/Seq.fGrad(iGradMPS);
    Grad(Gn).Amp(isLastCycle,:) = Grad(Gn).Amp(isLastCycle,:) .* ...
      cos(2*pi * Seq.fGrad(iGradMPS)/4 * bsxfun(@minus, Grad(Gn).Time(isLastCycle,1),Grad(Gn).Time(find(isLastCycle,1, 'first'),:))).^2;
  elseif Seq.tEndGrad + 2*Seq.tSampleGrad(iGradMPS) < Seq.tRep(1)
    Grad(Gn).Amp(end,:) = 0;
  end
  Grad(Gn).Amp(end,end) = 0;
end

% some important characteristics for the gradient programming:
% - the first Amp (Grad.Amp(1,:)) is copied to the beginning of each tRep,
% - the last used Amp per tRep is copied to the end of each tRep.
% - It's not (yet) allowed to have Grad ramp to the next tRep, only a constant value is allowed.


% HW.Grad(iDevice).AmpUnit={'V' 'V' 'V' 'V'};  % Change the Text at the plots to a proper unit
% HW.Grad(iDevice).AmpUnitScale=1./HW.Grad(iDevice).Amp2PaUout(HW.Grad(iDevice).xyzB);  % Change the scale to match the units

% while 1 % start sequence and loop until pressing Strg+C
  [~, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);
%   sleep(0.2); % restart sequence after less than 2 Seconds to avoid turning of the gradients by the auto Mute
% end
% MRIUnit.MRISequency.exctractArrayToFile(talker.mySequency.getCommandArray,'test.txt')
if Seq.plotData
  Channel = 1;  % FIXME: Add support for multiple acquisition channels?
  iAQ = find([SeqOut.AQ(:).Channel] == Channel & [SeqOut.AQ(:).Device] == iDevice, 1, 'first');

  % discard data from all but the used channels
  data = data(iAQ);

  HW.RX(iDevice).AmplitudeUnitScale=1/1000/data.Amplitude2Uin(SeqOut.nPrepare+1);
  HW.RX(iDevice).AmplitudeUnit='mV';
  plot_data_1D(HW, data_1D);

  clear dataOut
  dataOut.fAQ = SeqOut.AQ(iAQ).Frequency(1,Seq.nPrepare+1);
  dataOut.fSAQ = SeqOut.AQ(iAQ).fSample(1,Seq.nPrepare+1);

  switch Seq.AQmode
    case 0
      dataOut.time_all1D = linspace(data.time_all(1,1,Seq.nPrepare+1), ...
        data.time_all(end,1,Seq.nPrepare+1), ...
        round((data.time_all(end,1,Seq.nPrepare+1) - data.time_all(1,1,Seq.nPrepare+1)) * SeqOut.AQ(iAQ).fSample(1,Seq.nPrepare+1))+1).';
      dataOut.data1D = interp1(data.time_all(~isnan(data.time_all)), ...
        data.data(~isnan(data.time_all)), dataOut.time_all1D, 'spline') ...
        * data.Amplitude2Norm(Seq.nPrepare+1);
      hf = figure(5); clf(hf);
      hax(1) = subplot(3,1,1, 'Parent', hf);
      plot(hax, dataOut.time_all1D, dataOut.data1D);

      hax(2) = subplot(3,1,2, 'Parent', hf);
      dataOut.nAQ = length(dataOut.time_all1D(:));
      dataOut.tAQ = (dataOut.time_all1D(2)-dataOut.time_all1D(1))*dataOut.nAQ;
      dataOut.df = 1/dataOut.tAQ;
      dataOut.fstart = -dataOut.df*floor(dataOut.nAQ/2)+dataOut.fAQ;
      dataOut.fstop = dataOut.df*floor(dataOut.nAQ/2-0.5)+dataOut.fAQ;
      dataOut.fft = fftshift(ifft(dataOut.data1D));
      % myPeaks=[0;((diff(abs(data.fft(1:))).*diff(abs(data.fft(2:end))))<=0) & ((diff(abs(data.fft(1:end-1))))>=0);0];
      % plot(data.f_fft1_data_interp(data.fRoi)/fmax*1e6-1e6,abs(data.fft1_data_interp(data.fRoi)))
      % hold on
      % plot(data.f_fft1_data_interp(myPeaks&data.fRoi)/fmax*1e6-1e6,abs(data.fft1_data_interp(myPeaks&data.fRoi)),'rx')
      % % text(data.f_fft1_data_interp(myPeaks&data.fRoi)/fmax*1e6-1e6, ...
      % %   abs(data.fft1_data_interp(myPeaks&data.fRoi)), ...
      % %   ['\downarrow' num2str(data.f_fft1_data_interp(myPeaks&data.fRoi)/fmax*1e6-1e6,'%10.2f').'], ...
      % %   'VerticalAlignment', 'baseline', 'HorizontalAlignment', 'center');
      % text(data.f_fft1_data_interp(myPeaks&data.fRoi)/fmax*1e6-1e6, ...
      %   abs(data.fft1_data_interp(myPeaks&data.fRoi)), ...
      %   [num2str(data.f_fft1_data_interp(myPeaks&data.fRoi).'/fmax*1e6-1e6,'%10.2f')], ...
      %   'Rotation', 90);
      % hold off

      semilogy(hax(2), linspace(dataOut.fstart, dataOut.fstop, dataOut.nAQ)/SeqOut.fGrad(1), abs(dataOut.fft), 'LineWidth', 2);
      set(hax(2), 'XTick', 0:1:dataOut.fstop/SeqOut.fGrad(1));
      % set(hax(2), 'XTickMode', 'manual');
      set(hax(2), 'XGrid', 'on');

      hax(3) = subplot(3,1,3, 'Parent', hf);
      plot(hax(3), linspace(dataOut.fstart, dataOut.fstop, dataOut.nAQ)/SeqOut.fGrad(1), angle(dataOut.fft), 'LineWidth', 2)
      set(hax(3), 'XTick', 0:1:dataOut.fstop/SeqOut.fGrad(1));
      % set(hax(3), 'XTickMode', 'manual');
      set(hax(3), 'XGrid', 'on');
      title(hax(3), ['Seq.fGrad = ' num2str(SeqOut.fGrad) ' Hz']);
      xlabel(hax(3), '"Harmonics" of Seq.fGrad(1)');

      linkaxes(hax(2:3), 'x');

    case 1
      % FIXME: The evaluation doesn't seem to be right. It would probably
      % be better to take the average of all repetitions instead of
      % "eating" the gaps. (What does that do to the spectrum?)
      dataOut.nAQ = length(data.time_all(~isnan(data.time_all)));
      dataOut.tAQ = (1/dataOut.fSAQ)*dataOut.nAQ;
      dataOut.time_all1D = linspace(1/dataOut.fSAQ,dataOut.tAQ,dataOut.nAQ).';
      dataOut.df = 1/dataOut.tAQ;
      dataOut.fstart = -dataOut.df*floor(dataOut.nAQ/2) + dataOut.fAQ;
      dataOut.fstop = dataOut.df*floor(dataOut.nAQ/2-0.5) + dataOut.fAQ;
      dataOut.data1D = data.data(~isnan(data.time_all)) * data.Amplitude2Norm(Seq.nPrepare+1);
      dataOut.fft = (fftshift(ifft(dataOut.data1D)));


      hf = figure(5); clf(hf);
      hax(1) = subplot(3,1,1, 'Parent', hf);
      plot(hax(1), dataOut.time_all1D, dataOut.data1D);

      hax(2) = subplot(3,1,2, 'Parent', hf);
      semilogy(hax(2), linspace(dataOut.fstart, dataOut.fstop, dataOut.nAQ)/SeqOut.fGrad(1), abs(dataOut.fft), 'LineWidth', 2);
      % set(hax(2), 'XTick', 0:1:dataOut.fstop/SeqOut.fGrad(1));
      % set(hax(2), 'XTickMode', 'manual');
      set(hax(2), 'XGrid', 'on');

      hax(3) = subplot(3,1,3, 'Parent', hf);
      plot(hax(3), linspace(dataOut.fstart, dataOut.fstop, dataOut.nAQ)/SeqOut.fGrad(1), angle(dataOut.fft), 'LineWidth', 2);
      % set(hax(3), 'XTick', 0:1:dataOut.fstop/SeqOut.fGrad(1));
      % set(hax(3), 'XTickMode', 'manual');
      set(hax(3), 'XGrid', 'on');
      title(hax(3), ['Seq.fGrad = ' num2str(SeqOut.fGrad) ' Hz']);
      xlabel(hax(3), '"Harmonics" of Seq.fGrad(1)');

      linkaxes(hax(2:3), 'x');
      xlim(hax(3), [0, dataOut.fstop/SeqOut.fGrad(1)]);

    otherwise
      % dataOut.nAQ = length(data.time_all(~isnan(data.time_all)));
      % dataOut.tAQ = (1/dataOut.fSAQ)*dataOut.nAQ;
      % dataOut.time_all1D = linspace(1/dataOut.fSAQ,dataOut.tAQ,dataOut.nAQ).';
      % dataOut.df = 1/dataOut.tAQ;
      % dataOut.fstart = -dataOut.df*floor(dataOut.nAQ/2) + dataOut.fAQ;
      % dataOut.fstop = dataOut.df*floor(dataOut.nAQ/2-0.5) + dataOut.fAQ;
      % dataOut.data1D = data.data(~isnan(data.time_all)) * data.Amplitude2Norm(SeqOut.nPrepare+1);

      maxSamplesPerPeriod = round(max(SeqOut.AQ(iAQ).fSample(1,SeqOut.nPrepare+1)./SeqOut.fGrad));
      dataOut.time = data.time_of_tRep((1:maxSamplesPerPeriod)',SeqOut.nPrepare+1) ...
        - data.time_of_tRep(1,SeqOut.nPrepare+1);
      dataOut.data = reshape(data.data(:,:,SeqOut.nPrepare+1:end), ...
        maxSamplesPerPeriod, ...
        size(data.time_of_tRep,1)/size(dataOut.time,1), [])*data.Amplitude2Uin(SeqOut.nPrepare+1);
      dataOut.meanData = squeeze(mean(dataOut.data,2));
      dataOut.fft1_data = ifftshift(ifft(dataOut.meanData,[],1));
      dataOut.f_fft1_data = get_FFTGrid(SeqOut.AQ(iAQ).fSample(1), maxSamplesPerPeriod);


      hf = figure(5); clf(hf);
      hax(1) = subplot(3,1,1, 'Parent', hf);
      plot(hax(1), dataOut.time, reshape(dataOut.data,maxSamplesPerPeriod,[])*1000);
      hold(hax(1), 'on');
      plot(hax(1), ...
        dataOut.time([1;end]), ...
        ones(2,1)/data.Amplitude2Norm(SeqOut.nPrepare+1) ...
        * data.Amplitude2Uin(SeqOut.nPrepare+1) * 1000 ...
        * (1+double(SeqOut.AQ(iAQ).Frequency(SeqOut.nPrepare+1)==0)), ...  % SeqOut.AQ.Frequency==0 -> double amplitude in time signal
        '--k');
      ylabel(hax(1), 'Amp in mV');
      hold(hax(1), 'off');
      grid(hax(1), 'on');

      hax(2) = subplot(3,1,2, 'Parent', hf);
      semilogy(hax(2), dataOut.f_fft1_data, abs(dataOut.fft1_data)*1000, 'LineWidth', 2);
      % set(hax(2), 'XTick', 0:1:dataOut.fstop/SeqOut.fGrad(1));
      % set(hax(2), 'XTickMode', 'manual');
      ylabel(hax(2), 'Amp in mV');
      grid(hax(2), 'on');

      hax(3) = subplot(3,1,3, 'Parent', hf);
      plot(hax(3),  dataOut.f_fft1_data, angle(dataOut.fft1_data), 'LineWidth', 2);
      % set(hax(3), 'XTick', 0:1:dataOut.fstop/SeqOut.fGrad(1));
      % set(hax(3), 'XTickMode', 'manual');
      set(hax(3), 'XGrid', 'on');
      title(hax(3), ['Seq.fGrad = ' num2str(SeqOut.fGrad) ' Hz']);
      xlabel(hax(3), '"Harmonics" of Seq.fGrad(1)');
      ylabel(hax(3), 'Phase in RAD');

      linkaxes(hax(2:3), 'x');
      xlim(hax(3), [-10, 10]*max(SeqOut.fGrad));


  end
end

end

function ResetPaEnable(HW, newValue, iDevice)
HW.Grad(iDevice).PaEnable = newValue;
end
