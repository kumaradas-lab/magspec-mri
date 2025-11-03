function [data, SeqOut, data_1D] = sequence_EchoStandard(HW, Seq, SliceSelectUser, temp)
%% Acquire Spin Echo or CPMG Echo train
%
%  [data, SeqOut, data_1D] = sequence_EchoStandard(HW, Seq, SliceSelect)
%
%
% INPUT:
%
%   HW
%           HW structure or object (see LoadSystem).
%
%   Seq
%           structure containing the settings measurement sequence. All
%           settings are optional and default values are used if they are
%           omitted. Among others the following fields can be set:
%
%     tEcho
%             Echo time in seconds (default: 10e-3).
%
%     tEchoFirst
%             Echo time of first echo in seconds that might deviate from
%             Seq.tEcho (default: Seq.tEcho).
%
%     tEchoFirstTau2
%             Optionally, an excitation pulse can be added in front of the
%             actual CPMG echo train. This is the time in seconds before the
%             excitation that starts the CPMG echo train. (default: 0, i.e. no
%             additional pulse).
%
%     AQEcho
%             Relative part of Seq.tEcho during which the echo is acquired
%             (default: 0.8).
%
%     tAQEcho
%             Acquisition time at echoes between the inversion pulses in
%             seconds. If 0, only one sample is acquired. For negative values,
%             no samples are acquired. If Seq.AQEcho is also defined, it takes
%             precedence.
%
%     AQFID
%             Relative part of the time between excitation pulse and first
%             inversion pulse during which the FID is acquired (default:
%             Seq.AQEcho).
%
%     nEchos
%             Number of echoes in the echo train (default: 1).
%
%     tEchoTrain
%             Duration of echo train in seconds. If Seq.nEchos also is defined,
%             it takes precedence.
%
%     FlipPulse
%             Pulse shape function for the excitation pulse. See Pulse_... for
%             examples. (Default: @Pulse_Rect)
%
%     InvertPulse
%             Pulse shape function for the inversion pulses. See Pulse_... for
%             examples. (Default: @Pulse_Rect)
%
%     FirstInvertPulse
%             Pulse shape function for the excitaion pulse before the actual
%             CPMG echo train. See Pulse_... for  examples.
%             (Default: Seq.FlipPulse)
%
%     SlicePulse
%             Pulse shape function for the excitation pulse that is used instead
%             of Seq.FlipPulse if Seq.useSliceSelect is true. See Pulse_... for
%             examples. (Default: Seq.FlipPulse)
%
%     useSliceSelect
%             Boolean value to set if a slice is selected (default: false).
%
%     Function_Prepare_Measurement
%             Function handle that is executed before the measurement is started
%             with the following function signature
%               [HW, Seq, AQ, TX, Grad] = @(HW, Seq, AQ, TX, Grad);
%             It can be used to modify the pulse sequence after the (default)
%             pulse sequence was created.
%
%
%   SliceSelect
%           Structure with data for the slice selection. Among others, the
%           following fields can be used:
%
%     alfa
%             First rotation around axis x in radians.
%
%     phi
%             Second rotation around axis y' in radians.
%
%     theta
%             Third rotation around axis z'' in radians.
%
%     Center2OriginImage
%             1x3 vector with the distance measured from the center of the image
%             to the origin in the magnet (in image coordinates) in meters.
%
%     thickness
%             Thickness of the slice in meters.
%
%
% OUTPUT:
%
%   data
%           Structure with the sorted measurement data. See "get_data" for more
%           details.
%
%   SeqOut
%           Structure with the actually used sequence settings.
%
%   data_1D
%           Structure with the serialized measurement data. See "get_data_1D"
%           for more details.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% FIXME: Complete documentation (including phase cycling using Seq.SeqAverage).

% Check if first argument is a talker.
% Please do not use talker as first argument.
if isa(HW, 'PD.Talker') || ...
    (isstruct(HW) && ~isfield(HW, 'MMRT'))
  if nargin >= 2, HW = Seq; end
  if nargin >= 3, Seq = SliceSelectUser; end
  if nargin >= 4, SliceSelectUser = temp; end
  nargint = nargin-1;
else
  nargint = nargin;
end

if nargint < 3
  % default slice selection
  SliceSelectUser.alfa = 0*pi;  % 1st rotation about x axis
  SliceSelectUser.phi = 0*pi;  % 2nd rotation about rotated y axis
  SliceSelectUser.theta = -0.5*pi;  % 3rd rotation about rotated z axis
  SliceSelectUser.CenterRot = [0, 0, 0];
  SliceSelectUser.nRead = 1;
  SliceSelectUser.nPhase = 1;
  SliceSelectUser.sizeRead = Inf;  % 1e12;
  SliceSelectUser.sizePhase = Inf;  % 1e12;
  SliceSelectUser.thickness = Inf;  % 1e12;
end

if ~isempty(whos('global', 'SliceSelect'))
  % A global SliceSelect structure overrides the local argument
  global SliceSelect
else
  SliceSelect = SliceSelectUser;
end

% Create standard parameters if missing
if nargint==1, Seq = []; end
Seq = set_EmptyField(Seq, 'PreProcessSequence', 1);
Seq = set_EmptyField(Seq, 'StartSequence', 1);
Seq = set_EmptyField(Seq, 'PollPPGfast', 1);
Seq = set_EmptyField(Seq, 'GetRawData', 1);
Seq = set_EmptyField(Seq, 'PostProcessSequence', 1);

if Seq.PreProcessSequence
  % slice selection parameters
  if isemptyfield(Seq, 'useSliceSelect'), Seq.useSliceSelect = false; end
  if ~Seq.useSliceSelect
    % overwrite with default slice selection
    SliceSelect.alfa = 0*pi;  % 1st rotation about x axis
    SliceSelect.phi = 0*pi;  % 2nd rotation about rotated y axis
    SliceSelect.theta = 0.5*pi; % 3rd rotation about rotated z axis
    SliceSelect.CenterRot = [0, 0, 0];
    SliceSelect.nRead = 1;
    SliceSelect.nPhase = 1;
    SliceSelect.sizeRead = Inf;  % 1e12;
    SliceSelect.sizePhase = Inf;  % 1e12;
    SliceSelect.thickness = Inf;  % 1e12;
  end

  if isempty(whos('global', 'SliceSelect'))  % FIXME: Why this condition?
    SliceSelect.Center = SliceSelect.CenterRot([2,3,1]).*[1,-1,1];
    SliceSelect.CenterRot = SliceSelect.Center([3,1,2]).*[1,1,-1];
    [SliceSelect.normV(1), SliceSelect.normV(2), SliceSelect.normV(3)] = ...
      sph2cart(SliceSelect.theta, -SliceSelect.phi, 1);
    SliceSelect.R = SliceSelect.CenterRot * SliceSelect.normV.';
    SliceSelect.Rauf = SliceSelect.R * SliceSelect.normV;
    SliceSelect.RaufCenter = SliceSelect.CenterRot - SliceSelect.Rauf;
    SliceSelect.CenterRauf = SliceSelect.Rauf - SliceSelect.CenterRot;
    SliceSelect.CenterRaufImage = tpaRotate(SliceSelect.CenterRauf', ...
      -SliceSelect.alfa ,-SliceSelect.phi, -SliceSelect.theta)';
    SliceSelect.CenterRaufImage = SliceSelect.CenterRaufImage([3,2,1]).*[1,-1,1];
    SliceSelect.Center2OriginImage = [SliceSelect.CenterRaufImage([1,2]), SliceSelect.R];
  end

  SliceSelect.ReadOS = 1;
  SliceSelect.PhaseOS = 1;
  SliceSelect.Flip = pi/2;
  if isemptyfield(SliceSelect, 'iDevice'), SliceSelect.iDevice = 1; end  % index for device  % FIXME: Support multiple devices?
  if isemptyfield(SliceSelect, 'MaxGradAmpSlice')
    SliceSelect.MaxGradAmpSlice = HW.Grad(SliceSelect.iDevice).MaxAmpSlice;
  end

  Seq = set_EmptyField(Seq, 'tEcho', 10e-3);                                % Echo time [s]
  Seq = set_EmptyField(Seq, 'tEchoFirst', Seq.tEcho);                       % First Echo time [s]
  Seq = set_EmptyField(Seq, 'tEchoFirstTau2', 0);                           % temporal separation of First Echo Pulse  [s]
  Seq = set_EmptyField(Seq, 'tOffset', 0);                                  % Seq offset time [s]
  Seq = set_EmptyField(Seq, 'tRepMin', 0);                                  % minimum repetition time [s]
  Seq = set_EmptyField(Seq, 'plot', true);                                  % use plot_data_1D
  if ~isfield(Seq, 'plotAllHandle'),    Seq.plotAllHandle   = [];     end   % axis handle for plot_data_1D
  Seq = set_EmptyField(Seq, 'plotRaiseWindow', true);                       % raise and focus figure window in plot_data_1D
  Seq = set_EmptyField(Seq, 'plotTR', false);                               % use plot_data_1D_TR
  Seq = set_EmptyField(Seq, 'TXAmp', HW.TX(SliceSelect.iDevice).AmpDef);  % TX amplitude of inversion pulses in T
  Seq = set_EmptyField(Seq, 'TXAmpP90', Seq.TXAmp);                         % TX amplitude of excitation pulse in T
  if isemptyfield(Seq, 'p90')
    % duration of 90 degrees pulse at Seq.TXAmpP90 in seconds
    Seq.p90 = HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/Seq.TXAmpP90/2;
  end
  if isemptyfield(Seq, 'p180')
    % duration of 180 degrees pulse at Seq.TXAmp in seconds
    Seq.p180 = HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/Seq.TXAmp;
  end
  Seq = set_EmptyField(Seq, 'tAQEcho', 0.8 * Seq.tEcho);                    % acquisition time at echoes between the inversion pulses [0...x s (0 => only one sample, negative no sample)
  Seq = set_EmptyField(Seq, 'AQEcho', Seq.tAQEcho/Seq.tEcho);               % relative part of the time tEcho to be acquired between the inversions pulses [0...1[ (0 => only one sample)
  Seq = set_EmptyField(Seq, 'AQFID', Seq.AQEcho);                           % AQFID relative part of the time tEcho/2 to be acquired between the first pulse and the inversion Pulse [0...1[ (0 => only one sample, negative no sample)
  Seq = set_EmptyField(Seq, 'AQFIDGrid', 0);                                % AQFID start on 0s + fSample grid
  Seq = set_EmptyField(Seq, 'AQFrequency', HW.fLarmor);                     % mixing frequency of the receiver for all windows
  Seq = set_EmptyField(Seq, 'TXFrequency', HW.fLarmor);                     % frequency of the transmitter for all pulses
  Seq = set_EmptyField(Seq, 'AQPhaseOffset', 0);                            % phase offset of the receiver for all windows
  Seq = set_EmptyField(Seq, 'TXPhaseOffset', 0);                            % phase offset of the transmitter for all pulses
  Seq = set_EmptyField(Seq, 'tEchoTrain', Seq.tEcho*1);                     % duration of echo Train [s]
  Seq = set_EmptyField(Seq, 'nEchos', round(Seq.tEchoTrain/Seq.tEcho));     % Number of acquired echoes
  if ~isfield(Seq, 'tsT2T2'),           Seq.tsT2T2          = [];     end   % separation time between CPMG pulse trains
  Seq = set_EmptyField(Seq, 'tEchoTrainT2T2', Seq.nEchos*Seq.tEcho);        % duration of echo Train [s]
  Seq = set_EmptyField(Seq, 'nEchosT2T2', round(Seq.tEchoTrainT2T2/Seq.tEcho));  % number of Echoes of second CPMG pulse trains
  Seq = set_EmptyField(Seq, 'fSample', HW.RX(SliceSelect.iDevice).fSample/1250);  % Sampling rate of the AQ windows at the echoes
  Seq = set_EmptyField(Seq, 'fSampleFID', Seq.fSample);                     % Sampling rate of the AQ window at the FID
  Seq = set_EmptyField(Seq, 'Shim', zeros(1, HW.Grad(SliceSelect.iDevice).n));  % Shim add to HW.MagnetShim
  % set shim for additional gradient channels to zero
  Seq.Shim((numel(Seq.Shim)+1):HW.Grad(SliceSelect.iDevice).n) = 0;
  Seq = set_EmptyField(Seq, 'average', 1);  % Number of averages
  Seq = set_EmptyField(Seq, 'averageBreak', HW.FindFrequencyPause);  % Break between averages
  Seq = set_EmptyField(Seq, 'TXdelay', 0);  % Add delay to all TX pulses
  Seq = set_EmptyField(Seq, 'TXdelay90', 0);  % Add delay to the TX 90 degrees pulse
  Seq = set_EmptyField(Seq, 'InvertPulse', @Pulse_Rect);  % 180 degrees pulse function handle
  Seq = set_EmptyField(Seq, 'FlipPulse', @Pulse_Rect);  % 90 degrees pulse function handle
  Seq = set_EmptyField(Seq, 'FirstInvertPulse', Seq.FlipPulse);  % 90 degrees pulse function handle for separated refocus pulse
  Seq = set_EmptyField(Seq, 'SlicePulse', Seq.FlipPulse);  % 90 degrees pulse function handle
  Seq.tFlip = Seq.p90 * Seq.FlipPulse(HW, 'Amp');  % Calculate the length needed for the 90 degrees pulse
  Seq.tInvertFirst = Seq.p90 * Seq.FirstInvertPulse(HW, 'Amp');  % Calculate the length needed for the 90 degrees pulse
  Seq.tInvert = Seq.p180 * Seq.InvertPulse(HW, 'Amp');  % Calculate the length needed for the 180 degrees pulse
  Seq = set_EmptyField(Seq, 'TxFlipBW', 1/Seq.tFlip* Seq.FlipPulse(HW,'Time'));   % Bandwidth of the RF pulses
  Seq = set_EmptyField(Seq, 'TxFirstInvertBW', 1/Seq.tInvertFirst* Seq.FirstInvertPulse(HW,'Time')); % Bandwidth of the RF pulses
  Seq = set_EmptyField(Seq, 'TxFlipSteps', 51);                             % Number of pulses used for the TX shape
  Seq = set_EmptyField(Seq, 'fast', 0);                                     % If 1: The hole sequence will be put into the first Repetition time. Helps recuding echo time.
  Seq = set_EmptyField(Seq, 'RawData', 0);
  if Seq.AQFID >= 0, Seq.RawData = 0; end  % Raw data mode doesn't work when also acquiring the FID.
  Seq = set_EmptyField(Seq, 'EndSpoiltStart', Seq.tEcho);
  Seq = set_EmptyField(Seq, 'EndSpoilDuration', Seq.tEcho/2);
  Seq = set_EmptyField(Seq, 'EndSpoilAmp', [0,0,0]);
  Seq = set_EmptyField(Seq, 'StartSequence', 1);
  % structure with settings for damping coil after TX pulse
  Seq = set_EmptyField(Seq, 'DampCoil', struct());
  Seq.DampCoil = set_EmptyField(Seq.DampCoil, 'Enable', false);  % Damp coil before AQ window (after 180 degrees pulses)
  Seq.DampCoil = set_EmptyField(Seq.DampCoil, 'Delay', 0);  % Delay after 180 degrees pulse before damping signal
  Seq.DampCoil = set_EmptyField(Seq.DampCoil, 'Duration', 9e-6);  % Duration of damp signal
  Seq.DampCoil = set_EmptyField(Seq.DampCoil, 'DOChannel', 1);  % Digital Output channel that switches damping circuit
  Seq.DampCoil = set_EmptyField(Seq.DampCoil, 'DampP90', false);  % Also damp after 90 degrees pulse
  Seq.DampCoil = set_EmptyField(Seq.DampCoil, 'DurationP90', Seq.DampCoil.Duration);  % Duration of damp signal after 90 degrees pulse

  % prepare function with special support for prepolarization or DNP pulse
  if ~isfield(Seq, 'Function_Prepare_Measurement'), Seq.Function_Prepare_Measurement = []; end
  if isemptyfield(Seq, 'Prepol'), Seq.Prepol = struct(); end  % settings for prepolarization pulse
  if isemptyfield(Seq.Prepol, 'use')
    % add pre-polarization pulse to pulse program
    if isemptyfield(HW.DefSeqValues, {'Prepol', 'use'})
      Seq.Prepol.use = false;
    else
      Seq.Prepol.use = HW.DefSeqValues.Prepol.use;
    end
  end

  if isemptyfield(Seq, 'DNP'), Seq.DNP = struct(); end  % settings for DNP pulse
  if isemptyfield(Seq.DNP, 'use')
    % add DNP pulse to pulse program
    if isemptyfield(HW.DefSeqValues, {'DNP', 'use'})
      Seq.DNP.use = false;
    else
      Seq.DNP.use = HW.DefSeqValues.DNP.use;
    end
  end

  if Seq.DNP.use && Seq.Prepol.use
    warning('PD:sequence:EchoStandard:PolarizationConflict', ...
      'Both the prepolarization pulse and the DNP pulse are active. The prepolarization pulse takes precedence.');
  end

  if Seq.Prepol.use
    if ~exist('Prepol_Prepare_Sequence_Spin_Echo', 'file')
      error('PD:sequence_EchoStandard:noPrepolPrepare', ...
        'Pre-polarization pulse can only be added with the necessary prepare function.');
    else
      Seq.Function_Prepare_Measurement = @Prepol_Prepare_Sequence_Spin_Echo;
    end
  elseif Seq.DNP.use
    if ~exist('DNP_Prepare_Sequence_Spin_Echo', 'file')
      error('PD:sequence_EchoStandard:noDNPPrepare', ...
        'DNP pulse can only be added with the necessary prepare function.');
    else
      Seq.Function_Prepare_Measurement = @DNP_Prepare_Sequence_Spin_Echo;
    end
  end

  if ~isfield(Seq, 'SeqAverage'), Seq.SeqAverage = []; end  % use average
  if ~isempty(Seq.SeqAverage)
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'average', 1);  % number of average
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'averageBreak', HW.FindFrequencyPause); %
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'SaveSeqAverageData', 1);       %
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'TXPhaseInvertIncrement', 180); %
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'TXPhaseFlipIncrement', 0);     %
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'AQPhaseIncrement', 0);         %
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'RandomTXRXPhaseOffset', 0);    %
    if Seq.SeqAverage.RandomTXRXPhaseOffset
      [~, IX] = sort(rand(1, Seq.SeqAverage.average));
      Seq.SeqAverage.RandomTXRXPhaseOffset=linspace(0,360-360./Seq.SeqAverage.average,Seq.SeqAverage.average);
      Seq.SeqAverage.RandomTXRXPhaseOffset=Seq.SeqAverage.RandomTXRXPhaseOffset(IX);
      Seq.SeqAverage.TXPhaseFlipOffset=Seq.SeqAverage.RandomTXRXPhaseOffset;
      Seq.SeqAverage.TXPhaseInvertOffset=Seq.SeqAverage.RandomTXRXPhaseOffset;
      Seq.SeqAverage.AQPhaseOffset=Seq.SeqAverage.RandomTXRXPhaseOffset;
    else
      Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'TXPhaseInvertOffset', 0);            %
      Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'TXPhaseFlipOffset', 0);              %
      Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'AQPhaseOffset', 0);                  %
      if numel(Seq.SeqAverage.TXPhaseInvertOffset) < Seq.SeqAverage.average
        Seq.SeqAverage.TXPhaseInvertOffset  = repmat(Seq.SeqAverage.TXPhaseInvertOffset, 1, Seq.SeqAverage.average);
      end
      if numel(Seq.SeqAverage.TXPhaseFlipOffset) < Seq.SeqAverage.average
        Seq.SeqAverage.TXPhaseFlipOffset    = repmat(Seq.SeqAverage.TXPhaseFlipOffset, 1, Seq.SeqAverage.average);
      end
      if numel(Seq.SeqAverage.AQPhaseOffset) < Seq.SeqAverage.average
        Seq.SeqAverage.AQPhaseOffset        = repmat(Seq.SeqAverage.AQPhaseOffset, 1, Seq.SeqAverage.average);
      end
    end
  end

  if ~isfield(Seq, 'plotData1D'), Seq.plotData1D = []; end

  Seq.tInvertmax=1000e-6;
  Seq.tAQmax=Seq.tEcho/2;
  Seq.tTxSlicemax = max(HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/(HW.TX(SliceSelect.iDevice).AmpDef*0.9)*(pi/2)/pi * Seq.SlicePulse(HW, 'Amp'), ...
    1/(SliceSelect.MaxGradAmpSlice*HW.GammaDef/2/pi*SliceSelect.thickness) * Seq.SlicePulse(HW, 'Time'));
  Seq.TxSliceBW = 1/Seq.tTxSlicemax*Seq.SlicePulse(HW, 'Time');
  Seq.tSlice = Seq.tTxSlicemax + 2*HW.Grad(SliceSelect.iDevice).tRamp + 2*HW.Grad(SliceSelect.iDevice).tEC;
  Seq.tGrad = Seq.tSlice/2 + HW.Grad(SliceSelect.iDevice).tRamp/2;
  Seq.SliceSelect.nPhase3D = 0;
  t = Seq.fSample;
  Seq = createSeq(Seq, SliceSelect, HW);
  Seq.fSample = t;

  % Gradient not used; only shim
  for t = 1:HW.Grad(SliceSelect.iDevice).n
    Grad(t).Time = NaN;
    Grad(t).Amp = 0;
    Grad(t).Shim = Seq.Shim(t);
  end

  if Seq.useSliceSelect
    GradTime = cumsum([...
      0; ...  % slice start
      HW.Grad(SliceSelect.iDevice).tRamp; ...
      Seq.tTxSlicemax + 2*HW.Grad(SliceSelect.iDevice).tEC; ...
      HW.Grad(SliceSelect.iDevice).tRamp; ...  % slice stop
      HW.Grad(SliceSelect.iDevice).tEC; ...
      HW.Grad(SliceSelect.iDevice).tRamp; ...  % Grads start
      Seq.tGrad - 2*HW.Grad(SliceSelect.iDevice).tRamp; ...
      HW.Grad(SliceSelect.iDevice).tRamp; ...  % Grads stop
      0]);
    GradTime = GradTime - Seq.tSlice/2;
    GradTimeEnd = GradTime(end-1) + HW.Grad(SliceSelect.iDevice).tEC;
    GradTimeSart = GradTime(1) + HW.Grad(SliceSelect.iDevice).tEC;
    z = zeros(1, 1);
    o = ones(1, 1);
    for t = 1:3
      Grad(t).Time = [GradTime*o, [zeros(1,Seq.nEchos); NaN(8,Seq.nEchos)]];
      Grad(t).Amp=[[...
        z; ...
        Seq.AmpSlice(t)*o; ...
        Seq.AmpSlice(t)*o; ...
        z; ...
        z; ...
        Seq.AmpSliceD(t); ...
        Seq.AmpSliceD(t); ...
        z; ...
        z],               zeros(9, Seq.nEchos)];
      % tt = tt+1;
    end

    if any(Seq.EndSpoilAmp)
      if numel(Seq.EndSpoilAmp)==1, Seq.EndSpoilAmp(1:3)=Seq.EndSpoilAmp(1); end
        EndSpoilGradTime=cumsum([...
                            Seq.EndSpoiltStart; ...  % spoil start
                            HW.Grad(SliceSelect.iDevice).tRamp; ...
                            Seq.EndSpoilDuration - HW.Grad(Seq.iDevice).tRamp; ...
                            HW.Grad(SliceSelect.iDevice).tRamp; ...  % spoil stop
                            HW.Grad(SliceSelect.iDevice).tEC ...
                        ]);
        z = zeros(1, 1);
        o = ones(1, 1);
        for t = 1:3
          Grad(t).Time = [Grad(t).Time; [NaN(5,Seq.nEchos), EndSpoilGradTime*o]];
          Grad(t).Amp = [Grad(t).Amp; [NaN(5,Seq.nEchos), [z; ...
                                                           Seq.EndSpoilAmp(t)*o; ...
                                                           Seq.EndSpoilAmp(t)*o; ...
                                                           z;...
                                                           z]]];
        end
    end
  else
    GradTimeEnd = 0;
    GradTimeSart = 0;
  end

  Seq.tRep = [Seq.tEcho/2, Seq.tEcho*ones(1,Seq.nEchos)];  % generate the repetition times
  % set the sampling rates of the AQ windows
  AQ.fSample = HW.RX(SliceSelect.iDevice).fSample ./ ...
    round(HW.RX(SliceSelect.iDevice).fSample ./ [Seq.fSampleFID, Seq.fSample*ones(1,Seq.nEchos)]);
  if Seq.nEchos==1 && Seq.AQEcho >= 1
    % special case: if only one echo is generated, the AQ window can be much longer
    Seq.tRep = [Seq.tEcho/2, Seq.tEcho*ones(1,Seq.nEchos)*Seq.AQEcho+1000e-6+Seq.tEcho/2]; % extend the second repetition time
    AQShift = (Seq.AQEcho-1)*Seq.tEcho/2 + Seq.tEcho/2/8;               % shift AQ of echo
    AQ1 = 0;
  else
    AQShift = 0;
    if numel(AQ.fSample) >= 2
      Seq.AQEcho = round(Seq.AQEcho*Seq.tEcho*AQ.fSample(2))./AQ.fSample(2)./Seq.tEcho;
    end
    AQ1 = 1;
  end
  if Seq.nEchos==0
    % special case: only the FID is generated, the AQ window can be much longer
    Seq.tRep = max(Seq.tEcho/2*Seq.AQFID+1/AQ.fSample(1)*5+500e-6, Seq.tRepMin); % extend the repetition time
  else
    Seq.tRep(1) = Seq.tRep(1) + Seq.tEchoFirst/2 - Seq.tEcho/2 - Seq.tEchoFirstTau2/2;
    Seq.tRep(2) = Seq.tRep(2) + Seq.tEchoFirst/2 - Seq.tEcho/2 + Seq.tEchoFirstTau2/2;
  end

  if Seq.useSliceSelect;
    if any(Seq.EndSpoilAmp)
      Seq.tRep(end) = max(Seq.tRep(end), EndSpoilGradTime(end)-GradTimeSart+max(HW.Grad(SliceSelect.iDevice).TimeDelay)+1e-3);
      for t = 1:3
        [~, IX] = sort(Grad(t).Time, 1);
        IX = IX + cumsum(ones(size(Grad(t).Time)),2).*size(Grad(t).Time,1) - size(Grad(t).Time,1);
        Grad(t).Time(:) = Grad(t).Time(IX);
        Grad(t).Amp(:) = Grad(t).Amp(IX);
      end
      clear t IX
    end
  end

  AQ.Frequency = Seq.AQFrequency;  % set the mixing frequency of the receiver for all windows

  if Seq.useSliceSelect
    GradTime(end) = Seq.tRep(1) - Seq.tSlice/2;  %!!!!!!!!!!!!!!!
  end

  % TX
  if isscalar(Seq.TXPhaseOffset)
    Seq.TXPhaseOffset = repmat(Seq.TXPhaseOffset, 1, numel(Seq.tRep));
    TX.Repeat = [0, zeros(1,double(Seq.nEchos>0)), ones(1,max(0,Seq.nEchos-1))];
  end
  if ischar(Seq.TXPhaseOffset)
    if strcmp(Seq.TXPhaseOffset, 'rand')
      Seq.TXPhaseOffset = 360 * rand(1, numel(Seq.tRep));
    end
    if strcmp(Seq.TXPhaseOffset, 'linear')
      Seq.TXPhaseOffset = 360 * linspace(0, 1, numel(Seq.tRep));
    end
    if strcmp(Seq.TXPhaseOffset, 'alternate')
      Seq.TXPhaseOffset = [0, 90+90*(-1).^(1:numel(Seq.tRep)-1)];
    end

  end

  if Seq.useSliceSelect
    % get properties for 90 degrees slice selective excitation pulse
    pulseData1 = Seq.SlicePulse(HW, Seq.TXdelay+Seq.TXdelay90, Seq.TxSliceBW, ...
      pi*Seq.p90/(HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/Seq.TXAmpP90), Seq.TxFlipSteps, ...
      Seq.tTxSlicemax, Seq.fTxSlice, 180);
  else
    % get properties for 90 degrees excitation pulse
    pulseData1 = Seq.FlipPulse(HW, Seq.TXdelay+Seq.TXdelay90, Seq.TxFlipBW, ...
      pi*Seq.p90/(HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/Seq.TXAmpP90), Seq.TxFlipSteps, ...
      Seq.tFlip, Seq.TXFrequency, 180);
  end
  % get properties for 180 degrees inversion pulses
  pulseData2 = Seq.InvertPulse(HW, Seq.TXdelay, 1/Seq.tInvert*Seq.InvertPulse(HW,'Time'), ...
    pi*Seq.p180/(HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/Seq.TXAmp), 20, ...
    Seq.tInvert, Seq.TXFrequency, -90);

  if Seq.tEchoFirstTau2
    % get properties for second 90 degrees pulse (first excitation [T1])
    pulseData3 = Seq.FirstInvertPulse(HW, Seq.TXdelay, Seq.TxFirstInvertBW, ...
      pi*Seq.p90/(HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/Seq.TXAmpP90), ...
      Seq.TxFlipSteps, Seq.tInvertFirst, Seq.TXFrequency, -90);
    TX.Start = NaN(max([size(pulseData1.Start,1), size(pulseData2.Start,1), 2*size(pulseData3.Start,1)]), 1+Seq.nEchos);
  else
    TX.Start = NaN(max(size(pulseData1.Start,1), size(pulseData2.Start,1)), 1+Seq.nEchos);
  end

  % Prepare matrix
  TX.Duration = TX.Start;
  TX.Amplitude = TX.Start;
  TX.Frequency = TX.Start;
  TX.Phase = TX.Start;

  % Insert data of 90 degrees pulse in matrix
  TX.Duration(1:size(pulseData1.Start,1),1) = pulseData1.Duration;
  TX.Start(1:size(pulseData1.Start,1),1) = pulseData1.Start;
  TX.Amplitude(1:size(pulseData1.Start,1),1) = pulseData1.Amplitude;
  TX.Frequency(1:size(pulseData1.Start,1),1) = pulseData1.Frequency;
  TX.Phase(1:size(pulseData1.Start,1),1) = pulseData1.Phase + Seq.TXPhaseOffset(1);

  % Insert data of 180 degrees pulse in matrix
  if ~Seq.tEchoFirstTau2
    if Seq.nEchos >= 1
      TX.Duration(1:size(pulseData2.Start,1),2:end) = pulseData2.Duration*ones(1,Seq.nEchos);
      TX.Start(1:size(pulseData2.Start,1),2:end) = pulseData2.Start*ones(1,Seq.nEchos);
      TX.Amplitude(1:size(pulseData2.Start,1),2:end) = pulseData2.Amplitude*ones(1,Seq.nEchos);
      TX.Frequency(1:size(pulseData2.Start,1),2:end) = pulseData2.Frequency*ones(1,Seq.nEchos);
      TX.Phase(1:size(pulseData2.Start,1),2:end) = pulseData2.Phase*ones(1,Seq.nEchos)+repmat(Seq.TXPhaseOffset(2:end),size(pulseData2.Phase,1),1);
    end
  else
    if Seq.nEchos >= 1
      TX.Duration(1:size(pulseData3.Start,1),2) = pulseData3.Duration;
      TX.Start(1:size(pulseData3.Start,1),2) = pulseData3.Start;
      TX.Amplitude(1:size(pulseData3.Start,1),2) = pulseData3.Amplitude;
      TX.Frequency(1:size(pulseData3.Start,1),2) = pulseData3.Frequency;
      TX.Phase(1:size(pulseData3.Start,1),2) = pulseData3.Phase+repmat(Seq.TXPhaseOffset(2),size(pulseData3.Phase,1),1);

      TX.Duration(size(pulseData3.Start,1)+1:2*size(pulseData3.Start,1),2) = pulseData3.Duration;
      TX.Start(size(pulseData3.Start,1)+1:2*size(pulseData3.Start,1),2) = pulseData3.Start+Seq.tEchoFirstTau2;
      TX.Amplitude(size(pulseData3.Start,1)+1:2*size(pulseData3.Start,1),2) = pulseData3.Amplitude;
      TX.Frequency(size(pulseData3.Start,1)+1:2*size(pulseData3.Start,1),2) = pulseData3.Frequency;
      TX.Phase(size(pulseData3.Start,1)+1:2*size(pulseData3.Start,1),2) = pulseData3.Phase+repmat(Seq.TXPhaseOffset(2),size(pulseData3.Phase,1),1);
    end
    if Seq.nEchos > 1
      TX.Duration(1:size(pulseData2.Start,1),3:end) = pulseData2.Duration*ones(1,Seq.nEchos-1);
      TX.Start(1:size(pulseData2.Start,1),3:end) = pulseData2.Start*ones(1,Seq.nEchos-1);
      TX.Amplitude(1:size(pulseData2.Start,1),3:end) = pulseData2.Amplitude*ones(1,Seq.nEchos-1);
      TX.Frequency(1:size(pulseData2.Start,1),3:end) = pulseData2.Frequency*ones(1,Seq.nEchos-1);
      TX.Phase(1:size(pulseData2.Start,1),3:end) = pulseData2.Phase*ones(1,Seq.nEchos-1)+repmat(Seq.TXPhaseOffset(3:end),size(pulseData2.Phase,1),1);
    end

  end

  % Damp Coil
  if Seq.DampCoil.Enable
    if HW.TX(SliceSelect.iDevice).DampCoil.Enable
      warning('PD:sequence_EchoStandard:DoubleDamp', ...
        'Damping is activated in HW.TX.DampCoil. Are you sure you want to include additional damping on the sequence level?');
    end
    if ~isemptyfield(Seq, 'DigitalIO')
      warning('PD:sequence_EchoStandard:DampCoil', 'Damping Coil not supported if Seq.DigitalIO is already defined')
      Seq.DampCoil.Enable = false;
    else
      % switch Digital IO on after every pulse
      Seq.DigitalIO.SetTime = repelem(TX.Start + TX.Duration,2,1) + repmat([0; Seq.DampCoil.Duration],size(TX.Start)) + Seq.DampCoil.Delay;
      Seq.DigitalIO.SetValue = repmat([2.^(Seq.DampCoil.DOChannel-1); 0], size(TX.Start));
      Seq.DigitalIO.Repeat = [0, zeros(1,double(Seq.nEchos>1)), ones(1,max(0,Seq.nEchos-1))];
      if ~Seq.DampCoil.DampP90
        % Do not damp 90 degrees pulse
        Seq.DigitalIO.SetTime(:,1) = NaN;
        Seq.DigitalIO.SetValue(:,1) = NaN;
      else
        % Set special duration after 90 degrees pulse
        Seq.DigitalIO.SetTime(2,1) = Seq.DigitalIO.SetTime(1,1) + Seq.DampCoil.DurationP90;
      end
    end
  end

  % AQ
  if numel(Seq.AQPhaseOffset) == 1
    Seq.AQPhaseOffset = repmat(Seq.AQPhaseOffset,1,numel(Seq.tRep));
    AQ.Repeat = [0, zeros(1,double(Seq.nEchos>0)), ones(1,max(0,Seq.nEchos-1))];
  else
    AQ.Repeat = [0, zeros(1,max(0,Seq.nEchos))];
  end
  if ischar(Seq.AQPhaseOffset)
    if strcmp(Seq.AQPhaseOffset, 'rand')
      Seq.AQPhaseOffset = 360 * rand(1, numel(Seq.tRep));
    end
    if strcmp(Seq.AQPhaseOffset, 'linear')
      Seq.AQPhaseOffset = -360*linspace(0.25, 0.25*numel(Seq.tRep), numel(Seq.tRep));
    end
    if strcmp(Seq.AQPhaseOffset, 'same')
      Seq.AQPhaseOffset = Seq.TXPhaseOffset;
    end
    if strcmp(Seq.AQPhaseOffset, 'alternate')
      Seq.AQPhaseOffset = [0, 0+90*(-1).^(1:numel(Seq.tRep)-1)];
    end
  end


  AQ.Start = [max(pulseData1.Start(end) + pulseData1.Duration(end) + get_DeadTimeTX2RX(HW, AQ.fSample(1), SliceSelect.iDevice), GradTimeEnd), ...  % FID
    AQShift + Seq.tEcho/2 - Seq.tEcho*Seq.AQEcho/2*ones(1,Seq.nEchos)];  % start AQ window behind the dead time of TX pulses
  if Seq.AQFIDGrid, AQ.Start(1) = (ceil(AQ.Start(1).*AQ.fSample(1)-0.5)+0.5)./AQ.fSample(1); end
  AQ.nSamples = [round((Seq.tEcho/2*Seq.AQFID-AQ.Start(1))*AQ.fSample(1))*AQ1, ...  % FID
    round((Seq.tEcho*Seq.AQEcho)*AQ.fSample(end))*ones(1,Seq.nEchos)];  % calculate the number of samples
  AQ.ResetPhases = [1, zeros(1,Seq.nEchos)];  % reset the phase only once
  AQ.Phase = Seq.AQPhaseOffset;


  % AQFID relative part of the time tEcho/2 to be acquired between the first
  % pulse and the inversion pulse [0...1[
  % (0 => only one sample, negative no sample)
  if AQ.nSamples(1,1) < 1
    if Seq.AQFID == 0
      AQ.nSamples(1,1) = 1;
    else
      AQ.Start(1,1) = NaN;
      AQ.nSamples(1,1) = 0;
    end
  end

  % relative part of the time tEcho to be acquired between the inversion pulses [0...1[
  % (0 => only one sample)
  if size(AQ.nSamples,2) >= 2
    if any(AQ.nSamples(1,2:end)<1)
      AQ.nSamples(1,2:end) = 1;
      AQ.Start(1,2:end) = Seq.tEcho./2 - 0.5./AQ.fSample(1,2);
    end
  end

  if Seq.tEchoFirst~=Seq.tEcho && Seq.nEchos~=0
    AQ.Start(1,2) = AQ.Start(1,2) + Seq.tEchoFirst/2 - Seq.tEcho/2 + Seq.tEchoFirstTau2/2;
    if Seq.nEchos > 1
      AQ.Repeat(3) = 0;
    end
  end


  % optionally, execute prepare function (e.g. to add pre-polarization pulse)
  if ~isempty(Seq.Function_Prepare_Measurement)
    % temporarily set Seq.tEcho = 0
    SeqTemp = Seq;
    SeqTemp.tEcho = 0;
    SeqTemp.P90tReps = 1;
    [HW, SeqTemp, AQ, TX, Grad] = Seq.Function_Prepare_Measurement(HW, SeqTemp, AQ, TX, Grad);

    if isfield(SeqTemp, 'Prepol')
      Seq.Prepol = SeqTemp.Prepol;
    end
    if isfield(SeqTemp, 'DNP')
      Seq.DNP = SeqTemp.DNP;
    end
    if isfield(SeqTemp, 'DigitalIO')
     Seq.DigitalIO = SeqTemp.DigitalIO;
    end
  end


  if Seq.fast == 1
    % re-organize complete pulse program into one single tRep
    % add preview repetition times to start times and reshape to a n x 1 vector
    AQ.Start = (AQ.Start+([0,cumsum(Seq.tRep(1:end-1))])).';
    AQ.nSamples = AQ.nSamples.';
    AQ.fSample = AQ.fSample.';
    AQ.Phase = AQ.Phase.';
    AQ.ResetPhases = 1;
    AQ.Repeat = 0;

    % add preview repetition times to start times
    TX.Start = TX.Start + (ones(size(TX.Start,1),1)*([0,cumsum(Seq.tRep(1:end-1))]));

    % Reshape to a n x 1 vector
    TX.Start = TX.Start(:);
    TX.Duration = TX.Duration(:);
    TX.Frequency = TX.Frequency(:);
    TX.Phase = TX.Phase(:);
    TX.Amplitude = TX.Amplitude(:);
    TX.Repeat = 0;

    % Remove the NANs
    TX.Duration = TX.Duration(~isnan(TX.Start));
    TX.Frequency = TX.Frequency(~isnan(TX.Start));
    TX.Phase = TX.Phase(~isnan(TX.Start));
    TX.Amplitude = TX.Amplitude(~isnan(TX.Start));
    TX.Start = TX.Start(~isnan(TX.Start));

    Seq.tRep = sum(Seq.tRep)+100e-6;  % extend the repetition time
    if Seq.useSliceSelect
      for t = 1:HW.Grad(SliceSelect.iDevice).n
        Grad(t).Time = Grad(t).Time(:);
        Grad(t).Amp = Grad(t).Amp(:);
        Grad(t).Time = Grad(t).Time(~isnan(Grad(t).Time));
        Grad(t).Amp = Grad(t).Amp(~isnan(Grad(t).Time));
      end
    end
  else
    if isemptyfield(Seq, 'CLTime')
      if Seq.nEchos == 0, Seq.CLTime = 3e-6; end
      if Seq.nEchos == 1, Seq.CLTime = [3e-6, 1.448e-6]; end
      if Seq.nEchos >= 1
        if Seq.tEchoFirstTau2, Seq.CLTime(2) = 1.688e-6; end
      end
      if Seq.nEchos >= 2
        Seq.CLTime = [3e-6,1.448e-6,0400e-9+zeros(1,Seq.nEchos-1)];
        if Seq.tEchoFirst~=Seq.tEcho; Seq.CLTime(3) = 872e-9; end
      end

    else
      if numel(Seq.CLTime)==1, Seq.CLTime = Seq.CLTime+zeros(1,numel(Seq.tRep)); end
    end
  end

  if ~isempty(Seq.tsT2T2)
    TX.Duration = [TX.Duration, TX.Duration(:,end), TX.Duration(:,1), TX.Duration(:,1), repmat(TX.Duration(:,2),1,Seq.nEchosT2T2)];
    TX.Start = [TX.Start, TX.Start(:,end), TX.Start(:,1), TX.Start(:,1), repmat(TX.Start(:,2),1,Seq.nEchosT2T2)];
    TX.Amplitude = [TX.Amplitude, TX.Amplitude(:,end), TX.Amplitude(:,1), TX.Amplitude(:,1), repmat(TX.Amplitude(:,2),1,Seq.nEchosT2T2)];
    TX.Frequency = [TX.Frequency, TX.Frequency(:,end), TX.Frequency(:,1), TX.Frequency(:,1), repmat(TX.Frequency(:,2),1,Seq.nEchosT2T2)];
    TX.Phase = [TX.Phase, TX.Phase(:,end), TX.Phase(:,1)+180, TX.Phase(:,1), repmat(TX.Phase(:,2),1,Seq.nEchosT2T2)];
    TX.Repeat = [0, zeros(1,Seq.nEchos+3+Seq.nEchosT2T2)];

    AQ.Start = [AQ.Start, NaN(size(AQ.Start,1),1), NaN(size(AQ.Start,1),1), AQ.Start(:,1), repmat(AQ.Start(:,2),1,Seq.nEchosT2T2)];
    AQ.nSamples = [AQ.nSamples, NaN(size(AQ.nSamples,1),1), NaN(size(AQ.nSamples,1),1), AQ.nSamples(:,1), repmat(AQ.nSamples(:,2),1,Seq.nEchosT2T2)];
    % AQ.Frequency = [AQ.Frequency, NaN(size(AQ.Frequency,1),1), NaN(size(AQ.Frequency,1),1), AQ.Frequency(:,1), repmat(AQ.Frequency(:,2),1,Seq.nEchosT2T2)];
    AQ.fSample = [AQ.fSample, NaN(size(AQ.fSample,1),1), NaN(size(AQ.fSample,1),1), AQ.fSample(:,1), repmat(AQ.fSample(:,2),1,Seq.nEchosT2T2)];
    AQ.Phase = [AQ.Phase, NaN(size(AQ.Phase,1),1), NaN(size(AQ.Phase,1),1), AQ.Phase(:,1), repmat(AQ.Phase(:,2),1,Seq.nEchosT2T2)];
    AQ.ResetPhases = [1, zeros(1,Seq.nEchos+3+Seq.nEchosT2T2)];  % reset the phase only once
    AQ.Repeat = [0, zeros(1,Seq.nEchos+3+Seq.nEchosT2T2)];
    Seq.tRep = [Seq.tRep, Seq.tEcho/2, Seq.tsT2T2, Seq.tRep(1), repmat(Seq.tRep(2),1,Seq.nEchosT2T2)];
    Seq.CLTime = 5e-6 * [1, ones(1,Seq.nEchos+3+Seq.nEchosT2T2)];

    % FIXME: DampCoil
  end

  if ~isempty(Seq.SeqAverage)
    tempTXPhaseInvertIncrement = cumsum(repmat(Seq.SeqAverage.TXPhaseInvertIncrement,[1,Seq.nEchos,Seq.SeqAverage.average]),3)-Seq.SeqAverage.TXPhaseInvertIncrement;
    tempTXPhaseFlipIncrement = cumsum(repmat(Seq.SeqAverage.TXPhaseFlipIncrement,[1,1,Seq.SeqAverage.average]),3)-Seq.SeqAverage.TXPhaseFlipIncrement;
    tempAQPhaseIncrement = cumsum(repmat(Seq.SeqAverage.AQPhaseIncrement,[1,1+Seq.nEchos,Seq.SeqAverage.average]),3)-Seq.SeqAverage.AQPhaseIncrement(1);

    tempTXPhaseInvertOffset = repmat(permute(Seq.SeqAverage.TXPhaseInvertOffset,[1,3,2]),[1,Seq.nEchos,1]);
    tempTXPhaseFlipOffset = repmat(permute(Seq.SeqAverage.TXPhaseFlipOffset,[1,3,2]),[1,1,1]);
    tempAQPhaseOffset = repmat(permute(Seq.SeqAverage.TXPhaseFlipOffset,[1,3,2]),[1,1+Seq.nEchos,1]);

    Seq.tRep = repmat([Seq.tRep(1:end-1), Seq.tRep(end)+Seq.SeqAverage.averageBreak], 1, Seq.SeqAverage.average);
    Seq.tRep(end) = Seq.tRep(end)-Seq.SeqAverage.averageBreak;
    TX.Start = repmat(TX.Start, 1, Seq.SeqAverage.average);
    TX.Duration = repmat(TX.Duration, 1, Seq.SeqAverage.average);
    TX.Frequency = repmat(TX.Frequency, 1, Seq.SeqAverage.average);
    TX.Amplitude = repmat(TX.Amplitude, 1, Seq.SeqAverage.average);

    % Damp Coil
    if isfield(Seq, 'DigitalIO')
      if ~isemptyfield(Seq.DigitalIO, 'SetTime')
        Seq.DigitalIO.SetTime = repmat(Seq.DigitalIO.SetTime, 1, Seq.SeqAverage.average);
      end
      if ~isemptyfield(Seq.DigitalIO, 'SetValue')
        Seq.DigitalIO.SetValue = repmat(Seq.DigitalIO.SetValue, 1, Seq.SeqAverage.average);
      end
      if ~isemptyfield(Seq.DigitalIO, 'Repeat')
        Seq.DigitalIO.Repeat = repmat(Seq.DigitalIO.Repeat, 1, Seq.SeqAverage.average);
      end
    end

    % TX.Amplitude(1:end/2)=0;
    if isfield(TX, 'Repeat'), TX.Repeat = repmat(TX.Repeat, 1, Seq.SeqAverage.average); end

    AQ.Start = repmat(AQ.Start, 1, Seq.SeqAverage.average);
    AQ.fSample = repmat(AQ.fSample, 1, Seq.SeqAverage.average);

    if isempty(Seq.tsT2T2)
      TX.Phase = repmat(TX.Phase,1,Seq.SeqAverage.average) + repmat(reshape(cat(2,tempTXPhaseFlipIncrement+tempTXPhaseFlipOffset,tempTXPhaseInvertIncrement+tempTXPhaseInvertOffset),1,[]),size(TX.Phase,1),1);
      AQ.Phase = repmat(AQ.Phase,1,Seq.SeqAverage.average) + repmat(reshape(tempAQPhaseIncrement+tempAQPhaseOffset,1,[]),size(AQ.Start,1),1);
    else
      tempTXPhaseInvertIncrementT2T2=cumsum(repmat(Seq.SeqAverage.TXPhaseInvertIncrement,[1,Seq.nEchosT2T2,Seq.SeqAverage.average]),3)-Seq.SeqAverage.TXPhaseInvertIncrement;
      tempAQPhaseIncrementT2T2=cumsum(repmat(Seq.SeqAverage.AQPhaseIncrement,[1,1+Seq.nEchosT2T2,Seq.SeqAverage.average]),3)-Seq.SeqAverage.AQPhaseIncrement(1);

      tempTXPhaseInvertOffsetT2T2=repmat(permute(Seq.SeqAverage.TXPhaseInvertOffset,[1,3,2]),[1,Seq.nEchosT2T2,1]);
      tempAQPhaseOffsetT2T2=repmat(permute(Seq.SeqAverage.TXPhaseFlipOffset,[1,3,2]),[1,1+Seq.nEchosT2T2,1]);

      TX.Phase=repmat(TX.Phase,1,Seq.SeqAverage.average)+repmat(reshape(cat(2,tempTXPhaseFlipIncrement+tempTXPhaseFlipOffset,...
                                                                              tempTXPhaseInvertIncrement+tempTXPhaseInvertOffset,...
                                                                              tempTXPhaseFlipIncrement+tempTXPhaseFlipOffset,...
                                                                              tempTXPhaseFlipIncrement+tempTXPhaseFlipOffset,...
                                                                              tempTXPhaseFlipIncrement+tempTXPhaseFlipOffset,...
                                                                              tempTXPhaseInvertIncrementT2T2+tempTXPhaseInvertOffsetT2T2...
                                                                            ),1,[]),size(TX.Phase,1),1);
      AQ.Phase=repmat(AQ.Phase,1,Seq.SeqAverage.average)+repmat(reshape(cat(2,tempAQPhaseIncrement+tempAQPhaseOffset,...
                                                                              tempAQPhaseIncrement(:,1,:)+tempAQPhaseOffset(:,1,:),...
                                                                              tempAQPhaseIncrement(:,1,:)+tempAQPhaseOffset(:,1,:),...
                                                                              tempAQPhaseIncrementT2T2+tempAQPhaseOffsetT2T2...
                                                                            ),1,[]),size(AQ.Start,1),1);
    end
    AQ.ResetPhases = repmat(AQ.ResetPhases, 1, Seq.SeqAverage.average);
    %  AQ.ResetPhases(:)=1;
    AQ.nSamples = repmat(AQ.nSamples, 1, Seq.SeqAverage.average);
    if isfield(AQ,'Repeat'), AQ.Repeat = repmat(AQ.Repeat, 1, Seq.SeqAverage.average); end
    if isfield(Seq,'CLTime')
      if numel(Seq.CLTime)~=numel(Seq.tRep)
        Seq.CLTime=repmat(Seq.CLTime,1,Seq.SeqAverage.average);
      end
    end
    if size(AQ.Frequency,2) < size(AQ.Start,2);
      AQ.Frequency = repmat(AQ.Frequency(:,1), 1, size(AQ.Start,2)./Seq.SeqAverage.average);
    end
    AQ.Frequency = repmat(AQ.Frequency, 1, Seq.SeqAverage.average);
  end

  for t = 1:HW.Grad(SliceSelect.iDevice).n
    if numel(Grad(t).Time) ~= 1
      % if size(Grad(t).Time,2)<size(Seq.tRep,2); Grad(t).Time=repmat(AQ.Time(:,1),1,size(Seq.tRep,2));end;
      % if size(Grad(t).Amp,2)<size(Seq.tRep,2); Grad(t).Amp=repmat(AQ.Amp(:,1),1,size(Seq.tRep,2));end;
      % if size(Grad(t).Shim,2)<size(Seq.tRep,2); Grad(t).Shim=repmat(AQ.Shim(:,1),1,size(Seq.tRep,2));end;

      if ~isempty(Seq.SeqAverage)
        Grad(t).Time = repmat(Grad(t).Time, 1, Seq.SeqAverage.average);
        Grad(t).Amp = repmat(Grad(t).Amp, 1, Seq.SeqAverage.average);
        % Grad(t).Shim = repmat(Grad(t).Shim, 1, Seq.SeqAverage.average);
      end
    end
  end

  % prepare sequence
  StartSequence = Seq.StartSequence;
  PollPPGfast = Seq.PollPPGfast;
  GetRawData = Seq.GetRawData;
  PostProcessSequence = Seq.PostProcessSequence;
  Seq.StartSequence = 0;
  Seq.PollPPGfast = 0;
  Seq.GetRawData = 0;
  Seq.PostProcessSequence = 0;

  [~, Seq] = set_sequence(HW, Seq, AQ, TX, Grad);

  Seq.StartSequence = StartSequence;
  Seq.PollPPGfast = PollPPGfast;
  Seq.GetRawData = GetRawData;
  Seq.PostProcessSequence = PostProcessSequence;
  SeqOut = Seq;
  Seq.PreProcessSequence = 0;
else
  AQ = Seq.AQ;
  TX = Seq.TX;
  Grad = Seq.Grad;
end

data = struct();
data_1D = NaN;
if Seq.StartSequence || Seq.PollPPGfast || Seq.GetRawData || Seq.PostProcessSequence

  if ~Seq.RawData
    % run measurement
    [~, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);

    if ~isempty(Seq.SeqAverage)
      for iAQ = 1:numel(data)
        if isempty(data(iAQ).data) || (isscalar(data(iAQ).data) && isnan(data(iAQ).data))
          continue;
        end
        if Seq.SeqAverage.SaveSeqAverageData
          data(iAQ).SeqAverage.data = reshape(data(iAQ).data, ...
            [size(data(iAQ).data,1), size(SeqOut.AQ(iAQ).Start,1), size(data(iAQ).data,3)/Seq.SeqAverage.average, Seq.SeqAverage.average]);
        end
        data(iAQ).data = mean(reshape(data(iAQ).data, ...
          [size(data(iAQ).data,1), size(SeqOut.AQ(iAQ).Start,1), size(data(iAQ).data,3)/Seq.SeqAverage.average, Seq.SeqAverage.average]), 4);
        data(iAQ).time_all = reshape(data(iAQ).time_all, ...
          [size(data(iAQ).data,1), size(SeqOut.AQ(iAQ).Start,1), size(data(iAQ).data,3), Seq.SeqAverage.average]);
        data(iAQ).time_all = data(iAQ).time_all(:,:,:,1);
        data(iAQ).time_of_tRep = reshape(data(iAQ).time_of_tRep, ...
          [size(data(iAQ).data,1), size(SeqOut.AQ(iAQ).Start,1), size(data(iAQ).data,3), Seq.SeqAverage.average]);
        data(iAQ).time_of_tRep = data(iAQ).time_of_tRep(:,:,:,1);
      end
    end

    if Seq.PostProcessSequence
      % Plot
      if ishghandle(Seq.plot, 'figure') || ishghandle(Seq.plot, 'uipanel') || Seq.plot % && SeqOut.PostProcessSequence
        if isemptyfield(Seq, 'plotAllHandle'), Seq.plotAllHandle = Seq.plot; end
        SeqOut.plotAllHandle = plot_data_1D(HW, data_1D, Seq.plotAllHandle, Seq.plotRaiseWindow, Seq.plotData1D);
      end

      if Seq.plotTR % && SeqOut.PostProcessSequence
        plot_data_1D_TR(HW, data_1D);
      end
    end
  else
    % run measurement
    [RawData, SeqOut] = set_sequence(HW, Seq, AQ, TX, Grad);
    if ~isempty(SeqOut.SeqAverage)
      szRawData = [size(RawData,1), size(AQ.Start,1), SeqOut.nEchos, SeqOut.SeqAverage.average];
      if SeqOut.SeqAverage.SaveSeqAverageData
        data.SeqAverage.data = cat(3, nan([size(RawData,1), size(AQ.Start,1), 1, SeqOut.SeqAverage.average]), ...
          reshape(RawData, szRawData));
      end
      % FIXME: We should not average over the samples per Echo
      data.data = mean(mean(cat(3, nan([size(RawData,1), size(AQ.Start,1), 1, SeqOut.SeqAverage.average]), ...
        reshape(RawData, szRawData)), 4), 1);
    else
      data.data = permute([nan(size(RawData,1),1),RawData],[1 3 2]);
    end
    data.time_all = permute([0,linspace(SeqOut.tEcho,SeqOut.nEchos*SeqOut.tEcho,SeqOut.nEchos)], [1 3 2]);
    data.time_of_tRep = permute([0,SeqOut.tEcho+zeros(1,SeqOut.nEchos)], [1 3 2]);
    if Seq.nEchos>=2;
      if Seq.tEchoFirst~=Seq.tEcho;
        data.time_all(2:end)=data.time_all(2:end)+Seq.tEchoFirst-SeqOut.tEcho;
        data.time_of_tRep(2)=Seq.tEchoFirst+Seq.tEchoFirstTau2/2;
      end
    end
  end

end

end
