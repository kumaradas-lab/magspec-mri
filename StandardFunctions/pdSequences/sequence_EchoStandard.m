function [data, SeqOut, data_1D] = sequence_EchoStandard(HW, Seq, SliceSelectUser, temp)
%% Acquire FID, Spin Echo or CPMG Echo train
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
%             Pulse shape function for the preparation pulse before the actual
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
%     excitationFlipAngle
%             Flip angle of the excitation pulse in degrees (default: 90).
%
%     inversionFlipAngle
%             Flip angle of the inversion pulses in degrees (default: 180).
%
%     excitationPhase
%             Phase of the excitation pulse in degrees (default: 180).
%
%     inversionPhase
%             Phase of the inversion pulses in degrees (default: -90).
%
%     preparationPhase
%             Phase of the preparation pulse in degrees (default: -90).
%
%     Function_Prepare_Measurement
%             Function handle that is executed before the measurement is started
%             with the following function signature
%               [HW, Seq, AQ, TX, Grad] = @(HW, Seq, AQ, TX, Grad);
%             It can be used to modify the pulse sequence after the (default)
%             pulse sequence was created.
%
%     dualNuclear
%             Boolean value to switch on dual nucleus transmission and
%             acquisition (at HW.GammaDef and HW.GammaX). The connected hardware
%             and driver must be compatible for this mode. (Default: false)
%             If dualNuclear is set to true, the following fields can be used
%             (additionally to the ones without "X" suffix) for the HW.GammaX
%             pulses and acquisitions:
%             excitationFlipAngleX, inversionFlipAngleX, excitationPhaseX,
%             inversionPhaseX, preparationPhaseX, TXAmpP90X, p90X, p180X,
%             AQPhaseOffsetX, TXPhaseOffsetX, FlipPulseX, SlicePulseX,
%             InvertPulseX, FirstInvertPulseX
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
% (C) Copyright 2012-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

% FIXME: Complete documentation (including phase cycling using Seq.SeqAverage).

% Check if first argument is a talker.
% Please do not use talker as first argument.
if isa(HW, 'PD.Talker') || ...
    (isstruct(HW) && ~isfield(HW, 'MMRT'))
  warning('PD:sequence_EchoStandard:TalkerArgDeprecated', ...
    ['It looks like a talker object was passed as the first argument for ', ...
    'the function "%s".\n', ...
    'That syntax is deprecated and will be removed in a future version of OpenMatlab.'], ...
    mfilename());
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
  if isemptyfield(Seq, 'dualNuclear')
    % use mixed TX pulses and acquisitions
    Seq.dualNuclear = false;
  end
  if isemptyfield(Seq, 'excitationFlipAngle')
    Seq.excitationFlipAngle = 90;
  end
  if isemptyfield(Seq, 'inversionFlipAngle')
    Seq.inversionFlipAngle = 180;
  end
  if Seq.dualNuclear
    if isemptyfield(Seq, 'excitationFlipAngleX')
      Seq.excitationFlipAngleX = Seq.excitationFlipAngle;
    end
    if isemptyfield(Seq, 'inversionFlipAngleX')
      Seq.inversionFlipAngleX = Seq.inversionFlipAngle;
    end
  end
  if isemptyfield(Seq, 'excitationPhase')
    Seq.excitationPhase = 180;
  end
  if isemptyfield(Seq, 'inversionPhase')
    Seq.inversionPhase = -90;
  end
  if isemptyfield(Seq, 'preparationPhase')
    Seq.preparationPhase = -90;
  end
  if Seq.dualNuclear
    if isemptyfield(Seq, 'excitationPhaseX')
      Seq.excitationPhaseX = Seq.excitationPhase;
    end
    if isemptyfield(Seq, 'inversionPhaseX')
      Seq.inversionPhaseX = Seq.inversionPhase;
    end
    if isemptyfield(Seq, 'preparationPhaseX')
      Seq.preparationPhaseX = Seq.preparationPhase;
    end
  end
  if isemptyfield(Seq, 'TXAmp')
    % TX amplitude of inversion pulses in T
    if Seq.dualNuclear
      % Calculate amplitudes for both components such that sum equals
      % Def.PaUoutCalibrated and both pulses have the same duration.
      Seq.TXAmp = HW.TX(SliceSelect.iDevice).Def.PaUoutCalibrated(HW.TX(SliceSelect.iDevice).ChannelDef) / ...
        (1 / HW.TX(SliceSelect.iDevice).PaUout2Amplitude(HW.TX(SliceSelect.iDevice).ChannelDef) + ...
        HW.GammaDef / HW.GammaX / HW.TX(SliceSelect.iDevice).PaUout2AmplitudeX(HW.TX(SliceSelect.iDevice).ChannelDef));
      Seq.TXAmpX = Seq.TXAmp * HW.GammaDef / HW.GammaX;
      % FIXME: Would it be better to have pulse bandwidths that correspond to
      % the same slice thickness for a given slice gradient amplitude (gamma!!)?
      % Seq.TXAmp = HW.TX(SliceSelect.iDevice).Def.UoutCalibrated(HW.TX(SliceSelect.iDevice).ChannelDef) / ...
      %   (1 / HW.TX(SliceSelect.iDevice).PaUout2Amplitude(HW.TX(SliceSelect.iDevice).ChannelDef) + ...
      %    1 / HW.TX(SliceSelect.iDevice).PaUout2AmplitudeX(HW.TX(SliceSelect.iDevice).ChannelDef));
      % Seq.TXAmpX = Seq.TXAmp;
    else
      Seq.TXAmp = HW.TX(SliceSelect.iDevice).AmpDef;
    end
  end
  if isemptyfield(Seq, 'TXAmpP90')
    % TX amplitude of excitation pulse in T
    Seq.TXAmpP90 = Seq.TXAmp;
  end
  if isemptyfield(Seq, 'p90')
    % duration of 90 degrees pulse at Seq.TXAmpP90 in seconds
    Seq.p90 = pi/180 * Seq.excitationFlipAngle / HW.GammaDef / Seq.TXAmpP90;
  end
  if isemptyfield(Seq, 'p180')
    % duration of 180 degrees pulse at Seq.TXAmp in seconds
    Seq.p180 = pi/180 * Seq.inversionFlipAngle / HW.GammaDef / Seq.TXAmp;
  end
  if Seq.dualNuclear
    if isemptyfield(Seq, 'TXAmpP90X')
      % TX amplitude of excitation pulse at secondary frequency in T
      Seq.TXAmpP90X = Seq.TXAmpX;
    end
    if isemptyfield(Seq, 'p90X')
      Seq.p90X = pi/180 * Seq.excitationFlipAngleX / HW.GammaX / Seq.TXAmpP90X;
    end
    if isemptyfield(Seq, 'p180X')
      Seq.p180X = pi/180 * Seq.inversionFlipAngleX / HW.GammaX / Seq.TXAmpX;
    end
  end

  % acquisition time at echoes between the inversion pulses [0...x s (0 => only
  % one sample, negative no sample)
  Seq = set_EmptyField(Seq, 'tAQEcho', 0.8 * Seq.tEcho);
  % relative part of the time tEcho to be acquired between the inversions pulses
  % [0...1[ (0 => only one sample)
  Seq = set_EmptyField(Seq, 'AQEcho', Seq.tAQEcho/Seq.tEcho);
  % AQFID relative part of the time tEcho/2 to be acquired between the first
  % pulse and the inversion Pulse [0...1[ (0 => only one sample, negative no
  % sample)
  Seq = set_EmptyField(Seq, 'AQFID', Seq.AQEcho);
  Seq = set_EmptyField(Seq, 'AQFIDGrid', 0);            % AQFID start on 0s + fSample grid
  Seq = set_EmptyField(Seq, 'AQFrequency', HW.fLarmor); % mixing frequency of the receiver for all windows
  Seq = set_EmptyField(Seq, 'TXFrequency', HW.fLarmor); % frequency of the transmitter for all pulses
  Seq = set_EmptyField(Seq, 'AQPhaseOffset', 0);        % phase offset of the receiver for all windows
  Seq = set_EmptyField(Seq, 'TXPhaseOffset', 0);        % phase offset of the transmitter for all pulses
  if Seq.dualNuclear
    if isemptyfield(Seq, 'AQFrequencyX')
      Seq.AQFrequencyX = HW.fLarmorX;
    end
    if isemptyfield(Seq, 'TXFrequencyX')
      Seq.TXFrequencyX = HW.fLarmorX;
    end
    if Seq.dualNuclear
      if isemptyfield(Seq, 'AQPhaseOffsetX')
        Seq.AQPhaseOffsetX = Seq.AQPhaseOffset;
      end
      if isemptyfield(Seq, 'TXPhaseOffsetX')
        Seq.TXPhaseOffsetX = Seq.TXPhaseOffset;
      end
    end
  end
  Seq = set_EmptyField(Seq, 'tEchoTrain', Seq.tEcho*1);                     % duration of echo Train [s]
  Seq = set_EmptyField(Seq, 'nEchos', round(Seq.tEchoTrain/Seq.tEcho));     % Number of acquired echoes
  if ~isfield(Seq, 'tsT2T2'),           Seq.tsT2T2          = [];     end   % separation time between CPMG pulse trains
  Seq = set_EmptyField(Seq, 'tEchoTrainT2T2', Seq.nEchos*Seq.tEcho);        % duration of echo Train [s]
  Seq = set_EmptyField(Seq, 'nEchosT2T2', round(Seq.tEchoTrainT2T2/Seq.tEcho));  % number of Echoes of second CPMG pulse trains
  Seq = set_EmptyField(Seq, 'fSample', HW.RX(SliceSelect.iDevice).fSample/1250);  % Sampling rate of the AQ windows at the echoes
  Seq = set_EmptyField(Seq, 'fSampleFID', Seq.fSample);                     % Sampling rate of the AQ window at the FID
  if isemptyfield(Seq, 'SamplingFactor'), Seq.SamplingFactor = 1; end  % sampling factor for AQ windows at echoes
  if isemptyfield(Seq, 'SamplingFactorFID'), Seq.SamplingFactorFID = Seq.SamplingFactor; end  % sampling factor for AQ window at the FID
  Seq = set_EmptyField(Seq, 'Shim', zeros(1, HW.Grad(SliceSelect.iDevice).n));  % Shim add to HW.MagnetShim
  % set shim for additional gradient channels to zero
  Seq.Shim((numel(Seq.Shim)+1):HW.Grad(SliceSelect.iDevice).n) = 0;
  Seq = set_EmptyField(Seq, 'average', 1);  % Number of averages
  Seq = set_EmptyField(Seq, 'averageBreak', HW.FindFrequencyPause);  % Break between averages
  Seq = set_EmptyField(Seq, 'TXdelay', 0);  % Add delay to all TX pulses
  Seq = set_EmptyField(Seq, 'TXdelay90', 0);  % Add delay to the TX 90 degrees pulse
  if isemptyfield(Seq, 'FlipPulse')
    % 90 degrees pulse function handle for the excitation pulse
    Seq.FlipPulse =  @Pulse_Rect;
  end
  if isemptyfield(Seq, 'SlicePulse')
    % 90 degrees pulse function handle that is used if the measurement is slice
    % selective
    % FIXME: Why is this a separate property from Seq.FlipPulse
    Seq.SlicePulse =  Seq.FlipPulse;
  end
  if isemptyfield(Seq, 'InvertPulse')
    % 180 degrees pulse function handle for the inversion pulses
    Seq.InvertPulse =  @Pulse_Rect;
  end
  if isemptyfield(Seq, 'FirstInvertPulse')
    % 180 degrees pulse function handle for separated inversion pulse (preparation)
    Seq.FirstInvertPulse =  Seq.FlipPulse;
  end
  if Seq.dualNuclear
    if isemptyfield(Seq, 'FlipPulseX')
      % 90 degrees pulse function handle for the excitation pulse
      Seq.FlipPulseX =  @Pulse_Rect;
    end
    if isemptyfield(Seq, 'SlicePulseX')
      % 90 degrees pulse function handle that is used if the measurement is slice
      % selective
      % FIXME: Why is this a separate property from Seq.FlipPulse
      Seq.SlicePulseX =  Seq.FlipPulseX;
    end
    if isemptyfield(Seq, 'InvertPulseX')
      % 180 degrees pulse function handle for the inversion pulses
      Seq.InvertPulseX =  @Pulse_Rect;
    end
    if isemptyfield(Seq, 'FirstInvertPulse')
      % 180 degrees pulse function handle for separated inversion pulse (preparation)
      Seq.FirstInvertPulseX =  Seq.FlipPulseX;
    end
  end

  Seq.tFlip = Seq.p90 * Seq.FlipPulse(HW, 'Amp');  % Calculate the length needed for the 90 degrees pulse
  Seq.tInvertFirst = Seq.p90 * Seq.FirstInvertPulse(HW, 'Amp');  % Calculate the length needed for the 90 degrees pulse
  Seq.tInvert = Seq.p180 * Seq.InvertPulse(HW, 'Amp');  % Calculate the length needed for the 180 degrees pulse
  Seq = set_EmptyField(Seq, 'TxFlipBW', 1/Seq.tFlip* Seq.FlipPulse(HW,'Time'));   % Bandwidth of the RF pulses
  Seq = set_EmptyField(Seq, 'TxFirstInvertBW', 1/Seq.tInvertFirst* Seq.FirstInvertPulse(HW,'Time')); % Bandwidth of the RF pulses

  if Seq.dualNuclear
    % FIXME: Supporting differing pulse shape functions for X nucleus?
    % Calculate the length needed for the 90 degrees pulse
    Seq.tFlipX = Seq.p90X * Seq.FlipPulse(HW, 'Amp');
    % Calculate the length needed for the 90 degrees pulse
    Seq.tInvertFirstX = Seq.p90X * Seq.FirstInvertPulse(HW, 'Amp');
    % Calculate the length needed for the 180 degrees pulse
    Seq.tInvertX = Seq.p180X * Seq.InvertPulse(HW, 'Amp');
    if isemptyfield(Seq, 'TxFlipBWX')
      % Bandwidth of the RF pulses at X frequency
      Seq.TxFlipBWX = 1/Seq.tFlipX * Seq.FlipPulse(HW,'Time');
    end
    if isemptyfield(Seq, 'TxFirstInvertBWX')
      % Bandwidth of the inversion RF pulses at X frequeny
      Seq.TxFirstInvertBWX = 1/Seq.tInvertFirstX * Seq.FirstInvertPulse(HW,'Time');
    end
  end

  Seq = set_EmptyField(Seq, 'TxFlipSteps', 51);                             % Number of pulses used for the TX shape
  % If 1: The hole sequence will be put into the first Repetition time. Helps
  % recuding echo time.
  Seq = set_EmptyField(Seq, 'fast', 0);
%  Seq = set_EmptyField(Seq, 'UnblankTx2', 0);
  Seq = set_EmptyField(Seq, 'RawData', 0);
  if Seq.AQFID >= 0
    % Raw data mode doesn't work when also acquiring the FID.
    Seq.RawData = 0;
  end
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
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'average', 1);  % number of averages
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'averageBreak', HW.FindFrequencyPause); %
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'SaveSeqAverageData', 1);       %
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'TXPhaseInvertIncrement', 180); %
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'TXPhaseFlipIncrement', 0);     %
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'AQPhaseIncrement', 0);         %
    Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'RandomTXRXPhaseOffset', 0);    %
    if Seq.SeqAverage.RandomTXRXPhaseOffset
      [~, IX] = sort(rand(1, Seq.SeqAverage.average));
      Seq.SeqAverage.RandomTXRXPhaseOffset = ...
        linspace(0,360-360./Seq.SeqAverage.average,Seq.SeqAverage.average);
      Seq.SeqAverage.RandomTXRXPhaseOffset = Seq.SeqAverage.RandomTXRXPhaseOffset(IX);
      Seq.SeqAverage.TXPhaseFlipOffset = Seq.SeqAverage.RandomTXRXPhaseOffset;
      Seq.SeqAverage.TXPhaseInvertOffset = Seq.SeqAverage.RandomTXRXPhaseOffset;
      Seq.SeqAverage.AQPhaseOffset = Seq.SeqAverage.RandomTXRXPhaseOffset;
    else
      Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'TXPhaseInvertOffset', 0);            %
      Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'TXPhaseFlipOffset', 0);              %
      Seq.SeqAverage = set_EmptyField(Seq.SeqAverage, 'AQPhaseOffset', 0);                  %
      if numel(Seq.SeqAverage.TXPhaseInvertOffset) < Seq.SeqAverage.average
        Seq.SeqAverage.TXPhaseInvertOffset = repmat(Seq.SeqAverage.TXPhaseInvertOffset, 1, Seq.SeqAverage.average);
      end
      if numel(Seq.SeqAverage.TXPhaseFlipOffset) < Seq.SeqAverage.average
        Seq.SeqAverage.TXPhaseFlipOffset = repmat(Seq.SeqAverage.TXPhaseFlipOffset, 1, Seq.SeqAverage.average);
      end
      if numel(Seq.SeqAverage.AQPhaseOffset) < Seq.SeqAverage.average
        Seq.SeqAverage.AQPhaseOffset = repmat(Seq.SeqAverage.AQPhaseOffset, 1, Seq.SeqAverage.average);
      end
    end
  end

  if ~isfield(Seq, 'plotData1D'), Seq.plotData1D = []; end

  Seq.tInvertmax = 1000e-6;
  Seq.tAQmax = Seq.tEcho/2;
  % slice selective excitation pulse
  if Seq.dualNuclear
    % primary gamma
    Seq.tTxSlicemax = ...
      max(HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec / ...
          (HW.TX(SliceSelect.iDevice).AmpDef/(1 + HW.GammaDef / HW.GammaX)*0.9) * (pi/2)/pi * ...
          Seq.SlicePulse(HW, 'Amp'), ...
          1/(SliceSelect.MaxGradAmpSlice*HW.GammaDef/2/pi*SliceSelect.thickness) * ...
          Seq.SlicePulse(HW, 'Time'));
    % secondary gamma
    Seq.tTxSlicemax = max(Seq.tTxSlicemax, ...
      max(pi/HW.GammaX / ...
          (HW.TX(SliceSelect.iDevice).AmpDef/(1 + HW.GammaX / HW.GammaDef)*0.9) * (pi/2)/pi * ...
          Seq.SlicePulse(HW, 'Amp'), ...
          1/(SliceSelect.MaxGradAmpSlice*HW.GammaX/2/pi*SliceSelect.thickness) * ...
          Seq.SlicePulse(HW, 'Time')));
  else
    Seq.tTxSlicemax = max(HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/(HW.TX(SliceSelect.iDevice).AmpDef*0.9)*(pi/2)/pi * Seq.SlicePulse(HW, 'Amp'), ...
      1/(SliceSelect.MaxGradAmpSlice*HW.GammaDef/2/pi*SliceSelect.thickness) * Seq.SlicePulse(HW, 'Time'));
  end
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
  AQ.SamplingFactor = [Seq.SamplingFactorFID, Seq.SamplingFactor*ones(1,Seq.nEchos)];
  if Seq.nEchos == 1 && Seq.AQEcho >= 1
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
  if Seq.nEchos == 0
    % special case: only the FID is generated, the AQ window can be longer
    if HW.RX(SliceSelect.iDevice).ClampCoil.Enable
      clampPostset = HW.RX(SliceSelect.iDevice).ClampCoil.nSamplePostset/AQ.fSample(1) + ...
        HW.RX(SliceSelect.iDevice).ClampCoil.tPostset;
    else
      clampPostset = 0;
    end
    % adapt the repetition time to acquisition window end
    Seq.tRep = ...
      max(Seq.tEcho/2*Seq.AQFID ...
          + max(get_DeadTimeRX2TX(HW, AQ.fSample(1), SliceSelect.iDevice) , clampPostset) ...
          + 500e-6, ...  % FIXME: These additional 500 us seem arbitrary. Is there an explanation for this gap?
          Seq.tRepMin);
  else
    Seq.tRep(1) = Seq.tRep(1) + Seq.tEchoFirst/2 - Seq.tEcho/2 - Seq.tEchoFirstTau2/2;
    Seq.tRep(2) = Seq.tRep(2) + Seq.tEchoFirst/2 - Seq.tEcho/2 + Seq.tEchoFirstTau2/2;
  end

  if Seq.useSliceSelect
    if any(Seq.EndSpoilAmp)
      Seq.tRep(end) = max(Seq.tRep(end), ...
        EndSpoilGradTime(end) - GradTimeSart + max(HW.Grad(SliceSelect.iDevice).TimeDelay) + 1e-3);
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
  idxTX_common = 1:(Seq.dualNuclear+1);
  if isscalar(Seq.TXPhaseOffset)
    Seq.TXPhaseOffset = repmat(Seq.TXPhaseOffset, 1, numel(Seq.tRep));
    TX(1).Repeat = [0, zeros(1,double(Seq.nEchos>0)), ones(1,max(0,Seq.nEchos-1))];
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

  if Seq.dualNuclear
    if isscalar(Seq.TXPhaseOffsetX)
      Seq.TXPhaseOffsetX = repmat(Seq.TXPhaseOffsetX, 1, numel(Seq.tRep));
      TX(2).Repeat = [0, zeros(1,double(Seq.nEchos>0)), ones(1,max(0,Seq.nEchos-1))];
    end
    if ischar(Seq.TXPhaseOffsetX)
      if strcmp(Seq.TXPhaseOffsetX, 'rand')
        Seq.TXPhaseOffsetX = 360 * rand(1, numel(Seq.tRep));
      end
      if strcmp(Seq.TXPhaseOffsetX, 'linear')
        Seq.TXPhaseOffsetX = 360 * linspace(0, 1, numel(Seq.tRep));
      end
      if strcmp(Seq.TXPhaseOffsetX, 'alternate')
        Seq.TXPhaseOffsetX = [0, 90+90*(-1).^(1:numel(Seq.tRep)-1)];
      end
    end
  end

  if Seq.useSliceSelect
    % get properties for 90 degrees slice selective excitation pulse
    pulseData1 = Seq.SlicePulse(HW, Seq.TXdelay+Seq.TXdelay90, Seq.TxSliceBW, ...
      pi*Seq.p90/(HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/Seq.TXAmpP90), Seq.TxFlipSteps, ...
      Seq.tTxSlicemax, Seq.TXFrequency, Seq.excitationPhase);
  else
    % get properties for 90 degrees excitation pulse
    pulseData1 = Seq.FlipPulse(HW, Seq.TXdelay+Seq.TXdelay90, Seq.TxFlipBW, ...
      pi*Seq.p90/(HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/Seq.TXAmpP90), Seq.TxFlipSteps, ...
      Seq.tFlip, Seq.TXFrequency, Seq.excitationPhase);
  end
  % get properties for 180 degrees inversion pulses
  pulseData2 = Seq.InvertPulse(HW, Seq.TXdelay, 1/Seq.tInvert*Seq.InvertPulse(HW,'Time'), ...
    pi*Seq.p180/(HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/Seq.TXAmp), 20, ...
    Seq.tInvert, Seq.TXFrequency, Seq.inversionPhase);

  if Seq.tEchoFirstTau2
    % get properties for second 90 degrees pulse (first excitation [T1])
    pulseData3 = Seq.FirstInvertPulse(HW, Seq.TXdelay, Seq.TxFirstInvertBW, ...
      pi*Seq.p90/(HW.TX(SliceSelect.iDevice).Amp2FlipPiIn1Sec/Seq.TXAmpP90), ...
      Seq.TxFlipSteps, Seq.tInvertFirst, Seq.TXFrequency, Seq.preparationPhase);
    TX(1).Start = NaN(max([size(pulseData1.Start,1), size(pulseData2.Start,1), 2*size(pulseData3.Start,1)]), 1+Seq.nEchos);
  else
    TX(1).Start = NaN(max(size(pulseData1.Start,1), size(pulseData2.Start,1)), 1+Seq.nEchos);
  end

  if Seq.dualNuclear
    if Seq.useSliceSelect
      % get properties for 90 degrees slice selective excitation pulse
      pulseData1X = Seq.SlicePulseX(HW, Seq.TXdelay+Seq.TXdelay90, Seq.TxSliceBWX, ...
        pi*Seq.p90X/(pi/HW.GammaX/Seq.TXAmpP90X), Seq.TxFlipSteps, ...
        Seq.tTxSlicemax, Seq.TXFrequencyX, Seq.excitationPhaseX);
    else
      % get properties for 90 degrees excitation pulse
      pulseData1X = Seq.FlipPulseX(HW, Seq.TXdelay+Seq.TXdelay90, Seq.TxFlipBWX, ...
        pi*Seq.p90X/(pi/HW.GammaX/Seq.TXAmpP90X), Seq.TxFlipSteps, ...
        Seq.tFlipX, Seq.TXFrequencyX, Seq.excitationPhaseX);
    end
    % get properties for 180 degrees inversion pulses
    pulseData2X = Seq.InvertPulseX(HW, Seq.TXdelay, 1/Seq.tInvertX*Seq.InvertPulse(HW,'Time'), ...
      pi*Seq.p180X/(pi/HW.GammaX/Seq.TXAmpX), 20, ...
      Seq.tInvertX, Seq.TXFrequencyX, Seq.inversionPhaseX);

    if Seq.tEchoFirstTau2
      % get properties for second 90 degrees pulse (first excitation [T1])
      pulseData3X = Seq.FirstInvertPulseX(HW, Seq.TXdelay, Seq.TxFirstInvertBW, ...
        pi*Seq.p90X/(pi/HW.GammaX/Seq.TXAmpP90X), Seq.TxFlipSteps, ...
        Seq.tInvertFirstX, Seq.TXFrequencyX, Seq.preparationPhaseX);
      TX(2).Start = NaN(max([size(pulseData1X.Start,1), size(pulseData2X.Start,1), 2*size(pulseData3X.Start,1)]), 1+Seq.nEchos);
    else
      TX(2).Start = NaN(max(size(pulseData1X.Start,1), size(pulseData2X.Start,1)), 1+Seq.nEchos);
    end
  end

  % Prepare matrices
  [TX(idxTX_common).Duration] = deal(TX.Start);
  [TX(idxTX_common).Amplitude] = deal(TX.Start);
  [TX(idxTX_common).Frequency] = deal(TX.Start);
  [TX(idxTX_common).Phase] = deal(TX.Start);

  % Insert data of 90 degrees pulse in matrices
  TX(1).Duration(1:size(pulseData1.Start,1),1) = pulseData1.Duration;
  TX(1).Start(1:size(pulseData1.Start,1),1) = pulseData1.Start;
  TX(1).Amplitude(1:size(pulseData1.Start,1),1) = pulseData1.Amplitude;
  TX(1).Frequency(1:size(pulseData1.Start,1),1) = pulseData1.Frequency;
  TX(1).Phase(1:size(pulseData1.Start,1),1) = pulseData1.Phase + Seq.TXPhaseOffset(1);

  % Insert data of 180 degrees pulse in matrices
  if ~Seq.tEchoFirstTau2
    if Seq.nEchos >= 1
      TX(1).Duration(1:size(pulseData2.Start,1),2:end) = pulseData2.Duration * ones(1,Seq.nEchos);
      TX(1).Start(1:size(pulseData2.Start,1),2:end) = pulseData2.Start * ones(1,Seq.nEchos);
      TX(1).Amplitude(1:size(pulseData2.Start,1),2:end) = pulseData2.Amplitude * ones(1,Seq.nEchos);
      TX(1).Frequency(1:size(pulseData2.Start,1),2:end) = pulseData2.Frequency * ones(1,Seq.nEchos);
      TX(1).Phase(1:size(pulseData2.Start,1),2:end) = ...
        pulseData2.Phase*ones(1,Seq.nEchos) + ...
        repmat(Seq.TXPhaseOffset(2:end), size(pulseData2.Phase,1), 1);
    end
  else
    if Seq.nEchos >= 1
      TX(1).Duration(1:size(pulseData3.Start,1),2) = pulseData3.Duration;
      TX(1).Start(1:size(pulseData3.Start,1),2) = pulseData3.Start;
      TX(1).Amplitude(1:size(pulseData3.Start,1),2) = pulseData3.Amplitude;
      TX(1).Frequency(1:size(pulseData3.Start,1),2) = pulseData3.Frequency;
      TX(1).Phase(1:size(pulseData3.Start,1),2) = ...
        pulseData3.Phase + repmat(Seq.TXPhaseOffset(2), size(pulseData3.Phase,1), 1);

      TX(1).Duration(size(pulseData3.Start,1)+1:2*size(pulseData3.Start,1),2) = ...
        pulseData3.Duration;
      TX(1).Start(size(pulseData3.Start,1)+1:2*size(pulseData3.Start,1),2) = ...
        pulseData3.Start+Seq.tEchoFirstTau2;
      TX(1).Amplitude(size(pulseData3.Start,1)+1:2*size(pulseData3.Start,1),2) = ...
        pulseData3.Amplitude;
      TX(1).Frequency(size(pulseData3.Start,1)+1:2*size(pulseData3.Start,1),2) = ...
        pulseData3.Frequency;
      TX(1).Phase(size(pulseData3.Start,1)+1:2*size(pulseData3.Start,1),2) = ...
        pulseData3.Phase + repmat(Seq.TXPhaseOffset(2), size(pulseData3.Phase,1), 1);
    end
    if Seq.nEchos > 1
      TX(1).Duration(1:size(pulseData2.Start,1),3:end) = pulseData2.Duration * ones(1, Seq.nEchos-1);
      TX(1).Start(1:size(pulseData2.Start,1),3:end) = pulseData2.Start * ones(1,Seq.nEchos-1);
      TX(1).Amplitude(1:size(pulseData2.Start,1),3:end) = pulseData2.Amplitude * ones(1,Seq.nEchos-1);
      TX(1).Frequency(1:size(pulseData2.Start,1),3:end) = pulseData2.Frequency * ones(1,Seq.nEchos-1);
      TX(1).Phase(1:size(pulseData2.Start,1),3:end) = ...
        pulseData2.Phase*ones(1,Seq.nEchos-1) + ...
        repmat(Seq.TXPhaseOffset(3:end), size(pulseData2.Phase,1), 1);
    end
  end

  if Seq.dualNuclear
    % Insert data of 90 degrees pulse at X frequency in matrices
    TX(2).Duration(1:size(pulseData1X.Start,1),1) = pulseData1X.Duration;
    TX(2).Start(1:size(pulseData1X.Start,1),1) = pulseData1X.Start;
    TX(2).Amplitude(1:size(pulseData1X.Start,1),1) = pulseData1X.Amplitude;
    TX(2).Frequency(1:size(pulseData1X.Start,1),1) = pulseData1X.Frequency;
    TX(2).Phase(1:size(pulseData1X.Start,1),1) = pulseData1X.Phase + Seq.TXPhaseOffsetX(1);

    % Insert data of 180 degrees pulse in matrices
    if ~Seq.tEchoFirstTau2
      if Seq.nEchos >= 1
        TX(2).Duration(1:size(pulseData2X.Start,1),2:end) = pulseData2X.Duration * ones(1,Seq.nEchos);
        TX(2).Start(1:size(pulseData2X.Start,1),2:end) = pulseData2X.Start * ones(1,Seq.nEchos);
        TX(2).Amplitude(1:size(pulseData2X.Start,1),2:end) = pulseData2X.Amplitude * ones(1,Seq.nEchos);
        TX(2).Frequency(1:size(pulseData2X.Start,1),2:end) = pulseData2X.Frequency * ones(1,Seq.nEchos);
        TX(2).Phase(1:size(pulseData2X.Start,1),2:end) = ...
          pulseData2X.Phase * ones(1,Seq.nEchos) + ...
          repmat(Seq.TXPhaseOffsetX(2:end),size(pulseData2X.Phase,1),1);
      end
    else
      if Seq.nEchos >= 1
        TX(2).Duration(1:size(pulseData3X.Start,1),2) = pulseData3X.Duration;
        TX(2).Start(1:size(pulseData3X.Start,1),2) = pulseData3X.Start;
        TX(2).Amplitude(1:size(pulseData3X.Start,1),2) = pulseData3X.Amplitude;
        TX(2).Frequency(1:size(pulseData3X.Start,1),2) = pulseData3X.Frequency;
        TX(2).Phase(1:size(pulseData3X.Start,1),2) = ...
          pulseData3X.Phase + repmat(Seq.TXPhaseOffsetX(2), size(pulseData3X.Phase,1), 1);

        TX(2).Duration(size(pulseData3X.Start,1)+1:2*size(pulseData3X.Start,1),2) = ...
          pulseData3X.Duration;
        TX(2).Start(size(pulseData3X.Start,1)+1:2*size(pulseData3X.Start,1),2) = ...
          pulseData3X.Start + Seq.tEchoFirstTau2;
        TX(2).Amplitude(size(pulseData3X.Start,1)+1:2*size(pulseData3X.Start,1),2) = ...
          pulseData3X.Amplitude;
        TX(2).Frequency(size(pulseData3X.Start,1)+1:2*size(pulseData3X.Start,1),2) = ...
          pulseData3X.Frequency;
        TX(2).Phase(size(pulseData3X.Start,1)+1:2*size(pulseData3X.Start,1),2) = ...
          pulseData3X.Phase + repmat(Seq.TXPhaseOffsetX(2), size(pulseData3X.Phase,1), 1);
      end
      if Seq.nEchos > 1
        TX(2).Duration(1:size(pulseData2X.Start,1),3:end) = pulseData2X.Duration * ones(1,Seq.nEchos-1);
        TX(2).Start(1:size(pulseData2X.Start,1),3:end) = pulseData2X.Start * ones(1,Seq.nEchos-1);
        TX(2).Amplitude(1:size(pulseData2X.Start,1),3:end) = pulseData2X.Amplitude * ones(1,Seq.nEchos-1);
        TX(2).Frequency(1:size(pulseData2X.Start,1),3:end) = pulseData2X.Frequency * ones(1,Seq.nEchos-1);
        TX(2).Phase(1:size(pulseData2X.Start,1),3:end) = ...
          pulseData2X.Phase * ones(1,Seq.nEchos-1) + ...
          repmat(Seq.TXPhaseOffsetX(3:end), size(pulseData2X.Phase,1), 1);
      end
    end
  end


  % Damp Coil
  if Seq.DampCoil.Enable
    if HW.TX(SliceSelect.iDevice).DampCoil.Enable
      warning('PD:sequence_EchoStandard:DoubleDamp', ...
        'Damping is activated in HW.TX.DampCoil. Are you sure you want to include additional damping on the sequence level?');
    end
    if ~isemptyfield(Seq, 'DigitalIO')
      warning('PD:sequence_EchoStandard:DampCoil', ...
        'Damping Coil not supported if Seq.DigitalIO is already defined');
      Seq.DampCoil.Enable = false;
    elseif Seq.dualNuclear
      warning('PD:sequence_EchoStandard:DampCoil', ...
        'Damping Coil in sequence level not supported for dual nuclear measurements. Use HW.TX.DampCoil.Enable instead.');
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
  if isscalar(Seq.AQPhaseOffset)
    Seq.AQPhaseOffset = repmat(Seq.AQPhaseOffset, 1, numel(Seq.tRep));
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
  if Seq.dualNuclear
    if isscalar(Seq.AQPhaseOffsetX)
      Seq.AQPhaseOffsetX = repmat(Seq.AQPhaseOffsetX, 1, numel(Seq.tRep));
    else
      AQ.Repeat = [0, zeros(1,max(0,Seq.nEchos))];
    end
    if ischar(Seq.AQPhaseOffsetX)
      if strcmp(Seq.AQPhaseOffsetX, 'rand')
        Seq.AQPhaseOffsetX = 360 * rand(1, numel(Seq.tRep));
      end
      if strcmp(Seq.AQPhaseOffsetX, 'linear')
        Seq.AQPhaseOffsetX = -360*linspace(0.25, 0.25*numel(Seq.tRep), numel(Seq.tRep));
      end
      if strcmp(Seq.AQPhaseOffsetX, 'same')
        Seq.AQPhaseOffsetX = Seq.TXPhaseOffsetX;
      end
      if strcmp(Seq.AQPhaseOffsetX, 'alternate')
        Seq.AQPhaseOffsetX = [0, 0+90*(-1).^(1:numel(Seq.tRep)-1)];
      end
    end
  end


  endPulse1 = pulseData1.Start(end) + pulseData1.Duration(end);
  if Seq.dualNuclear
    endPulse1 = max(endPulse1, pulseData1X.Start(end) + pulseData1X.Duration(end));
  end
  AQ.Start = [max(endPulse1 + get_DeadTimeTX2RX(HW, AQ.fSample(1), SliceSelect.iDevice), GradTimeEnd), ...  % FID
    AQShift + Seq.tEcho/2 - Seq.tEcho*Seq.AQEcho/2*ones(1,Seq.nEchos)];  % start AQ window behind the dead time of TX pulses
  if Seq.AQFIDGrid, AQ.Start(1) = (ceil(AQ.Start(1).*AQ.fSample(1)-0.5)+0.5)./AQ.fSample(1); end
  AQ.nSamples = [round((Seq.tEcho/2*Seq.AQFID-AQ.Start(1))*AQ.fSample(1))*AQ1, ...  % FID
    round((Seq.tEcho*Seq.AQEcho)*AQ.fSample(end))*ones(1,Seq.nEchos)];  % calculate the number of samples
  AQ.ResetPhases = [1, zeros(1,Seq.nEchos)];  % reset the phase only once
  AQ.Phase = Seq.AQPhaseOffset;
  if Seq.dualNuclear
    AQ.PhaseX = Seq.AQPhaseOffsetX;
  end

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

  if HW.RX(SliceSelect.iDevice).ClampCoil.Enable
    % extent last tRep for clamp coil signal (if it contains AQ window)
    lastAQ = find(~isnan(AQ.Start(:,end)), 1, 'last');
    Seq.tRep(end) = max(Seq.tRep(end), ...
      AQ.Start(lastAQ,end) + AQ.nSamples(lastAQ,end)/AQ.fSample(lastAQ,end) + HW.RX(SliceSelect.iDevice).ClampCoil.tPostset + 0.1e-3);
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
    AQ.SamplingFactor = AQ.SamplingFactor.';
    AQ.ResetPhases = 1;
    AQ.Repeat = 0;
    if Seq.dualNuclear
      AQ.PhaseX = AQ.PhaseX.';
    end


    for iTX = idxTX_common
      % add preview repetition times to start times
      TX(iTX).Start = TX(iTX).Start + (ones(size(TX(iTX).Start,1),1)*([0,cumsum(Seq.tRep(1:end-1))]));

      % Reshape to a n x 1 vector
      TX(iTX).Start = TX(iTX).Start(:);
      TX(iTX).Duration = TX(iTX).Duration(:);
      TX(iTX).Frequency = TX(iTX).Frequency(:);
      TX(iTX).Phase = TX(iTX).Phase(:);
      TX(iTX).Amplitude = TX(iTX).Amplitude(:);
      TX(iTX).Repeat = 0;

      % Remove the NANs
      TX(iTX).Duration = TX(iTX).Duration(~isnan(TX(iTX).Start));
      TX(iTX).Frequency = TX(iTX).Frequency(~isnan(TX(iTX).Start));
      TX(iTX).Phase = TX(iTX).Phase(~isnan(TX(iTX).Start));
      TX(iTX).Amplitude = TX(iTX).Amplitude(~isnan(TX(iTX).Start));
      TX(iTX).Start = TX(iTX).Start(~isnan(TX(iTX).Start));
    end

    % Grad (might be set by prepare function or Seq.useSliceSelect)
    for iGrad = 1:numel(Grad)
      if ~isemptyfield(Grad(iGrad), 'Time')
        Grad(iGrad).Time = Grad(iGrad).Time + (ones(size(Grad(iGrad).Time,1),1)*([0,cumsum(Seq.tRep(1:end-1))]));
        Grad(iGrad).Time = Grad(iGrad).Time(:);
        Grad(iGrad).Amp = Grad(iGrad).Amp(:);
        Grad(iGrad).Amp = Grad(iGrad).Amp(~isnan(Grad(iGrad).Time));
        Grad(iGrad).Time = Grad(iGrad).Time(~isnan(Grad(iGrad).Time));
      end
    end

    % Digital IO (might be set by prepare function)
    if ~isemptyfield(Seq, 'DigitalIO')
      for iIO = 1:numel(Seq.DigitalIO)
        if ~isemptyfield(Seq.DigitalIO(iIO), 'SetTime')
          % append NaN values for trailing tReps without digital output changes
          if size(Seq.DigitalIO(iIO).SetTime, 2) < numel(Seq.tRep)
            Seq.DigitalIO(iIO).SetTime(:,(size(Seq.DigitalIO(iIO).SetTime, 2)+1):numel(Seq.tRep)) = NaN;
          end
          if size(Seq.DigitalIO(iIO).SetValue, 2) < numel(Seq.tRep)
            Seq.DigitalIO(iIO).SetValue(:,(size(Seq.DigitalIO(iIO).SetValue, 2)+1):numel(Seq.tRep)) = NaN;
          end
          Seq.DigitalIO(iIO).SetTime = Seq.DigitalIO(iIO).SetTime + (ones(size(Seq.DigitalIO(iIO).SetTime,1),1)*([0,cumsum(Seq.tRep(1:end-1))]));
          Seq.DigitalIO(iIO).SetTime = Seq.DigitalIO(iIO).SetTime(:);
          Seq.DigitalIO(iIO).SetValue = Seq.DigitalIO(iIO).SetValue(:);
          Seq.DigitalIO(iIO).SetValue = Seq.DigitalIO(iIO).SetValue(~isnan(Seq.DigitalIO(iIO).SetTime));
          Seq.DigitalIO(iIO).SetTime = Seq.DigitalIO(iIO).SetTime(~isnan(Seq.DigitalIO(iIO).SetTime));
        end
      end
    end

    Seq.tRep = sum(Seq.tRep)+100e-6;  % extend the repetition time
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
    for iTX = idxTX_common
      TX(iTX).Duration = ...
        [TX(iTX).Duration, TX(iTX).Duration(:,end), TX(iTX).Duration(:,1), ...
         TX(iTX).Duration(:,1), repmat(TX(iTX).Duration(:,2),1,Seq.nEchosT2T2)];
      TX(iTX).Start = ...
        [TX(iTX).Start, TX(iTX).Start(:,end), TX(iTX).Start(:,1), ...
         TX(iTX).Start(:,1), repmat(TX(iTX).Start(:,2),1,Seq.nEchosT2T2)];
      TX(iTX).Amplitude = ...
        [TX(iTX).Amplitude, TX(iTX).Amplitude(:,end), TX(iTX).Amplitude(:,1), ...
         TX(iTX).Amplitude(:,1), repmat(TX(iTX).Amplitude(:,2),1,Seq.nEchosT2T2)];
      TX(iTX).Frequency = ...
        [TX(iTX).Frequency, TX(iTX).Frequency(:,end), TX(iTX).Frequency(:,1), ...
        TX(iTX).Frequency(:,1), repmat(TX(iTX).Frequency(:,2),1,Seq.nEchosT2T2)];
      TX(iTX).Phase = ...
        [TX(iTX).Phase, TX(iTX).Phase(:,end), TX(iTX).Phase(:,1)+180, ...
         TX(iTX).Phase(:,1), repmat(TX(iTX).Phase(:,2),1,Seq.nEchosT2T2)];
      TX(iTX).Repeat = [0, zeros(1,Seq.nEchos+3+Seq.nEchosT2T2)];
    end

    AQ.Start = ...
      [AQ.Start, NaN(size(AQ.Start,1),1), NaN(size(AQ.Start,1),1), ...
       AQ.Start(:,1), repmat(AQ.Start(:,2),1,Seq.nEchosT2T2)];
    AQ.nSamples = ...
      [AQ.nSamples, NaN(size(AQ.nSamples,1),1), NaN(size(AQ.nSamples,1),1), ...
       AQ.nSamples(:,1), repmat(AQ.nSamples(:,2),1,Seq.nEchosT2T2)];
    % AQ.Frequency = ...
    %   [AQ.Frequency, NaN(size(AQ.Frequency,1),1), NaN(size(AQ.Frequency,1),1), ...
    %    AQ.Frequency(:,1), repmat(AQ.Frequency(:,2),1,Seq.nEchosT2T2)];
    AQ.fSample = ...
      [AQ.fSample, NaN(size(AQ.fSample,1),1), NaN(size(AQ.fSample,1),1), ...
       AQ.fSample(:,1), repmat(AQ.fSample(:,2),1,Seq.nEchosT2T2)];
    AQ.SamplingFactor = ...
      [AQ.SamplingFactor, NaN(size(AQ.SamplingFactor,1),1), NaN(size(AQ.SamplingFactor,1),1), ...
       AQ.SamplingFactor(:,1), repmat(AQ.SamplingFactor(:,2),1,Seq.nEchosT2T2)];
    AQ.Phase = ...
      [AQ.Phase, NaN(size(AQ.Phase,1),1), NaN(size(AQ.Phase,1),1), ...
       AQ.Phase(:,1), repmat(AQ.Phase(:,2),1,Seq.nEchosT2T2)];
    AQ.ResetPhases = [1, zeros(1,Seq.nEchos+3+Seq.nEchosT2T2)];  % reset the phase only once
    AQ.Repeat = [0, zeros(1,Seq.nEchos+3+Seq.nEchosT2T2)];
    if Seq.dualNuclear
      AQ.PhaseX = ...
        [AQ.PhaseX, NaN(size(AQ.PhaseX,1),1), NaN(size(AQ.PhaseX,1),1), ...
         AQ.PhaseX(:,1), repmat(AQ.PhaseX(:,2),1,Seq.nEchosT2T2)];
    end
    Seq.tRep = [Seq.tRep, Seq.tEcho/2, Seq.tsT2T2, Seq.tRep(1), repmat(Seq.tRep(2),1,Seq.nEchosT2T2)];
    Seq.CLTime = 5e-6 * [1, ones(1,Seq.nEchos+3+Seq.nEchosT2T2)];

    % FIXME: DampCoil
  end

  if ~isempty(Seq.SeqAverage)
    if Seq.fast
      tempTXPhaseInvertIncrement = ...
        cumsum(repmat(Seq.SeqAverage.TXPhaseInvertIncrement, [Seq.nEchos,1,Seq.SeqAverage.average]), 3) ...
        - Seq.SeqAverage.TXPhaseInvertIncrement;
      tempAQPhaseIncrement = ...
        cumsum(repmat(Seq.SeqAverage.AQPhaseIncrement, [1+Seq.nEchos,1,Seq.SeqAverage.average]), 3) ...
        - Seq.SeqAverage.AQPhaseIncrement(1);
    else
      tempTXPhaseInvertIncrement = ...
        cumsum(repmat(Seq.SeqAverage.TXPhaseInvertIncrement, [1,Seq.nEchos,Seq.SeqAverage.average]), 3) ...
        - Seq.SeqAverage.TXPhaseInvertIncrement;
      tempAQPhaseIncrement = ...
        cumsum(repmat(Seq.SeqAverage.AQPhaseIncrement, [1,1+Seq.nEchos,Seq.SeqAverage.average]), 3) ...
        - Seq.SeqAverage.AQPhaseIncrement(1);
    end
    tempTXPhaseFlipIncrement = ...
      cumsum(repmat(Seq.SeqAverage.TXPhaseFlipIncrement, [1,1,Seq.SeqAverage.average]), 3) ...
      - Seq.SeqAverage.TXPhaseFlipIncrement;

    if Seq.fast
      tempTXPhaseInvertOffset = repmat(permute(Seq.SeqAverage.TXPhaseInvertOffset,[1,3,2]), [Seq.nEchos,1,1]);
      tempAQPhaseOffset = repmat(permute(Seq.SeqAverage.TXPhaseFlipOffset,[1,3,2]), [1+Seq.nEchos,1,1]);
    else
      tempTXPhaseInvertOffset = repmat(permute(Seq.SeqAverage.TXPhaseInvertOffset,[1,3,2]), [1,Seq.nEchos,1]);
      tempAQPhaseOffset = repmat(permute(Seq.SeqAverage.TXPhaseFlipOffset,[1,3,2]), [1,1+Seq.nEchos,1]);
    end
    tempTXPhaseFlipOffset = repmat(permute(Seq.SeqAverage.TXPhaseFlipOffset,[1,3,2]), [1,1,1]);

    Seq.tRep = repmat([Seq.tRep(1:end-1), Seq.tRep(end)+Seq.SeqAverage.averageBreak], 1, Seq.SeqAverage.average);
    Seq.tRep(end) = Seq.tRep(end)-Seq.SeqAverage.averageBreak;
    for iTX = idxTX_common
      TX(iTX).Start = repmat(TX(iTX).Start, 1, Seq.SeqAverage.average);
      TX(iTX).Duration = repmat(TX(iTX).Duration, 1, Seq.SeqAverage.average);
      TX(iTX).Frequency = repmat(TX(iTX).Frequency, 1, Seq.SeqAverage.average);
      TX(iTX).Amplitude = repmat(TX(iTX).Amplitude, 1, Seq.SeqAverage.average);
    end

    % Digital IO (might be set by prepare function or damp coil)
    if isfield(Seq, 'DigitalIO')
      for iDevice = 1:numel(Seq.DigitalIO)
        if ~isemptyfield(Seq.DigitalIO(iDevice), 'SetTime')
          Seq.DigitalIO(iDevice).SetTime = ...
            repmat(Seq.DigitalIO(iDevice).SetTime, 1, Seq.SeqAverage.average);
        end
        if ~isemptyfield(Seq.DigitalIO(iDevice), 'SetValue')
          Seq.DigitalIO(iDevice).SetValue = ...
            repmat(Seq.DigitalIO(iDevice).SetValue, 1, Seq.SeqAverage.average);
        end
        if ~isemptyfield(Seq.DigitalIO(iDevice), 'Repeat')
          Seq.DigitalIO(iDevice).Repeat = ...
            repmat(Seq.DigitalIO(iDevice).Repeat, 1, Seq.SeqAverage.average);
        end
      end
    end

    % TX.Amplitude(1:end/2)=0;
    if isfield(TX(1), 'Repeat')
      for iTX = idxTX_common
        TX(iTX).Repeat = repmat(TX(iTX).Repeat, 1, Seq.SeqAverage.average);
      end
    end
    AQ.Start = repmat(AQ.Start, 1, Seq.SeqAverage.average);
    AQ.fSample = repmat(AQ.fSample, 1, Seq.SeqAverage.average);
    AQ.SamplingFactor = repmat(AQ.SamplingFactor, 1, Seq.SeqAverage.average);

    if isempty(Seq.tsT2T2)
       for iTX = idxTX_common
        if Seq.fast
          TX(iTX).Phase = repmat(TX(iTX).Phase, 1, Seq.SeqAverage.average) + ...
            reshape(cat(1, ...
                        tempTXPhaseFlipIncrement + tempTXPhaseFlipOffset, ...
                        tempTXPhaseInvertIncrement + tempTXPhaseInvertOffset), [], Seq.SeqAverage.average);
        else
          TX(iTX).Phase = repmat(TX(iTX).Phase, 1, Seq.SeqAverage.average) + ...
            repmat(reshape(cat(2, ...
                               tempTXPhaseFlipIncrement + tempTXPhaseFlipOffset, ...
                               tempTXPhaseInvertIncrement + tempTXPhaseInvertOffset), 1, []), ...
                           size(TX(iTX).Phase,1),1);
        end
      end

      if Seq.fast
        AQ.Phase = repmat(AQ.Phase, 1, Seq.SeqAverage.average) + ...
          reshape(tempAQPhaseIncrement+tempAQPhaseOffset, [], Seq.SeqAverage.average);
      else
        AQ.Phase = repmat(AQ.Phase, 1, Seq.SeqAverage.average) + ...
          repmat(reshape(tempAQPhaseIncrement+tempAQPhaseOffset, 1, []), size(AQ.Start,1), 1);
      end

      if Seq.dualNuclear
        % FIXME: Do we need to support a different phase cycle for the X
        % frequency?
        if Seq.fast
          AQ.PhaseX = repmat(AQ.PhaseX, 1, Seq.SeqAverage.average) + ...
            reshape(tempAQPhaseIncrement+tempAQPhaseOffset, [], Seq.SeqAverage.average);
        else
          AQ.PhaseX = repmat(AQ.PhaseX, 1, Seq.SeqAverage.average) + ...
            repmat(reshape(tempAQPhaseIncrement+tempAQPhaseOffset, 1, []), size(AQ.Start,1), 1);
        end
      end

    else
      tempTXPhaseInvertIncrementT2T2 = ...
        cumsum(repmat(Seq.SeqAverage.TXPhaseInvertIncrement, ...
                      [1,Seq.nEchosT2T2,Seq.SeqAverage.average]), 3) - ...
        Seq.SeqAverage.TXPhaseInvertIncrement;
      tempAQPhaseIncrementT2T2 = ...
        cumsum(repmat(Seq.SeqAverage.AQPhaseIncrement, [1,1+Seq.nEchosT2T2,Seq.SeqAverage.average]), 3) - ...
        Seq.SeqAverage.AQPhaseIncrement(1);

      tempTXPhaseInvertOffsetT2T2=repmat(permute(Seq.SeqAverage.TXPhaseInvertOffset,[1,3,2]),[1,Seq.nEchosT2T2,1]);
      tempAQPhaseOffsetT2T2=repmat(permute(Seq.SeqAverage.TXPhaseFlipOffset,[1,3,2]),[1,1+Seq.nEchosT2T2,1]);

      for iTX = idxTX_common
        TX(iTX).Phase = repmat(TX(iTX).Phase, 1, Seq.SeqAverage.average) + ...
          repmat(reshape(cat(2, tempTXPhaseFlipIncrement+tempTXPhaseFlipOffset, ...
                                tempTXPhaseInvertIncrement+tempTXPhaseInvertOffset, ...
                                tempTXPhaseFlipIncrement+tempTXPhaseFlipOffset, ...
                                tempTXPhaseFlipIncrement+tempTXPhaseFlipOffset, ...
                                tempTXPhaseFlipIncrement+tempTXPhaseFlipOffset, ...
                                tempTXPhaseInvertIncrementT2T2+tempTXPhaseInvertOffsetT2T2...
                            ), 1, []), size(TX(iTX).Phase,1), 1);
      end
      AQ.Phase = repmat(AQ.Phase, 1, Seq.SeqAverage.average) + ...
        repmat(reshape(cat(2, tempAQPhaseIncrement+tempAQPhaseOffset, ...
                              tempAQPhaseIncrement(:,1,:)+tempAQPhaseOffset(:,1,:), ...
                              tempAQPhaseIncrement(:,1,:)+tempAQPhaseOffset(:,1,:), ...
                              tempAQPhaseIncrementT2T2+tempAQPhaseOffsetT2T2...
                          ), 1, []), size(AQ.Start,1),1);
      if Seq.dualNuclear
        AQ.PhaseX = repmat(AQ.PhaseX, 1, Seq.SeqAverage.average) + ...
          repmat(reshape(cat(2, tempAQPhaseIncrement+tempAQPhaseOffset, ...
                                tempAQPhaseIncrement(:,1,:)+tempAQPhaseOffset(:,1,:), ...
                                tempAQPhaseIncrement(:,1,:)+tempAQPhaseOffset(:,1,:), ...
                                tempAQPhaseIncrementT2T2+tempAQPhaseOffsetT2T2...
                            ), 1, []), size(AQ.Start,1),1);
      end
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

    % Grad (might have been added by prepare function)
    for iGrad = 1:numel(Grad)
      if ~isemptyfield(Grad(iGrad), 'Time') ...
          && size(Grad(iGrad).Time, 2) == numel(Seq.tRep)/Seq.SeqAverage.average
        Grad(iGrad).Time = repmat(Grad(iGrad).Time, 1, Seq.SeqAverage.average);
        Grad(iGrad).Amp = repmat(Grad(iGrad).Amp, 1, Seq.SeqAverage.average);
      end
    end
  end


  if Seq.dualNuclear
    AQ.FrequencyX = AQ.Frequency;
    AQ.FrequencyX(~isnan(AQ.FrequencyX)) = Seq.AQFrequencyX;

    % Select power sync frequency for which the DC-DC peaks are possibly
    % furthest away for both AQ frequencies.
    [~, iMin] = min(abs(mod(Seq.AQFrequency ./ HW.MMRT(SliceSelect.iDevice).PowerSyncFrequencyList, 1) - 0.5) ...
                    + abs(mod(Seq.AQFrequencyX ./ HW.MMRT(SliceSelect.iDevice).PowerSyncFrequencyList, 1) - 0.5));
    AQ.PowerSyncFrequency = HW.MMRT.PowerSyncFrequencyList(iMin);
    % AQ.AutoPowerSyncFrequency = 0;
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

%% Fixme
%  if Seq.UnblankTx2
%    UnblankTx2Offset=1000e-3;
%    TX(2).Start=TX(1).Start(1,:)-[UnblankTx2Offset,zeros(1,size(TX(1).Start,2)-1)];
%    TX(2).Duration=(Seq.tRep-50e-6)+[UnblankTx2Offset,zeros(1,size(TX(1).Start,2)-1)];
%    TX(2).Amplitude=0;
%    TX(2).BlankPostset=500e-6;
%  end

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
      iAQ = find([SeqOut.AQ(:).Device] == SliceSelect.iDevice, 1, 'first');
      iData = find([data(:).device] == SliceSelect.iDevice);

      for iAQX = 1:numel(iData)
        if isempty(data(iAQX).data) || (isscalar(data(iAQX).data) && isnan(data(iAQX).data))
          continue;
        end
        if Seq.SeqAverage.SaveSeqAverageData
          data(iAQX).SeqAverage.data = reshape(data(iAQX).data, ...
            [size(data(iAQX).data,1), size(SeqOut.AQ(iAQ).Start,1), ...
             size(data(iAQX).data,3)/Seq.SeqAverage.average, Seq.SeqAverage.average]);
        end
        data(iAQX).data = mean(reshape(data(iAQX).data, ...
          [size(data(iAQX).data,1), size(SeqOut.AQ(iAQ).Start,1), ...
           size(data(iAQX).data,3)/Seq.SeqAverage.average, Seq.SeqAverage.average]), 4);
        data(iAQX).time_all = reshape(data(iAQX).time_all, ...
          [size(data(iAQX).data,1), size(SeqOut.AQ(iAQ).Start,1), ...
           size(data(iAQX).data,3), Seq.SeqAverage.average]);
        data(iAQX).time_all = data(iAQX).time_all(:,:,:,1);
        data(iAQX).time_of_tRep = reshape(data(iAQX).time_of_tRep, ...
          [size(data(iAQX).data,1), size(SeqOut.AQ(iAQ).Start,1), ...
           size(data(iAQX).data,3), Seq.SeqAverage.average]);
        data(iAQX).time_of_tRep = data(iAQX).time_of_tRep(:,:,:,1);
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
    if Seq.dualNuclear
      error('PD:sequence_EchoStandard:RawDataDual', ...
        'RawData mode currently not supported in dual nuclear measurements.');
    end
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
