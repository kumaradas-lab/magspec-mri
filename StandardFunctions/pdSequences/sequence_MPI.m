function [SeqOut, dataOut, data, data_1D] = sequence_MPI(HW, Seq)
%% MPI experiment with travelling FFL (field free line)
%
%   [SeqOut, dataOut, data, data_1D] = sequence_MPI(HW, Seq)
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
%     tRep
%           Repetition time in case of multiple (segmented?) acquisitions in s.
%           (Default: HW.MPI.tRep or approximate minimum repetition time)
%
%     tAQ
%           Duration of the acquisition in seconds.
%           (Default: HW.MPI.tAQ or 20e-3)
%
%     kf1
%           Number of periods in the acquisition time for the rotational
%           frequency. If it is a scalar, the set value applies to all gradients
%           set in Seq.GradMPS. It must have the same number of elements as
%           Seq.GradMPS otherwise. Set to 0 for no modulation in this direction
%           (i.e., no rotation of the FFL).
%           (Default: HW.MPI.kf1 or 1)
%
%     kf2
%           Number of periods in the acquisition time for the travel
%           (longitudinal) frequency. If it is a scalar, the set value applies
%           to all gradients set in Seq.GradMPS. It must have the same number of
%           elements as Seq.GradMPS otherwise.
%           (Default: HW.MPI.kf2 or 11)
%
%     kf3
%           Number of periods in the acquisition time for the swipe frequency.
%           It is emitted on an rf channel (#5).
%           (Default: HW.MPI.kf3 or 111)
%
%     GradMPS
%           Vector with the indices for the gradients to be used where x = 1,
%           y = 2, z = 3, B = 4. (Default: 4)
%
%     ampGrad
%           Amplitude of the gradient fields in Tesla/meter. If it is a scalar,
%           the set value applies to all gradients set in Seq.GradMPS. It must
%           have the same number of elements as Seq.GradMPS otherwise.
%           Seq.ampGrad, Seq.IinGrad, and Seq.UinGrad are mutually exclusive.
%           (Default in case Seq.IinGrad and Seq.UinGrad are not set: 0.01)
%
%     IinGrad
%           Amplitude of the current through the gradient coils in Ampere. If it
%           is a scalar, the set value applies to all gradients set in
%           Seq.GradMPS. It must have the same number of elements as Seq.GradMPS
%           otherwise. Seq.ampGrad, Seq.IinGrad, and Seq.UinGrad are mutually
%           exclusive. This property can be used in case the corresponding
%           gradient channel is current controlled (e.g. the DC-600 gradient
%           amplifier). (Default: not set)
%
%     UinGrad
%           Amplitude of the voltage at the gradient coils in Volts. If it is a
%           scalar, the set value applies to all gradients set in Seq.GradMPS.
%           It must have the same number of elements as Seq.GradMPS otherwise.
%           Seq.ampGrad, Seq.IinGrad, and Seq.UinGrad are mutually exclusive.
%           This property can be used in case the corresponding gradient channel
%           is voltage controlled (e.g. the "Grad" of the drive-l).
%           (Default: not set)
%
%     tRampUpGrad
%           Duration of the (sin^2) ramp up of the gradient amplitude before the
%           acquisition in seconds. If it is a scalar, the set value applies to
%           all gradients set in Seq.GradMPS. It must have the same number of
%           elements as Seq.GradMPS otherwise.
%           (Default: HW.MPI.tRampUpGrad or tAQ/kf1)
%
%     tRampDownGrad
%           Duration of the (sin^2) ramp down of the gradient amplitude after
%           the acquisition in seconds. If it is a scalar, the set value applies
%           to all gradients set in Seq.GradMPS. It must have the same number of
%           elements as Seq.GradMPS otherwise.
%           (Default: HW.MPI.tRampDownGrad or tAQ/kf1)
%
%     tRampUpHoldGrad
%           Hold time of the gradient amplitude after the ramp up and before the
%           acquisition in seconds. If it is a scalar, the set value applies to
%           all gradients set in Seq.GradMPS. It must have the same number of
%           elements as Seq.GradMPS otherwise.
%           (Default: HW.MPI.tRampUpHoldGrad or
%           (AQmode>1)*tAQ/kf1 - tRampUpGrad)
%
%     tRampDownHoldGrad
%           Hold time of the gradient amplitude after the acquisition and before
%           the ramp down in seconds. If it is a scalar, the set value applies
%           to all gradients set in Seq.GradMPS. It must have the same number of
%           elements as Seq.GradMPS otherwise.
%           (Default: HW.MPI.tRampDownHoldGrad or
%           (AQmode>0)*tAQ/kf1 - tRampDownGrad)
%
%     tSampleGrad
%           Sample rate for generating the gradient pulses in seconds. If it is
%           a scalar, the set value applies to all gradients set in Seq.GradMPS.
%           It must have the same number of elements as Seq.GradMPS otherwise.
%           (Default: HW.MPI.tSampleGrad or as fast as possible on a
%           6e-6 seconds grid while not exceeding 1000 gradient points per tRep)
%
%     phaseGrad
%           Phase of the sinusoidal gradient pulse in radians. Applies to the
%           oscillation corresponding to kf2. If it is a scalar, the set value
%           applies to all gradients set in Seq.GradMPS. It must have the same
%           number of elements as Seq.GradMPS otherwise.
%           (Default: HW.MPI.phaseGrad or 0)
%
%     phaseGradEnvelope
%           Phase in radians of the sinusoidal envelope for the oscillating
%           gradient signal (corresponding to kf1). Only applies if Seq.kf1 is
%           not empty. If it is a scalar, the set value applies to all gradients
%           set in Seq.GradMPS. It must have the same number of elements as
%           Seq.GradMPS otherwise.
%           (Default: HW.MPI.phaseGradEnvelope or [])
%
%     nPhaseGrad
%           Number of phase increment steps for the gradient channels. This
%           value applies to all gradients set in Seq.GradMPS. (Default: 1)
%
%     phaseGradIncrement
%           Increment for phaseGrad (in radians). This value applies to all
%           gradients set in Seq.GradMPS. (Default: 2*pi/Seq.nPhaseGrad)
%
%     ampSwipe
%           Deflection of FFL in swipe direction in Tesla/meter.
%           (Default: 1)
%
%     ampTX
%           rf amplitude (channel #5, kf3) in Tesla.
%           (Default: Seq.ampSwipe*HW.MPI.swipeRadius)
%
%     tRampUpTX
%           Duration of the (linear) ramp up of the rf signal (#5) before the
%           acquisition in seconds. (Default: HW.MPI.tRampUpTX or 0)
%
%     tRampDownTX
%           Duration of the (linear) ramp down of the rf signal (#5) after the
%           acquisition in seconds. (Default: HW.MPI.tRampDownTX or 0)
%
%     tRampUpHoldTX
%           Hold time of the rf amplitude after the ramp up and before the
%           acquisition in seconds.
%           (Default: HW.MPI.tRampUpHoldTX or (AQmode>1)*tAQ/kf1 - tRampUpTX)
%
%     tRampDownHoldTX
%           Hold time of the rf amplitude after the acquisition and before the
%           ramp down in seconds.
%           (Default: HW.MPI.tRampDownHoldTX or (AQmode>0)*tAQ/kf1 - tRampDownTX)
%
%     phaseTX
%           Phase of rf pulse in radians. (Default: 0)
%
%     nPhaseTX
%           Number of phase increment steps for the rf pulse. The phase is
%           incremented each time the phase incrementation of the phase for the
%           gradient channels has completed a full cycle. (Default: 1)
%
%     phaseTXIncrement
%           Increment for phaseTX (in radians). (Default: 2*pi ./ Seq.nPhaseTX)
%
%     channelTX
%           Used TX channel. (Default: 2)
%
%     invertOnOtherChannelTX
%           If all TX properties are scalars and this parameter is true, send
%           the same signal that is emitted on Seq.channelTX inverted (i.e.,
%           with 180 degrees phase difference) on the other channel.
%           (Default: false)
%
%     invertOnOtherChannelTX_phase
%           Additional phase in degrees for the inverted signal that is emitted
%           on the other channel. This only applies if
%           Seq.invertOnOtherChannelTX is true and can apply.
%           (Default: 0)
%
%     invertOnOtherChannelTX_factor
%           Factor for the amplitude of the inverted signal that is emitted on
%           the other channel. This only applies if Seq.invertOnOtherChannelTX
%           is true and can apply. (Default: 0)
%
%     phaseCH14
%           Additional (potentially hardware-specific) phase of the sinusoidal
%           gradient pulse in radians. This phase is additive to Seq.phaseGrad
%           (see above).
%           (Default: HW.MPI.phaseCH14 or 0)
%
%     phaseCH5
%           Additional (potentially hardware-specific) phase of the rf pulse
%           in radians. This phase is additive to Seq.phaseTX (see above).
%           (Default: HW.MPI.phaseCH5 or 0)
%
%     nPrepare
%           Number of tReps without acquisition (but with the MPI gradient
%           pulses as set) before the tReps with acquisition. (Default: 1)
%
%     nMeasurements
%           Number of tReps with acquisition. (Default: prod(Seq.nPhaseGrad) )
%
%     Gain
%           Acquisition gain (must be between HW.RX.GainMin and HW.RX.GainMax).
%           (Default:
%           get_RX_Amplitude(HW, 1/HW.MPI.RX_Uin_max, 'Uin', 'AQ', AQ_Gain)
%           @ 0 Hz center frequency or HW.RX.GainMax/10)
%
%     fSampleAQ
%           Sampling frequency of the acquisition in samples per second.
%           (Default: HW.MPI.fSampleAQ or 2.5e6)
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
%           periods of the slowest MPS frequency are not acquired. For a value
%           <0, the acquisition is omitted.
%
%     AQDataMode
%           Select data mode (bit width) for transmission of received signal
%           between console and PC.
%           (Default: 3*(HW.MMRT.FPGA_Firmware >= 20221129) )
%
%     GetDataEarly
%           Add packet end command to each tRep with acquisition to receive data
%           as soon as it has been acquired.
%           (Default: false)
%
%     Device
%           Index of the connected MRT device that is used for the experiment.
%           (Default: 1)
%
%     trigger
%           Structure with settings for an optional trigger signal containing
%           the following fields:
%
%       use
%             Boolean value that indicates whether a trigger signal should be
%             added to the pulse program.
%             (Default: HW.MPI.trigger.use)
%
%       duration
%             Duration of the trigger signal in seconds.
%             (Default: HW.MPI.trigger.duration)
%
%       offset
%             Time of the rising edge of the trigger signal relative to the
%             start of the CH5 rampup in seconds.
%             (Default: HW.MPI.trigger.offset)
%
%       output
%             Value of the digital output channels that is set while the trigger
%             signal is high.
%             (Default: HW.MPI.trigger.output)
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
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------


%% FIXME:
% * Add option for "repeated" line travel in acquisition window?
% * Add options to "segment" FFL trajectory over multiple acquisition windows?
% * Add options for different ramp types for resonant coils (e.g.,
%   pre-emphasis).
% * Restore AQmode=0? (Probably, not working currently.)


%% input check
if nargin ~= 2
  error('PD:sequence_MPI:nargin', 'Number of input arguments must be 2.');
end

if isemptyfield(Seq, 'PreProcessSequence'),  Seq.PreProcessSequence = true;  end
if isemptyfield(Seq, 'StartSequence'),  Seq.StartSequence = true;  end
if isemptyfield(Seq, 'PollPPGfast'),  Seq.PollPPGfast = true;  end
if isemptyfield(Seq, 'GetRawData'),  Seq.GetRawData = true;  end
if isemptyfield(Seq, 'PostProcessSequence'),  Seq.PostProcessSequence = true;  end

if ~isemptyfield(Seq, 'AQ'),  AQ = Seq.AQ;  end
if ~isemptyfield(Seq, 'TX'),  TX = Seq.TX;  end
if ~isemptyfield(Seq, 'Grad'),  Grad = Seq.Grad;  end

if Seq.PreProcessSequence
  %% default parameters
  if isemptyfield(Seq, 'Device'), Seq.Device = 1; end  % index for used MRT
  if numel(HW.MMRT) < Seq.Device
    error('PD:sequence_MPI:InvalidDevice', ...
      'Seq.Device (%d) higher than number of connected MRT devices (%d).\n', ...
      Seq.Device, numel(HW.MMRT));
  end
  if isemptyfield(Seq, 'GradMPS'), Seq.GradMPS = 4; end  % index of gradient connected to MPS coil

  % gradient field amplitude parameters
  if isemptyfield(Seq, 'ampGrad') && isemptyfield(Seq, 'IinGrad') && isemptyfield(Seq, 'UinGrad')
    % gradient field amplitude in T/m
    Seq.ampGrad = 0.01 * ones(size(Seq.GradMPS));
  end
  if ~isemptyfield(Seq, 'ampGrad') && ~isemptyfield(Seq, 'IinGrad') && ~isemptyfield(Seq, 'UinGrad')
    error('PD:sequence_MPI:AmpIinUinSet', ...
      'Seq.AmpGrad, Seq.IinGrad, and Seq.UinGrad are mutually exclusive.');
  end
  if ~isemptyfield(Seq, 'IinGrad') ...
      && any(HW.Grad(Seq.Device).PaCurrentControlled(HW.Grad(Seq.Device).xyzB(Seq.GradMPS)) ~= 1)
    error('PD:sequence_MPI:IinNotCurrentControlled', ...
      'Seq.IinGrad can only be used if the corresponding channels are current controlled.');
  end
  if ~isemptyfield(Seq, 'UinGrad') ...
      && any(HW.Grad(Seq.Device).PaCurrentControlled(HW.Grad(Seq.Device).xyzB(Seq.GradMPS)) ~= 0)
    error('PD:sequence_MPI:UinNotVoltageControlled', ...
      'Seq.UinGrad can only be used if the corresponding channels are voltage controlled.');
  end

  if isemptyfield(Seq, 'phaseGrad')
    % phase of gradient oscillation in rad
    if isempty(HW.MPI.phaseGrad)
      Seq.phaseGrad = 0;
    else
      Seq.phaseGrad = HW.MPI.phaseGrad;
    end
  end
  if isemptyfield(Seq, 'phaseCH14')
    % additional phase of gradient oscillation in rad
    if isempty(HW.MPI.phaseCH14)
      Seq.phaseCH14 = 0;
    else
      Seq.phaseCH14 = HW.MPI.phaseCH14;
    end
  end
  if isemptyfield(Seq, 'nPhaseGrad'), Seq.nPhaseGrad = 1; end  % number of phase steps for gradient channels
  if isemptyfield(Seq, 'phaseGradIncrement'), Seq.phaseGradIncrement = 2*pi ./ Seq.nPhaseGrad; end  % phase increment of gradient oscillation in rad
  if isemptyfield(Seq, 'RampGrad'), Seq.RampGrad = 1; end  % ramp gradient amplitude at start and end of tRep for given number of periods

  if isemptyfield(Seq, 'phaseGradEnvelope')
    % phase of gradient envelope in rad
    if isempty(HW.MPI.phaseGrad)
      Seq.phaseGradEnvelope = 0;
    else
      Seq.phaseGradEnvelope = HW.MPI.phaseGradEnvelope;
    end
  end

  % acquisition parameters
  if isemptyfield(Seq, 'tAQ')
    % duration of acquisition
    if isempty(HW.MPI.tAQ)
      Seq.tAQ = 20e-3;
    else
      Seq.tAQ = HW.MPI.tAQ;
    end
  end
  if isemptyfield(Seq, 'AQmode')
    % 0 for smallest gap, 1 for one slow period at end,
    % 2: customizable
    Seq.AQmode = 2;
  end
  omitAQ = false;
  if Seq.AQmode < 0
    Seq.AQmode = abs(Seq.AQmode);
    if Seq.AQmode < 1
       Seq.AQmode = 0;
    end
    omitAQ = true;
  end
  if isemptyfield(Seq, 'AQDataMode')
    % data mode (bit width) for transmission of received signal between console and PC
    Seq.AQDataMode = 3*(HW.MMRT(Seq.Device).FPGA_Firmware >= 20221129);
  end
  if isemptyfield(Seq, 'GetDataEarly')
    % get data as early as possible after each acquisition
    Seq.GetDataEarly = false;
  end

  if isemptyfield(Seq, 'fAQ'), Seq.fAQ = 0; end  % AQ mixing frequency in Hz
  if isemptyfield(Seq, 'Gain')
    % AQ Gain (must be between HW.RX.GainMin and HW.RX.GainMax)
    if isempty(HW.MPI.RX_Uin_max)
      Seq.Gain = HW.RX(Seq.Device).GainMax / 10;
    else
      AQ_Gain.Frequency = Seq.fAQ;  % temporary structure
      Seq.Gain = get_RX_Amplitude(HW, 1/HW.MPI.RX_Uin_max, 'Uin', 'AQ', AQ_Gain);
    end
  end
  if isemptyfield(Seq, 'fSampleAQ')
    % AQ sample rate in Hz
    if isempty(HW.MPI.fSampleAQ)
      Seq.fSampleAQ = 2.5e6;
    else
      Seq.fSampleAQ = HW.MPI.fSampleAQ;
    end
  end
  % round to integer reduction factor
  Seq.fSampleAQ = HW.RX(Seq.Device).fSample ...
    ./ round(HW.RX(Seq.Device).fSample ./ Seq.fSampleAQ);
  % round to integer number of samples
  Seq.tAQ = round(Seq.tAQ .* Seq.fSampleAQ) ./ Seq.fSampleAQ;

  % gradient field modulation parameters
  if isemptyfield(Seq, 'kf1')
    % number of periods during AQ for rotational frequency (f1, "envelope" Grad)
    if isempty(HW.MPI.kf1)
      Seq.kf1 = 1;
    else
      Seq.kf1 = HW.MPI.kf1;
    end
  end
  if isemptyfield(Seq, 'kf2')
    % number of periods during AQ for travel frequency (f2, "fast" Grad)
    if isempty(HW.MPI.kf2)
      Seq.kf2 = 11;
    else
      Seq.kf2 = HW.MPI.kf2;
    end
  end
  if isemptyfield(Seq, 'kf3')
    % number of periods during AQ for swipe frequency (f3, TX)
    if isempty(HW.MPI.kf3)
      Seq.kf3 = 111;
    else
      Seq.kf3 = HW.MPI.kf3;
    end
  end

  % expand gradient properties to all used gradient channels if needed
  nGrads = numel(Seq.GradMPS);
  if numel(Seq.kf1) < nGrads, Seq.kf1 = Seq.kf1(1) * ones(size(Seq.GradMPS)); end
  if numel(Seq.kf2) < nGrads, Seq.kf2 = Seq.kf2(1) * ones(size(Seq.GradMPS)); end

  if Seq.kf1 == 0
    Seq.fGradEnvelope = [];  % no envelope
  else
    Seq.fGradEnvelope = Seq.kf1/Seq.tAQ;  % frequency of "envelope" in Hertz
  end
  Seq.fGrad = Seq.kf2/Seq.tAQ;  % travel frequency in Hertz
  Seq.fTX = Seq.kf3/Seq.tAQ;  % rf frequency at TX in Hertz

  if isemptyfield(Seq, 'tRampUpGrad')
    % ramp up duration for gradient amplitude in seconds
    if isempty(HW.MPI.tRampUpGrad)
      Seq.tRampUpGrad = 1./Seq.fGradEnvelope;
    else
      Seq.tRampUpGrad = HW.MPI.tRampUpGrad;
    end
  end
  if isemptyfield(Seq, 'tRampDownGrad')
    % ramp up duration for gradient amplitude in seconds
    if isempty(HW.MPI.tRampDownGrad)
      Seq.tRampDownGrad = 1./Seq.fGradEnvelope;
    else
      Seq.tRampDownGrad = HW.MPI.tRampDownGrad;
    end
  end
  if isemptyfield(Seq, 'tRampUpHoldGrad')
    % hold duration for gradient amplitude after ramp up in seconds
    if isempty(HW.MPI.tRampUpHoldGrad)
      Seq.tRampUpHoldGrad = (Seq.AQmode>1)/Seq.fGradEnvelope - Seq.tRampUpGrad;
    else
      Seq.tRampUpHoldGrad = HW.MPI.tRampUpHoldGrad;
    end
  end
  if isemptyfield(Seq, 'tRampDownHoldGrad')
    % hold duration for gradient amplitude before ramp down in seconds
    if isempty(HW.MPI.tRampDownHoldGrad)
      Seq.tRampDownHoldGrad = (Seq.AQmode>0)/Seq.fGradEnvelope - Seq.tRampDownGrad;
    else
      Seq.tRampDownHoldGrad = HW.MPI.tRampDownHoldGrad;
    end
  end

  if numel(Seq.tRampUpGrad) < nGrads
    Seq.tRampUpGrad = Seq.tRampUpGrad(1) * ones(size(Seq.GradMPS));
  end
  if numel(Seq.tRampDownGrad) < nGrads
    Seq.tRampDownGrad = Seq.tRampDownGrad(1) * ones(size(Seq.GradMPS));
  end
  if numel(Seq.tRampUpHoldGrad) < nGrads
    Seq.tRampUpHoldGrad = Seq.tRampUpHoldGrad(1) * ones(size(Seq.GradMPS));
  end
  if numel(Seq.tRampDownHoldGrad) < nGrads
    Seq.tRampDownHoldGrad = Seq.tRampDownHoldGrad(1) * ones(size(Seq.GradMPS));
  end

  if isemptyfield(Seq, 'tSampleGrad')
    % sample rate for generating gradient pulses in seconds
    if isempty(HW.MPI.tSampleGrad)
      Seq.tSampleGrad = ceil(...
        max(Seq.tRampUpGrad + Seq.tRampUpHoldGrad + Seq.tAQ ...
            + Seq.tRampDownHoldGrad + Seq.tRampDownGrad) / 1000 / 6e-6) * 6e-6;
    else
      Seq.tSampleGrad = HW.MPI.tSampleGrad;
    end
  end

  Seq.startGrad = -Seq.tRampUpGrad - Seq.tRampUpHoldGrad;

  Seq.tEndGrad = Seq.tAQ + (Seq.AQmode>0).*max(Seq.tRampDownHoldGrad+Seq.tRampDownGrad);  % end time of oscillations on gradient channels in seconds

  % rf field parameters
  Seq.startTX = [];  % start of rf pulse(s) in seconds
  if ~isfield(Seq, 'ampSwipe'), Seq.ampSwipe = []; end  % deflection in swipe direction in Tesla/meter
  if ~isfield(Seq, 'ampTX'), Seq.ampTX = []; end  % rf amplitude at TX in Tesla
  if ~isfield(Seq, 'phaseTX'), Seq.phaseTX = []; end  % phase of rf pulse in radians
  if ~isfield(Seq, 'phaseCH5'), Seq.phaseCH5 = []; end   % additional phase correction in rad for channel 5 (rf channel)
  if ~isfield(Seq, 'nPhaseTX'), Seq.nPhaseTX = []; end  % number of phase steps for rf transmission
  if ~isfield(Seq, 'phaseTXIncrement'), Seq.phaseTXIncrement = []; end  % number of phase steps for rf transmission
  if ~isfield(Seq, 'tRampUpTX'), Seq.tRampUpTX = []; end  % ramp up duration for rf transmission in seconds
  if ~isfield(Seq, 'tRampDownTX'), Seq.tRampDownTX = []; end  % ramp down duration for rf transmission in seconds
  if ~isfield(Seq, 'tRampUpHoldTX'), Seq.tRampUpHoldTX = []; end  % hold duration for rf transmission after ramp up in seconds
  if ~isfield(Seq, 'tRampDownHoldTX'), Seq.tRampDownHoldTX = []; end  % hold duration for rf transmission before ramp downin seconds
  if ~isfield(Seq, 'channelTX'), Seq.channelTX = []; end  % channel for rf transmission
  if isemptyfield(Seq, 'invertOnOtherChannelTX')
    % send inverted signal additionally on the other TX channel
    Seq.invertOnOtherChannelTX = false;
  end
  if isemptyfield(Seq, 'invertOnOtherChannelTX_phase')
    % additional phase for the signal on the inverted channel in degrees
    Seq.invertOnOtherChannelTX_phase = 0;
  end
  if isemptyfield(Seq, 'invertOnOtherChannelTX_factor')
    % factor for amplitude on inverted channel
    Seq.invertOnOtherChannelTX_factor = 0;
  end

  nTX = max([numel(Seq.fTX), numel(Seq.ampSwipe), numel(Seq.ampTX), numel(Seq.phaseTX), ...
    numel(Seq.tRampUpTX), numel(Seq.tRampDownTX), ...
    numel(Seq.tRampUpHoldTX), numel(Seq.tRampDownHoldTX)]);
  if nTX > 0
    if isempty(Seq.fTX), Seq.fTX = 10e-3; end
    if isempty(Seq.channelTX), Seq.channelTX = 2; end
    if isempty(Seq.ampSwipe)
      Seq.ampSwipe = 1;
    end
    if isempty(Seq.ampTX)
      Seq.ampTX = Seq.ampSwipe * HW.MPI.swipeRadius;
    end
    if isempty(Seq.phaseTX), Seq.phaseTX = 0; end
    if isemptyfield(Seq, 'phaseCH5')
      % additional phase of rf pulse in rad
      if isempty(HW.MPI.phaseCH5)
        Seq.phaseCH5 = 0;
      else
        Seq.phaseCH5 = HW.MPI.phaseCH5;
      end
    end
    if isempty(Seq.nPhaseTX), Seq.nPhaseTX = 1; end
    if isempty(Seq.phaseTXIncrement), Seq.phaseTXIncrement = 2*pi ./ Seq.nPhaseTX; end
    if isempty(Seq.tRampUpTX)
      if isempty(HW.MPI.tRampUpTX)
        Seq.tRampUpTX = 0;
      else
        Seq.tRampUpTX = HW.MPI.tRampUpTX;
      end
    end
    if isempty(Seq.tRampDownTX)
      if isempty(HW.MPI.tRampDownTX)
        Seq.tRampDownTX = 0;
      else
        Seq.tRampDownTX = HW.MPI.tRampDownTX;
      end
    end
    if isempty(Seq.tRampUpHoldTX)
      if isempty(HW.MPI.tRampUpHoldTX)
        Seq.tRampUpHoldTX = (Seq.AQmode>1)/Seq.fGradEnvelope - Seq.tRampUpTX;
      else
        Seq.tRampUpHoldTX = HW.MPI.tRampUpHoldTX;
      end
    end
    if isempty(Seq.tRampDownHoldTX)
      if isempty(HW.MPI.tRampDownHoldTX)
        Seq.tRampDownHoldTX = (Seq.AQmode>0)/Seq.fGradEnvelope - Seq.tRampDownTX;
      else
        Seq.tRampDownHoldTX = HW.MPI.tRampDownHoldTX;
      end
    end
    if isempty(Seq.startTX)
      if Seq.AQmode >= 2
        Seq.startTX = -Seq.tRampUpTX - Seq.tRampUpHoldTX;
      else
        Seq.startTX = -(Seq.AQmode-1)/Seq.fGradEnvelope;
      end
    end

    if numel(Seq.startTX) < nTX, Seq.startTX = ones(1, nTX) * Seq.startTX(1); end
    if numel(Seq.fTX) < nTX, Seq.fTX = ones(1, nTX) * Seq.fTX(1); end
    if numel(Seq.ampSwipe) < nTX, Seq.ampSwipe = ones(1, nTX) * Seq.ampSwipe(1); end
    if numel(Seq.ampTX) < nTX, Seq.ampTX = ones(1, nTX) * Seq.ampTX(1); end
    if numel(Seq.phaseTX) < nTX, Seq.phaseTX = ones(1, nTX) * Seq.phaseTX(1); end
    if numel(Seq.phaseCH5) < nTX, Seq.phaseCH5 = ones(1, nTX) * Seq.phaseCH5(1); end
    if numel(Seq.nPhaseTX) < nTX, Seq.nPhaseTX = ones(1, nTX) * Seq.nPhaseTX(1); end
    if numel(Seq.phaseTXIncrement) < nTX, Seq.phaseTXIncrement = ones(1, nTX) * Seq.phaseTXIncrement(1); end
    if numel(Seq.tRampUpTX) < nTX, Seq.tRampUpTX = ones(1, nTX) * Seq.tRampUpTX(1); end
    if numel(Seq.tRampDownTX) < nTX, Seq.tRampDownTX = ones(1, nTX) * Seq.tRampDownTX(1); end
    if numel(Seq.tRampUpHoldTX) < nTX, Seq.tRampUpHoldTX = ones(1, nTX) * Seq.tRampUpHoldTX(1); end
    if numel(Seq.tRampDownHoldTX) < nTX, Seq.tRampDownHoldTX = ones(1, nTX) * Seq.tRampDownHoldTX(1); end
    if numel(Seq.channelTX) < nTX
      error('PD:sequence_MPI:UnknownChannelTX', ...
        'Seq.channelTX must be defined if multiple rf channels are used.');
    end
  end

  % timing parameters
  if isemptyfield(Seq, 'nPrepare'), Seq.nPrepare = 1; end  % number of tReps without AQ before tReps with AQ
  if isemptyfield(Seq, 'nMeasurements'), Seq.nMeasurements = Seq.nPhaseGrad * Seq.nPhaseTX; end  % number of tReps with AQ

  if any(Seq.phaseGradIncrement ~= 0)
    Seq.CLTime = 100e-6;
  else
    Seq.CLTime = 0.8e-6;
  end
  if isemptyfield(Seq, 'tRep')
    if isempty(HW.MPI.tRep)
      % estimate minimum duration for one tRep
      Seq.tRep = ...
        max(Seq.tEndGrad + max(HW.Grad(Seq.Device).TimeDelay(Seq.GradMPS))*(Seq.AQmode~=0), ...  % maximum end of ramp down gradient channels
            Seq.tAQ + max(Seq.tRampDownHoldTX + Seq.tRampDownTX + HW.TX(Seq.Device).Latenz + HW.TX(Seq.Device).BlankPostset)*(Seq.AQmode~=0)) ...  % maximum end of ramp down TX channels
        + Seq.CLTime*(Seq.AQmode~=0) ...
        + ((Seq.nPrepare+Seq.nMeasurements)>0) ...  % additionally if there are multiple tReps
          * (max(max(Seq.tRampUpTX + Seq.tRampUpHoldTX), ...  % maximum ramp up TX channels
                 max(Seq.tRampUpGrad + Seq.tRampUpHoldGrad)) ...  % maximum ramp up gradient channels
             + 1e-3);  % tOffset "buffer" (see below)
    else
      Seq.tRep = HW.MPI.tRep;
    end
  end

  % trigger
  if isemptyfield(Seq, {'trigger', 'use'})
    % Boolean value to indicate whether trigger signal should be used
    Seq.trigger.use = HW.MPI.trigger.use;
  end
  if isemptyfield(Seq, {'trigger', 'duration'})
    % duration of the trigger signal in seconds
    Seq.trigger.duration = HW.MPI.trigger.duration;
  end
  if isemptyfield(Seq, {'trigger', 'offset'})
    % time of the rising edge of the trigger relative to the start of the CH5
    % rampup in seconds
    Seq.trigger.offset = HW.MPI.trigger.offset;
  end
  if isemptyfield(Seq, {'trigger', 'output'})
    % digital output value for trigger
    Seq.trigger.output = HW.MPI.trigger.output;
  end

  if isemptyfield(Seq, 'average'), Seq.average = 1; end  % number of averages
  if isemptyfield(Seq, 'averageBreak'), Seq.averageBreak = []; end  % break between averages in seconds (must be > 1 ms)
  if isemptyfield(Seq, 'plotData'), Seq.plotData = 1; end  % plot acquired data

  % expand gradient properties to all used gradient channels if needed
  nGrads = numel(Seq.GradMPS);
  if ~isemptyfield(Seq, 'ampGrad') && (numel(Seq.ampGrad) < nGrads)
    Seq.ampGrad = Seq.ampGrad(1) * ones(size(Seq.GradMPS));
  end
  if ~isemptyfield(Seq, 'UinGrad') && (numel(Seq.UinGrad) < nGrads)
    Seq.UinGrad = Seq.UinGrad(1) * ones(size(Seq.GradMPS));
  end
  if ~isemptyfield(Seq, 'IinGrad') && (numel(Seq.IinGrad) < nGrads)
    Seq.IinGrad = Seq.IinGrad(1) * ones(size(Seq.GradMPS));
  end
  if numel(Seq.tSampleGrad) < nGrads, Seq.tSampleGrad = Seq.tSampleGrad(1) * ones(size(Seq.GradMPS)); end
  if numel(Seq.phaseGrad) < nGrads, Seq.phaseGrad = Seq.phaseGrad(1) * ones(size(Seq.GradMPS)); end
  if numel(Seq.phaseCH14) < nGrads, Seq.phaseCH14 = Seq.phaseCH14(1) * ones(size(Seq.GradMPS)); end
  if ~isscalar(Seq.phaseGradIncrement)
    error('PD:sequence_MPI:NonScalarPhaseGradIncrement', ...
        'Seq.phaseGradIncrement must have a scalar value.');
  end
  if ~isempty(Seq.fGradEnvelope)
    if numel(Seq.fGradEnvelope) < nGrads
      Seq.fGradEnvelope = Seq.fGradEnvelope(1) * ones(size(Seq.GradMPS));
    end
    if numel(Seq.phaseGradEnvelope) < nGrads
      Seq.phaseGradEnvelope = Seq.phaseGradEnvelope(1) * ones(size(Seq.GradMPS));
    end
  end


  %% pulse program
  Seq.tRep = Seq.tRep * ones(1, Seq.nPrepare+Seq.nMeasurements);  % tRep

  AQ.fSample = Seq.fSampleAQ;
  AQ.Start = [NaN(1, Seq.nPrepare), ...
    zeros(1, Seq.nMeasurements)];
  AQ.nSamples = round(Seq.tAQ * AQ.fSample);
  AQ.DataMode = Seq.AQDataMode;
  AQ.Device = Seq.Device;
  AQ.GetData = false(size(Seq.tRep));
  if Seq.GetDataEarly
    AQ.GetData(Seq.nPrepare+1:end) = true;
  end

  % Seq.MissingSamples = Seq.tEndGrad*AQ.fSample(1)-AQ.nSamples(1);
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

  Seq.tEndAQ = AQ.nSamples / AQ.fSample;

  oldPaEnable = HW.Grad(Seq.Device).PaEnable;
  guard = onCleanup(@() ResetPaEnable(HW, oldPaEnable, Seq.Device));
  HW.Grad(Seq.Device).PaEnable = 1;

  Seq.tOffset = max(max(-Seq.startTX), max(-Seq.startGrad)) + 1e-3;

  % rf transmission
  % FIXME: Correctly handle signal on multiple TX channels.
  [phaseIncGradAll, phaseIncTXAll] = ...
    ndgrid(Seq.phaseGradIncrement *(0:1:(Seq.nPhaseGrad-1)), ...
           Seq.phaseTXIncrement *(0:1:(Seq.nPhaseTX-1)));
  [TX(1:max(1,nTX)).Start] = deal(NaN);
  [TX(1:max(1,nTX)).Device] = Seq.Device;

  % round to integer reduction factor
  % FIXME: That is the actually used rf frequency anyway.  As a consequence of
  %        this rounding, the period of the rf pulse might not be an integer
  %        divider of the AQ duration.  That might lead to "frequency leaking" in
  %        the spectrum.  Is that an issue?
  FrequencyGridTX = HW.TX(Seq.Device).fSample / 2^HW.TX(Seq.Device).DdsPicBits;
  Seq.fTX = round(Seq.fTX / FrequencyGridTX) * FrequencyGridTX;

  for iTX = 1:nTX
    TX(iTX).Channel = Seq.channelTX(iTX);
    % reference phase at start of acquisition window
    phaseTX = Seq.phaseTX(iTX) + Seq.phaseCH5(iTX) - Seq.tOffset*Seq.fTX(iTX)*2*pi;

    % linear ramp up or ramp down
    numSegments = 51;  % FIXME: Make configurable?
    tRampUp = linspace(0, Seq.tRampUpTX(iTX), numSegments).';
    tRampDown = linspace(0, Seq.tRampDownTX(iTX), numSegments).';
    TX(iTX).Start = ...
      [Seq.startTX(iTX)+tRampUp; ...
       Seq.tEndAQ+Seq.tRampDownHoldTX(iTX)+tRampDown(1:end-1)-100e-9];
    TX(iTX).Amplitude = ...
      [tRampUp(1:end-1)/Seq.tRampUpTX(iTX)*Seq.ampTX(iTX); ...
       Seq.ampTX(iTX)
       flipud(tRampDown(1:end-1))/Seq.tRampDownTX(iTX)*Seq.ampTX(iTX)];
    TX(iTX).Frequency = repmat(Seq.fTX(iTX), size(TX(iTX).Amplitude,1), 1);
    TX(iTX).Phase = repmat(rad2deg(phaseTX), size(TX(iTX).Amplitude,1), 1);
    TX(iTX).Duration = ...
      [repmat(Seq.tRampUpTX(iTX)/(numSegments-1), numSegments-1, 1); ...
       Seq.tRampUpHoldTX(iTX)+Seq.tEndAQ+Seq.tRampDownHoldTX(iTX)-100e-9; ...
       repmat(Seq.tRampDownTX(iTX)/(numSegments-1), numSegments-1, 1);];

    % remove duplicate time points
    iRemove = [diff(TX(iTX).Start) < 4/HW.TX(Seq.Device).fSample; false];
    % keep first of duplicates
    iRemove([diff(iRemove)>0; false]) = false;
    % remove segments with NaN amplitude
    iRemove(isnan(TX(iTX).Amplitude)) = true;
    % iRemove(1,:) = false;
    TX(iTX).Start(iRemove) = [];
    TX(iTX).Amplitude(iRemove) = [];
    TX(iTX).Frequency(iRemove) = [];
    TX(iTX).Phase(iRemove) = [];
    TX(iTX).Duration(iRemove) = [];

    % repeat (with phase increment) for all tReps
    TX(iTX).Start = repmat(TX(iTX).Start, 1, numel(Seq.tRep));
    TX(iTX).Amplitude = repmat(TX(iTX).Amplitude, 1, numel(Seq.tRep));
    TX(iTX).Frequency = repmat(TX(iTX).Frequency, 1, numel(Seq.tRep));
    TX(iTX).Phase = repmat(TX(iTX).Phase, 1, numel(Seq.tRep));
    phaseInc = [(-Seq.nPrepare:1:-1)*Seq.phaseTXIncrement, ...
      reshape(phaseIncTXAll, 1, [])]/pi*180;
    TX(iTX).Phase = bsxfun(@plus, TX(iTX).Phase, phaseInc);
    TX(iTX).Duration = repmat(TX(iTX).Duration, 1, numel(Seq.tRep));
  end

  if nTX == 1 && Seq.invertOnOtherChannelTX
    TX(2) = TX(1);
    TX(2).Channel = mod(TX(2).Channel, 2) + 1;
    TX(2).Phase = TX(2).Phase + 180 + Seq.invertOnOtherChannelTX_phase;
    TX(2).Amplitude = TX(2).Amplitude * Seq.invertOnOtherChannelTX_factor;
  else
    % TX(2).Amplitude = NaN;
    % TX(2).Start = NaN;
  end

  offsetGrad_device = sum([HW.Grad(1:(Seq.Device-1)).n]);
  % differential output voltage at the Grad output (static offset)
  %               |
  %               v
  Grad(offsetGrad_device+1).Shim = 0.0 ./ HW.Grad(Seq.Device).Amp2LoadUin(HW.Grad(Seq.Device).x); % the HW.Grad.x indexing ist used, because HW.Grad.Amp2PaUout is DAC channel indexed and Grad(1).Shim is [x y z B0] = [1 2 3 4] indexed.
  Grad(offsetGrad_device+2).Shim = 0.0 ./ HW.Grad(Seq.Device).Amp2LoadUin(HW.Grad(Seq.Device).y); % so if you don't use HW.Grad.x indexing, you get in trouble if HW.Grad.x ~= 1 or HW.Grad.y ~= 2 or HW.Grad.z ~= 3 or HW.Grad.B ~= 4
  Grad(offsetGrad_device+3).Shim = 0.0 ./ HW.Grad(Seq.Device).Amp2LoadUin(HW.Grad(Seq.Device).z);
  Grad(offsetGrad_device+4).Shim = 0.0 ./ HW.Grad(Seq.Device).Amp2LoadUin(HW.Grad(Seq.Device).B);

  % differential output voltage at the Grad output (dynamic)
  % Generate sine wave on selected channels
  for iGradMPS = 1:numel(Seq.GradMPS(:))
    Gn = Seq.GradMPS(iGradMPS);
    tEndGrad = Seq.tRampUpGrad(iGradMPS) + Seq.tRampUpHoldGrad(iGradMPS) ...
      + Seq.tEndAQ + Seq.tRampDownHoldGrad(iGradMPS) + Seq.tRampDownGrad(iGradMPS);
    tGrid = (0:Seq.tSampleGrad(iGradMPS):tEndGrad) + Seq.startGrad(iGradMPS);
    if (tGrid(end) - Seq.tEndGrad) == 0, tGrid(end) = []; end
    % Grad(Gn).Time = [linspace(0, Seq.tEndGrad, round(Seq.tEndGrad/Seq.tSampleGrad+1)).' + HW.Grad.TimeDelay ];
    phaseInc = [(-Seq.nPrepare:1:-1)*Seq.phaseGradIncrement, ...
      reshape(phaseIncGradAll, 1, [])];
    Grad(Gn).Time = (tGrid.') * ones(1, size(Seq.tRep,2));
    Grad(Gn).Amp = sin(bsxfun(@plus, 2*pi*(Grad(Gn).Time(:,1))*Seq.fGrad(iGradMPS), ...
      phaseInc+Seq.phaseGrad(iGradMPS)+Seq.phaseCH14(iGradMPS)-AQ.Start(1,Seq.nPrepare+1)*Seq.fGrad(iGradMPS)*2*pi));
    if ~isempty(Seq.fGradEnvelope)
      Grad(Gn).Amp = bsxfun(@times, Grad(Gn).Amp, ...
        sin(Grad(Gn).Time(:,1) * 2*pi * Seq.fGradEnvelope(iGradMPS) ...
            + Seq.phaseGradEnvelope(iGradMPS) ...
            - AQ.Start(1,Seq.nPrepare+1)*Seq.fGradEnvelope(iGradMPS)*2*pi));
    end
    if isfield(Seq, 'UinGrad')
      Grad(Gn).Amp = Seq.UinGrad(iGradMPS) * Grad(Gn).Amp * HW.Grad.LoadUin2Amp(iGradMPS);
    elseif isfield(Seq, 'IinGrad')
      Grad(Gn).Amp = Seq.IinGrad(iGradMPS) * Grad(Gn).Amp * HW.Grad.LoadIin2Amp(iGradMPS);
    else
      Grad(Gn).Amp = Seq.ampGrad(iGradMPS) * Grad(Gn).Amp;
    end
    % Grad(Gn).Repeat = [0, ones(1, Seq.nMeasurements+Seq.nPrepare-1), 0];

    if Seq.tRampUpGrad(iGradMPS) > 0
      % ramp up amplitude at start of tRep
      isFirstCycle = (Grad(Gn).Time(:,1) - Grad(Gn).Time(1,1)) <= Seq.tRampUpGrad(iGradMPS);
      Grad(Gn).Amp(isFirstCycle,:) = Grad(Gn).Amp(isFirstCycle,:) .* ...
        sin(2*pi * 1/Seq.tRampUpGrad(iGradMPS)/4 * bsxfun(@minus, Grad(Gn).Time(isFirstCycle,:), Grad(Gn).Time(1,:))).^2;
    end
    if Seq.tRampDownGrad(iGradMPS) > 0
      % ramp down amplitude at end of tRep
      isLastCycle = (Grad(Gn).Time(end,1) - Grad(Gn).Time(:,1)) <= Seq.tRampDownGrad(iGradMPS);
      Grad(Gn).Amp(isLastCycle,:) = Grad(Gn).Amp(isLastCycle,:) .* ...
        cos(2*pi * 1/Seq.tRampDownGrad(iGradMPS)/4 * bsxfun(@minus, Grad(Gn).Time(isLastCycle,1),Grad(Gn).Time(find(isLastCycle,1, 'first'),:))).^2;
    elseif Seq.tEndGrad + 2*Seq.tSampleGrad(iGradMPS) < Seq.tRep(1)
      Grad(Gn).Amp(end,:) = 0;
    end
    Grad(Gn).Amp(end,end) = 0;
  end

  if omitAQ
    AQ.Start(:) = NaN;
  end

  if Seq.trigger.use
    Seq.DigitalIO.SetTime = NaN(2, numel(Seq.tRep));
    % timing is relative to start of CH5 signal
    Seq.DigitalIO.SetTime(1,:) = TX(1).Start(1,:) + Seq.trigger.offset;
    Seq.DigitalIO.SetTime(2,:) = Seq.DigitalIO.SetTime(1,:) + Seq.trigger.duration;
    Seq.DigitalIO.SetValue = Seq.DigitalIO.SetTime;
    Seq.DigitalIO.SetValue(1,:) = Seq.trigger.output;
    Seq.DigitalIO.SetValue(2,:) = 0;
  end

  % some important characteristics for the gradient programming:
  % - the first Amp (Grad.Amp(1,:)) is copied to the beginning of each tRep,
  % - the last used Amp per tRep is copied to the end of each tRep.
  % - It's not (yet) allowed to have Grad ramp to the next tRep, only a constant value is allowed.
end

if Seq.PreProcessSequence || Seq.StartSequence || Seq.PollPPGfast || Seq.GetRawData
  % HW.Grad.AmpUnit={'V' 'V' 'V' 'V'};  % Change the Text at the plots to a proper unit
  % HW.Grad.AmpUnitScale=1./HW.Grad.Amp2PaUout(HW.Grad.xyzB); % Change the scale to match the units


  % while 1 % start sequence and loop until pressing Strg+C
  if nargout > 3 || Seq.plotData
    [~, SeqOut, data, data_1D] = set_sequence(HW, Seq, AQ, TX, Grad);
  else
    [~, SeqOut, data] = set_sequence(HW, Seq, AQ, TX, Grad);
  end
  %   sleep(0.2); % restart sequence after less than 2 Seconds to avoid turning off the gradients by the auto Mute
  % end
  % MRIUnit.MRISequency.exctractArrayToFile(talker.mySequency.getCommandArray,'test.txt')
end

if ~Seq.PostProcessSequence || (isscalar(data.data) && isnan(data.data))
  % return early
  dataOut = [];
  return;
end

if Seq.plotData
  HW.RX(SeqOut.Device).AmplitudeUnitScale = 1/1000/data(SeqOut.Device).Amplitude2Uin(SeqOut.nPrepare+1);
  HW.RX(SeqOut.Device).AmplitudeUnit = 'mV';
  plot_data_1D(HW, data_1D(SeqOut.Device));

  clear dataOut
  iAQ = find([SeqOut.AQ.Device] == SeqOut.Device, 1, 'first');
  dataOut.fAQ = SeqOut.AQ(iAQ).Frequency(1,Seq.nPrepare+1);
  dataOut.fSAQ = SeqOut.AQ(iAQ).fSample(1,Seq.nPrepare+1);

  switch Seq.AQmode
    case 0
      dataOut.time_all1D = linspace(data(1).time_all(1,1,Seq.nPrepare+1), ...
        data(1).time_all(end,1,Seq.nPrepare+1), ...
        round((data(1).time_all(end,1,Seq.nPrepare+1) - data(1).time_all(1,1,Seq.nPrepare+1)) * SeqOut.AQ(iAQ).fSample(1,Seq.nPrepare+1))+1).';
      dataOut.data1D = interp1(data(1).time_all(~isnan(data(1).time_all)), ...
        data(1).data(~isnan(data(1).time_all)), dataOut.time_all1D, 'spline') ...
        * data(1).Amplitude2Norm(Seq.nPrepare+1);
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
      dataOut.nAQ = length(data(1).time_all(~isnan(data(1).time_all)));
      dataOut.tAQ = (1/dataOut.fSAQ)*dataOut.nAQ;
      dataOut.time_all1D = linspace(1/dataOut.fSAQ,dataOut.tAQ,dataOut.nAQ).';
      dataOut.df = 1/dataOut.tAQ;
      dataOut.fstart = -dataOut.df*floor(dataOut.nAQ/2) + dataOut.fAQ;
      dataOut.fstop = dataOut.df*floor(dataOut.nAQ/2-0.5) + dataOut.fAQ;
      dataOut.data1D = data(1).data(~isnan(data(1).time_all)) * data(1).Amplitude2Norm(Seq.nPrepare+1);
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

      % wrap into "equivalent" repetitions (i.e., multiple passes of the same trajectory)
      nSamples = size(data(1).data, 1);
      maxSamplesPerPeriod = nSamples ...
        / gcd(nSamples, gcd(gcd(min(Seq.kf1), min(Seq.kf2)), min(Seq.kf3)));
      dataOut.time = data(1).time_of_tRep((1:maxSamplesPerPeriod)',SeqOut.nPrepare+1) ...
        - data(1).time_of_tRep(1,SeqOut.nPrepare+1);
      dataOut.data = reshape(data(1).data(:,:,SeqOut.nPrepare+1:end), ...
        maxSamplesPerPeriod, ...
        size(data(1).time_of_tRep,1)/size(dataOut.time,1), [])*data(1).Amplitude2Uin(SeqOut.nPrepare+1);
      dataOut.meanData = squeeze(mean(dataOut.data,2));
      dataOut.fft1_data = ifftshift(ifft(dataOut.meanData,[],1));
      dataOut.f_fft1_data = get_FFTGrid(SeqOut.AQ(iAQ).fSample(1), maxSamplesPerPeriod);


      % hf = figure(5); clf(hf);
      % hax(1) = subplot(3,1,1, 'Parent', hf);
      % plot(hax(1), dataOut.time, reshape(dataOut.data,maxSamplesPerPeriod,[])*1000);
      % hold(hax(1), 'on');
      % plot(hax(1), ...
      %   dataOut.time([1;end]), ...
      %   ones(2,1)/data(1).Amplitude2Norm(SeqOut.nPrepare+1) ...
      %   * data(1).Amplitude2Uin(SeqOut.nPrepare+1) * 1000 ...
      %   * (1+double(SeqOut.AQ(iAQ).Frequency(SeqOut.nPrepare+1)==0)), ...  % SeqOut.AQ.Frequency==0 -> double amplitude in time signal
      %   '--k');
      % ylabel(hax(1), 'Amp in mV');
      % hold(hax(1), 'off');
      % grid(hax(1), 'on');
      %
      % hax(2) = subplot(3,1,2, 'Parent', hf);
      % semilogy(hax(2), dataOut.f_fft1_data, abs(dataOut.fft1_data)*1000, 'LineWidth', 2);
      % % set(hax(2), 'XTick', 0:1:dataOut.fstop/SeqOut.fGrad(1));
      % % set(hax(2), 'XTickMode', 'manual');
      % ylabel(hax(2), 'Amp in mV');
      % grid(hax(2), 'on');
      %
      % hax(3) = subplot(3,1,3, 'Parent', hf);
      % plot(hax(3),  dataOut.f_fft1_data, angle(dataOut.fft1_data), 'LineWidth', 2);
      % % set(hax(3), 'XTick', 0:1:dataOut.fstop/SeqOut.fGrad(1));
      % % set(hax(3), 'XTickMode', 'manual');
      % set(hax(3), 'XGrid', 'on');
      % title(hax(3), ['Seq.fGrad = ' num2str(SeqOut.fGrad) ' Hz']);
      % xlabel(hax(3), '"Harmonics" of Seq.fGrad(1)');
      % ylabel(hax(3), 'Phase in RAD');
      %
      % linkaxes(hax(2:3), 'x');
      % xlim(hax(3), [-10, 10]*max(SeqOut.fGrad));

  end
else
  dataOut = [];
end

end


function ResetPaEnable(HW, newValue, iDevice)
%% function to (re)set HW.Grad.PaEnable
% used in onCleanup callback
HW.Grad(iDevice).PaEnable = newValue;
end
