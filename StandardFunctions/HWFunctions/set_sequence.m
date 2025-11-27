function [raw_data_AQs, Seq, data_S_AQs_TRs, data_1D, data_Cell, averages] = set_sequence(HW, Seq, AQ, TX, Grad)
%% Send pulse program to MMRT and receive data
%
%   [raw_data_AQs, Seq, data_S_AQs_TRs, data_1D, data_Cell, averages] = ...
%       set_sequence(HW, Seq, AQ, TX, Grad)
%
%
% This function is the main interface between Matlab and the MMRT for sending
% and receiving data. The pulse program is defined by the structures AQ, TX and
% Grad (see the manual for further details). Additionally, fields in the
% structure Seq are used to determine the pulse program (see below).
%
%
%
% INPUT:
%
%   HW
%         Structure or PD.HWClass object as created when executing LoadSystem.
%
%
%   Seq
%         Structure with pulse program settings. Amongst others, the following
%         fields can be used (default values apply if the fields are omitted or
%         empty):
%
%     tRep
%           Row vector with the duration of repetition blocks in seconds. Pulse
%           programs can be divided into repetition blocks. The commands for
%           each repetition block are loaded at the beginning of the repetition
%           block (see CLTime below). The remainder of the tRep is available for
%           pulse program activity.
%
%     tOffset
%           Row vector with the offsets of the tReps in seconds. These offsets
%           must be non-negative. This allows having pulse program elements at
%           negative times for any actions in the tReps. A tOffset effectively
%           reduces the active time in the previous tRep (i.e. there cannot be
%           any actions close to the end of the previous tRep). If the values in
%           tOffset are too small to account for all pulse program elements, it
%           is automatically extended to the minimum necessary value.
%           (Default: 0)
%
%     CLTime
%           Row vector with the command load times of each tRep in seconds.
%           Before each tRep, the commands have to be loaded by the MMRT. During
%           this time, the state of the MMRT cannot change. This effectively
%           limits the active time of the previous tRep. By default, the minimum
%           time necessary is calculated and used.
%
%     Reinitialize
%           Boolean value. If true, the MMRT is re-initialized before the
%           sequence is started. I.e. among others, the clocks necessary for
%           timing are re-started and the gradient amplifiers are switched on in
%           case they have been muted before. (default: HW.ReInit)
%
%     tRepInit
%           Double scalar with the time in seconds for the tRep of the
%           initialization package. (Default: HW.tRepInit)
%
%     StartSequenceTime
%           Double scalar with the time when the sequence should be started in
%           Matlab "now" format. The actual start of the sequence might be
%           slightly later depending on the system load. If omitted or empty,
%           the sequence is started as soon as possible. In the output, this
%           field contains the (approximate) Matlab time immediately before
%           starting the sequence. (Default input value: [])
%
%     TimeFromLastSequence
%           Double scalar with the time in seconds measured from the end of the
%           last tRep in the previous sequence to t=0 in the first tRep of the
%           current sequence. Timing between sequences with this value is exact
%           with the precission of the quartz in the MMRT (or an external
%           source). This is only effective if Seq.Reinitialize is false.
%           Additionally, you might need to pass the value of Seq.LoopCountEnd
%           of the previous sequence as Seq.LoopCountStart of the current
%           sequence. By default, the sequence starts as soon as possible
%           (Seq.TimeFromLastSequence = []).
%
%     TimeToNextSequence
%           Double scalar with the time in seconds measured from the end of the
%           last tRep in the current sequence to t=0 in the first tRep of the
%           next sequence.
%
%     average
%           Scalar with the number of averages using the same pulse program. The
%           measurement is repeated Seq.average times and the averaged values of
%           all measurements are returned. The layout of the returned data is as
%           if there was only one measurement. (Default: 1)
%
%     averageRepetitionTime
%           Double scalar with the time in seconds measured from the start of
%           the first tRep to the start of the first tRep of the next averaging
%           step. See Seq.averageBreak. (Default: [])
%
%     averageBreak
%           Double scalar with the time in seconds measured from the end of the
%           last tRep to the start of the next averaging step taking the tOffset
%           of the first tRep into account. If Seq.averageBreak and
%           Seq.averageRepetitionTime are both set, they must match. If only one
%           of the two is set, the other one will be set accordingly. (Default
%           if neither are set: 0)
%
%     fLarmor
%           If no frequency is defined in the AQ or TX structures (see below),
%           default to this value for the mixer or rf frequency, respectively.
%           (Default: HW.fLarmor)
%
%     DigitalIO
%           Structure with settings for the digital output on the Digital I/O
%           port of the drive-l. The following fields are used:
%       SetTime
%             2d double matrix with the time in seconds for the digital output
%             action. Each column corresponds to one tRep.
%       SetValue
%             2d double matrix with the 6bit value for the digital output. The
%             size must match Seq.DigitalIO.SetTime.
%
%     tRepIsTriggered
%             Boolean vector the same size as Seq.tRep which indicates which
%             tReps are waiting for an external trigger (high level at digital
%             IN1). If this is a vector smaller than Seq.tRep, zeros are
%             appended. (Default: false(size(Seq.tRep)))
%
%     tRepTriggerTimeOut
%             Row vector the same size as Seq.tRep with the timeout in seconds
%             for the corresponding trigger. Only those values where
%             Seq.tRepIsTriggered is true have any effect. If this is a scalar,
%             the timeout applies to all triggered tReps. (Default: 5)
%
%     tRepTriggerDebounceTime
%             Row vector the same size as Seq.tRep with the debounce time for
%             the trigger event in seconds. That means the trigger must be in
%             high state at least for this amount of time before the tRep starts
%             executing. If this is a scalar, the debounce time applies to all
%             triggered tReps. (Default: HW.MMRT.TriggerDebounceTime)
%
%     SyncMaster
%             Boolean value to indicate whether the (primary) device will wait
%             for a trigger signal at digital IN1. The trigger sent at digital
%             OUT6 at the start of the pulse program can be used for that.
%             (Default: 0)
%
%     SyncDevices
%             Boolean value to indicate whether secondary devices will be
%             synchronized. If true, all secondary devices wait for a trigger
%             at digital IN1 to start their pulse programs at the same time. The
%             trigger sent (by the primary device) at digital OUT6 at the start
%             of its pulse program can be used for that. (Default: 0)
%
%     PPGStartTimeOut
%             Maximum time in seconds between sending the start signal for the
%             pulse program to the device and the actual start of the pulse
%             program. This time-out might need to be set to a higher value if
%             it is expected that the time until receiving a trigger signal
%             might be longer than the default. See, e.g., Seq.SyncMaster.
%             (Default: 5)
%
%     CommandWindowOutput
%           Boolean value. If false, (most) output to the command window is
%           suppressed. (Default: true)
%
%     plotSeq
%           Row vector with the numbers of the gradient channels to be included
%           in the pulse program plot. "0" meaning no gradient channels to be
%           included in that plot. "[]" suppresses the pulse program plot
%           completely. See documentation of function "plotSeq" for further
%           details. (Default: [])
%
%     plotSeqTR
%           Same as Seq.plotSeq. But the pulse program plot is wrapped at each
%           tRep. See documentation of function "plotSeqTR" for further details.
%           (Default: [])
%
%     PostProcessingFcn
%           Function handle for overloading the post-processing function. By
%           default, that is @Postprocessing which calls get_data and
%           get_data_1D. The function prototype for this handle must be
%           compatible to:
%           [Seq, data_1D, data_S_AQs_TRs] = ...
%             @(nargout_set_sequence, Seq, raw_data_AQs, num_averages, ...
%               Channel, iDevice)
%           (Default: HW.PostProcessingFcn)
%
%     AuxCommands
%           Optional matrix with the same number of colums as Seq.tRep that
%           contains FPGA commands (64-bit unsigned integers) that are prepended
%           to the commands for each respective tRep.
%           (Default: [])
%
%     AuxCommandsInit
%           Optional column vector that contains FPGA commands (64-bit unsigned
%           integers) that are prepended to the commands for the initialization
%           package.
%           (Default: [])
%
%
%   The AQ, TX, and Grad structures contain fields with matrices that act as
%   "timetables" for defining the pulse program. Each column of these matrices
%   corresponds to one tRep in the pulse program. If settings are provided for a
%   single tRep only, they are repeated in all tReps. See the manual for further
%   details.
%
%   AQ
%       Structure (or array of structures) with settings for the acquisition
%       windows in the pulse program. Among others, it might have the following
%       optional fields:
%
%     Start
%         Array with start times of the acquisition window in seconds (relative
%         to the t=0 time of the corresponding tRep).
%
%     nSamples
%         Array with number of samples in the acquisition windows.
%
%     fSample
%         Array with sample frequencies of the acquisition windows in Hertz.
%         This frequency might be rounded to the closest possible sample
%         frequency realizable by the mixer. (I.e., the system frequency
%         HW.RX.fSample divided by an integer in the range from
%         HW.RX.CIC_Decimation_Min to HW.RX.CIC_Decimation_Max.) Alternatively,
%         it can be set to the system frequency to get the raw data (without
%         mixing).
%
%     Frequency
%         Array with mixer frequencies for the acquisition windows in Hertz.
%
%     Phase
%         Array with phases for mixing the signal in degrees.
%
%     SamplingFactor
%         Array with (software) sampling factors that can be used to reach
%         sample frequencies that are lower than the one that can be realized
%         with the firmware mixer (see AQ.fSample). If this value is 1, no
%         additional (software) downsampling is used. If set to 0, the lowest
%         possible sampling factor is chosen to come close to the set sample
%         frequency (AQ.fSample). If set to Inf, a higher sampling factor might
%         be chosen if it allows to reach the set sample frequency *exactly*. If
%         that isn't possible, it acts the same as 0. (Default: 0)
%
%     ResetPhases
%         Row vector with Boolean values to indicate tReps for which the phase
%         reference for TX pulses and RX mixer is re-started.
%
%     Gain
%         Scalar with receiver gain (1/AQ.Gain = max input Amplitude).
%
%     Repeat
%         Row vector with Boolean values to indicate that all actions from the
%         previous tRep are repeated identically. That might reduce the time
%         needed for loading the commands between tReps (CLTime).
%
%     Device
%         Scalar with the index identifying the device if multiple devices are
%         connected simultaneously.
%
%     FrequencyX
%         Array with the secondary mixer frequency for each acquisition window
%         in Hertz. This can be used to mix the received signal to two different
%         center frequencies, e.g., when running dual frequency experiments at
%         1H and at an X nucleus resonance.
%         NOTE: This feature might not be supported by all configurations.
%         (Default: [], i.e., only the mixer for the primary frequency is used)
%
%     PhaseX
%         Array with phases for mixing the signal at the secondary frequency in
%         degrees.
%         NOTE: This feature might not be supported by all configurations.
%         (Default: [], i.e., only the mixer for the primary frequency is used)
%
%     DataMode
%         NOTE: This is an experimental feature that is not supported by all
%               configurations. Use with care!
%         Select a mode for the bit-width of the data transferred via USB from
%         the console to the PC. If AQ.FrequencyX is set, DataMode 1 is used
%         unconditionally. For raw data mode (AQ.fSample == HW.RX.fSample),
%         DataMode 0 is used unconditionally.
%         0: 1 AQ channel (mixer) with 48 bits for the real and imaginary parts,
%            respectively.
%         1: 2 AQ channels (mixers) with each 24 bits for the real and imaginary
%            parts, respectively.
%         2: 1 AQ channel (mixer) with 24 bits for the real and imaginary parts,
%            respectively. This allows to approximately double the maximum
%            sampling rate compared to data mode 0.
%         3: 1 AQ channel (mixer) with 24 bits for the real part (discarding the
%            imaginary part of the mixer output). This allows to approximately
%            quadruple the maximum sampling rate compared to data mode 0.
%         (Default: HW.RX(AQ.Device).DataMode)
%
%     GetDigitalInTime
%         NOTE: This is an experimental feature that is not supported by all
%               configurations. Use with care!
%         Get the last three time stamps for the rising and falling edge at the
%         digital input channels 1, 2, and 3 of the console. (See "get_data" for
%         more details.)
%         (Default: 0)
%
%
%   TX
%       Structure (or array of structures) with settings for rf pulses in the
%       pulse program. Among others, it might have the following optional
%       fields:
%
%     Channel
%         Scalar with an index for the output channel (1: TRx, 2: Tx2).
%
%     Start
%         Array with start times of the rf pulses in seconds (relative to the
%         t=0 time of the corresponding tRep).
%
%     Frequency
%         Array with frequencies of the rf pulse in Hertz.
%
%     Duration
%         Array with durations of the rf pulses in seconds.
%
%     Amplitude
%         Array with rf pulse amplitudes. If calibrated correctly, these
%         amplitudes are in Tesla by default.
%
%     Phase
%         Array with phases of the rf pulses in degrees.
%
%     BlankOffset
%         Row vector with offset times for unblanking the transmitter before an
%         rf pulse in seconds.
%
%     BlankPostset
%         Row vector with postset times for unblanking the transmitter after an
%         rf pulse in seconds.
%
%     AmplitudeDC
%         Vector with amplitudes of DC output during rf pulse per tRep. The same
%         calibration values as for TX.Amplitude are used. If this has any
%         non-zero elements, restrictions apply on switching between different
%         output ports for signals generated by one DDs.
%         (Default: no DC output during rf pulses)
%
%     SkipAmplitudeCheck
%         Optional array with the same size as TX.Amplitude with Boolean values
%         that indicate whether any checks on the amplitude should be skipped
%         for a given pulse (segment).
%         (Default: [])
%
%     Repeat
%         Row vector with Boolean values to indicate that all actions from the
%         previous tRep are repeated identically. That might reduce the time
%         needed for loading the commands between tReps (CLTime).
%
%     Device
%         Scalar with the index identifying the device if multiple devices are
%         connected simultaneously.
%
%
%   Grad
%       Array of structures with settings for gradient pulses in the pulse
%       program. By default, the order of this array corresponds to the x, y,
%       and z gradient followed by the fourth gradient channel. Among others,
%       they might have the following optional fields:
%
%     Time
%         Array with times of gradient points in seconds (relative to the
%         t=0 time of the corresponding tRep).
%
%     Amp
%         Array with amplitudes of the gradient points in T/m. The gradient
%         amplitude is linearly interpolated between these points.
%
%     Shim
%         Row vector with shim values for each tRep in T/m.
%
%     Repeat
%         Row vector with Boolean values to indicate that all actions from the
%         previous tRep are repeated identically. That might reduce the time
%         needed for loading the commands between tReps (CLTime).
%
%
%
% OUTPUT:
%
%   raw_data_AQs
%         2d double (or single) matrix with the (averaged) raw data. The size of
%         the array is: max. nSamples x nAQs (number of AQ windows in all tReps)
%
%
%   Seq
%         Same as the input structure Seq. Some of the values might be changed
%         to the actually used ones. Additionally, some fields are added.
%         Amongst the added fields are:
%
%     SequenceTime
%           Total duration of the sequence in seconds starting at the tOffset of
%           the first tRep of the first averaging step and ending at the end of
%           the last tRep of the last averaging step.
%
%     HW
%           Structure that retains the values set in the HW object or structure
%           as they were used for the measurement.
%
%
%   data_S_AQs_TRs
%         Structure with the data sorted in 3d arrays with the dimensions
%         [nSamples per AQ window x nAQs per tRep x ntReps]
%         See the function "get_data" for further details.
%
%
%   data_1D
%         Structure with the data sorted into column vectors where the samples
%         from different acquisition windows are separated by NaN. See the
%         function "get_data_1D" for further details.
%
%
%   data_Cell
%         For systems with multiple simultaneous acquisition frequencies
%         (HW.RX.n > 1), the data for each frequency is returned in the format
%         of raw_data_AQs in separate cells. For most system setups, this output
%         is empty and can be ignored.
%
%
%   averages
%         2d double matrix with the same layout as raw_data_AQs with the number
%         of averages received for each sample (so far).
%
%
%
% LIST OF ADDITIONAL IMPORTANT FIELDS:
%
%% Seq parameters
%     Seq.CorrectAQWindowPhase
%     Seq.CommandWindowOutput
%% Plot
%     Seq.PlotCurrentSquareTime
%     Seq.PlotGradCoilTemperature
%% Seq output
%         Seq.CLTime_ok
%         Seq.CommandsPerTrep
%         Seq.GradPowerAll
%         Seq.GradTimeAll
%         Seq.GradCoilMeanPower
%         Seq.GradCoilPowerDissipation
%         Seq.GradCoilThermalChange
%         Seq.GradCoilTemperatur
%% Other structures
%%    Seq.AQSlice
%     nRead
%     nPhase
%     HzPerPixMin
%     sizeRead
%     sizePhase
%     thickness
%     excitationPulse: @Pulse_...
%     inversionPulse: @Pulse_...
%     PhaseOS
%     TurboFactor
%     alfa
%     phi
%     theta
%     plotkSpace
%     plotImage
%     plotPhase
%     AmplitudeUnit
%     AmplitudeUnitScale
%     ReadCoordinate
%     PhaseCoordinate
%     SliceCoordinate
%     SliceCoordinateInvert
%     plotImageHandle
%     UseAQWindow
%     thicknessInversion
%     UsetRep
%     ReadOS
%     ReadOSUsedForImage
%     sizePhaseSpoil
%     SpoilCoordinate
%     TurboBreak
%     Center2OriginImage
%     sizePhase3D
%     CenterRot
%     AreaCoil
%     nPhase3D
%     VoxelVolume
%%    HW.
%         RootPath
%         MagnetShimPath
%         TX: [1x1 struct]
%         NetworkCalPath
%         MMRT: [1x1 struct]
%         Function_Prepare_Measurement
%         Function_Before_Measurement
%         Function_While_Measurement
%         Function_After_Measurement
%         Function_Abort_Measurement
%         Function_Run_Drawnow_While_Measurement
%         FindFrequencyPlot
%         FindFrequencyPause
%         FindFrequencyPauseFID
%         FindFrequencySweep: [1x1 struct]
%         Network: [1x1 struct]
%         FindShim: [1x1 struct]
%         ReInit
%         Grad: [1x1 struct]
%         DICOM: [1x1 struct]
%         Gamma: [1x1 struct]
%         GammaDef
%         GammaX
%         FindFrequencyGamma
%         Constant: [1x1 struct]
%         RelativCableSpeed
%         AddCableLength
%         TempR50
%         RX: [1x1 struct]
%         DigitalIO: [1x1 struct]
%         tRepInit
%         B0
%         MagnetShim: [1x4 double]
%         TX2RXdeadTime
%         RX2TXdeadTime
%         NetworkCalPathStandard
%         fLarmor
%         NetworkCal: [1x1 struct]
%         tFlip180Def
%         tFlip90Def
%         tFlip1Def
%         tFlip180min
%         tFlip90min
%         tFlip1min
%           LOOP_COUNT_END
%         RX.WarningNonValidData
%         RX.NonValidData
%%    AQ.
%         Start
%         nSamples
%         fSample
%         Frequency
%         Phase
%         ResetPhases
%         Gain
%         Repeat
%         AutoPowerSyncFrequency
%         PowerSyncFrequency
%           Amplitude2Norm
%           Amplitude2VMmrt
%           Amplitude2Vadc
%           Amplitude2Uin
%           Amplitude2Adc
%           Norm2Amplitude
%           VGAGain
%           PreamplifierVoltage
%           setPreamplifier
%           rawData2Norm
%           Latenz
%           AmplitudeCalGain
%           PhaseCal
%%    TX.
%         Channel
%         Start
%         Frequency
%         Duration
%         Amplitude
%         Phase
%         BlankOffset
%         BlankPostset
%         Repeat
%         setSwitchMode
%         setNetworkMode
%         setNetworkForward
%           StartEmpty
%           AmplitudeCorr
%           PaUout2AmplitudeCalibrationGain
%           PaUout2AmplitudeCalibrationPhaseOffset
%           PaUoutCalibrated
%           PaUoutCalibratedPhase
%           Uout2PaUoutCalibrationGain
%           Uout2PaUoutCalibrationPhaseOffset
%           UoutCalibrated
%           UoutCalibratedPhase
%           MmrtUout2UoutCalibrationGain
%           MmrtUout2UoutCalibrationPhaseOffset
%           MmrtUoutCalibrated
%           MmrtUoutCalibratedPhase
%           Norm2MmrtUoutCalibrationGain
%           Norm2MmrtUoutCalibrationPhaseOffset
%           NormCalibrated
%           NormCalibratedPhase
%           PhaseCal
%           DAC
%           RelativeDACQuantizationError
%           TxBits
%           AmpNorm
%           Latenz
%           TXBlankAfterLatenz
%           MIXmod
%%    Grad.
%         Time
%         Amp
%         Shim
%         Repeat
%           Current
%           CurrentTime
%           Power
%           MeanCurrentSquare
%           CurrentSquareDissipation
%           CurrentSquareTimeChange
%           PowerOk
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

end


%#function CIC_corr
%#function GetRawData_PPGfast
%#function PD.IgnoreWarning
%#function PD.MRISequence
%#function PD.TXMaxDef
%#function PD.Talker
%#function Poll_PPGfast
%#function Start_PPGfast
%#function add_DigitalIO
%#function check_CoilTemperature
%#function check_GradientFuse
%#function create_PPGfast
%#function get_PowerSyncFrequency
%#function hex2uint64
%#function isemptyfield
%#function plotSeq
%#function plotSeqTR
%#function set_AQfast
%#function set_DigitalIO
%#function set_Gradfast
%#function set_HFfast
%#function sleep

