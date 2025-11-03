function [raw_data_AQs, Seq, data_S_AQs_TRs, data_1D, data_Cell, averages] = set_sequence(HW, Seq, AQ, TX, Grad)
%% Send pulse program to MMRT and receive data
%
%   [raw_data_AQs, Seq, data_S_AQs_TRs, data_1D, data_Cell, averages] = set_sequence(HW, Seq, AQ, TX, Grad)
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
%         Structure or PD.HW object as created when executing LoadSystem.
%
%
%   Seq
%         Structure with pulse program settings. Amongst others, the following
%         fields can be used (default values apply if the fields are omitted or
%         empty):
%
%     tRep
%           Row vector with the duration of the tReps in seconds.
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
%
%   AQ, TX, Grad
%         Structures defining the acquisition windows, rf pulses and gradient
%         pulses in the pulse program. See the manual for further details.
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
% (C) Copyright 2012-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------

end


%#function CIC_corr
%#function GetRawData_PPGfast
%#function PD.Talker
%#function Poll_PPGfast
%#function Start_PPGfast
%#function add_DigitalIO
%#function check_CoilTemperature
%#function check_GradientFuse
%#function create_PPGfast
%#function hex2uint64
%#function isemptyfield
%#function plotSeq
%#function plotSeqTR
%#function set_AQfast
%#function set_DigitalIO
%#function set_Gradfast
%#function set_HFfast
%#function sleep

