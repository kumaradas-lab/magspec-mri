classdef RX < handleHidden
%% Class storing properties for the receive channel
%
%   RXobj = PD.RX(HW)
%
% This object manages receiption properties of connected mobile MRT device(s).
% A user shouldn't manually create an object of this type. Instead it will be
% created automatically as a property of "HW" while executing "LoadSystem".
%
%
% INPUT:
%   HW
%       Object of class PD.HWClass
%
%
% PROPERTIES:
%
%   fSample
%       Sampling frequency of the RX unit in Hertz. Don't change this property
%       manually. Instead use the function "set_system_frequency".
%
%   ClockDivisor
%       Integer scalar with a divisor with respect to the reference clock of
%       the mobile MRT (1 GHz). Don't change this property manually. Instead
%       use the function "set_system_frequency".
%
%   ClockDelay
%       Delay for the system clock in seconds. An end user shouldn't need to
%       change this property. Changing it from its default value can lead to
%       undefined behavior.
%
%   DdsPicBits
%       Bit width of the DDS PIC register. Changing it from its default value
%       can lead to undefined behavior.
%
%   DdsPofBits
%       Bit width of the DDS POF register. Changing it from its default value
%       can lead to undefined behavior.
%
%   ADCBits
%       Bit width of the RX ADC. Changing it from its default value can lead
%       to undefined behavior.
%
%   AdcUinMax
%       Maximum input voltage of ADC in Volts. Changing it from its default
%       value can lead to undefined behavior.
%
%   AdcUin2Norm
%       Conversion factor from input voltage of the ADC to its normalized
%       output in 1/Volts. Changing it from its default value can lead to
%       undefined behavior.
%
%   Norm2Adc
%       Conversion factor from normalized output to digitized output of ADC.
%       Changing it from its default value can lead to undefined behavior.
%
%   VgaDacBits
%       Number of bits for the control voltage of the variable gain amplifier.
%       Changing it from its default value can lead to undefined behavior.
%
%   VgaDacUoutMax
%       Maximum for the control voltage of the variable gain amplifier in
%       Volts. Changing it from its default value can lead to undefined
%       behavior.
%
%   MmrtUin2AdcUin
%       Linear amplification factor from MMRT input voltage to the input
%       voltage of the ADC. Changing it from its default value can lead to
%       undefined behavior.
%
%   NormWarnSaturation
%       If the normalized amplitude exceeds this level, a warning is emitted.
%
%   CIC_Decimation_Min
%       Minimum decimation factor of the CIC filter.
%
%   CIC_Decimation_Max
%       Maximum decimation factor of the CIC filter.
%
%   CIC_N
%       Number of stages of the cascaded integrator-comb (CIC) filter.
%       Changing it from its default value can lead to undefined behavior.
%
%   CIC_M
%       Number of samples per stage of the cascaded integrator-comb (CIC)
%       filter. Changing it from its default value can lead to undefined
%       behavior.
%
%   n
%       Number of RX inputs of the connected device. Changing it from its
%       default value can lead to undefined behavior.
%
%   nSampleLatenz
%       Latency of acquisition in number of decimated samples.
%
%   nSampleRXLatenz
%       Latency of acquisition in number of system frequency samples. Changing
%       it from its default value can lead to undefined behavior.
%
%   nSampleLatenzRaw
%       Latency of acquisition in number of system frequency samples in raw
%       data mode. Changing it from its default value can lead to undefined
%       behavior.
%
%   CutSamples
%       Number of decimated samples that are cut off from each start of
%       acquistion windows. Changing it from its default value can lead to
%       undefined behavior.
%
%   nSamplesExtra
%       Number of decimated samples that are additionally acquired for each
%       acquisition window (due to HW.RX.CutSamples). Changing it from its
%       default value can lead to undefined behavior.
%
%   nSamplesRXExtra
%       Number of system frequency samples that are additionally acquired for
%       each acquisition window. Changing it from its default value can lead
%       to undefined behavior.
%
%   nSamplesExtraCIC
%       Number of decimated samples that are additionally acquired to
%       "initiate" the CIC filter. Changing it from its default value can lead
%       to undefined behavior.
%
%   nSamplesRXExtraCIC
%       Number of system frequency samples that are additionally acquired to
%       "initiate" the CIC filter. Changing it from its default value can lead
%       to undefined behavior.
%
%   nSamplesExtraSamplingFactor
%       Number of software downsampled samples that are additionally acquired
%       for each acquisition window. This only applies if software
%       downsampling is used. Changing it from its default value can lead to
%       undefined behavior.
%
%   nSamplesPhaseLatency
%       Number of decimated samples that are used to correct the phase of the
%       acquired signal due to latency. Changing it from its default value can
%       lead to undefined behavior.
%
%   fSampleName
%       String with a (short) description that is used in plots that show the
%       sampling frequency of the acquisition windows.
%
%   fSampleUnit
%       String with the unit symbol for sampling frequency used in plots.
%
%   fSampleUnitScale
%       Scalar value with the scaling factor that is used for the sampling
%       frequency in plots.
%
%   NonValidData
%       Scalar value that is used to replace ivalid data in the array with
%       the acquired signal.
%
%   WarningNonValidData
%       Boolean value to indicate if a warning should be emitted in case
%       invalid data was received. That could happen if the USB-connection
%       between the console and the PC was too slow to transfer all data
%       before the buffer in the console overflowed.
%
%   Uin2VgaUin
%       Conversion factor from voltage at the TRx port to the input of the
%       variable gain amplifier in V/V. Changing it from its default value can
%       lead to undefined behavior.
%
%   Gain2VGAVolts
%       Function handle with the prototype VGAVolts = @(Gain, HW). This
%       function is used to calculate the VGA control voltage necessary to set
%       it to the requested gain.
%
%   VGAGainMin
%       Minimum gain of the used variable gain amplifier. Changing it from its
%       default value can lead to undefined behavior.
%
%   VGAGainMax
%       Maximum gain of the used variable gain amplifier. Changing it from its
%       default value can lead to undefined behavior.
%
%   VGAGainDef
%       Default gain of the used variable gain amplifier.
%
%   Rin
%       Input impedance of the TRx port in Ohm. Changing it from its default
%       value can lead to undefined behavior.
%
%   VgaUout2MmrtUin
%       Conversion factor from the output of the variable gain amplifier to
%       the MMRT input in V/V. Changing it from its default value can lead to
%       undefined behavior.
%
%   LnaSN
%       Serial number of the currently loaded low-noise amplifier (LNA). This
%       is an informational value only.
%
%   DataClassSingle
%       Boolean value to indicate whether the received signal is converted to
%       double precision floating-point numbers (false) or single precision
%       floating-point numbers (true). Using single precision floating-point
%       number can reduce the amount of memory needed to hold the received
%       data. However, not all Matlab operations are available for single
%       precision input.
%
%   DataMode
%       Integer indicating the default mode for the transfer of the data
%       between the console and the PC. The default mode can be overridden by
%       setting AQ.DataMode. Depending on the loaded FPGA firmware the
%       following data modes are available:
%         0: 2x48 bits for the real and imaginary parts of each decimated
%            sample of a single mixer.
%         1: two mixers (for dual-frequency acquisitions) with 2x24 bits for
%            the real and imaginary parts of each decimated sample of each
%            mixer.
%         2: 2x24 bits for the real and imaginary parts of each decimated
%            sample of a single mixer. This allows to approximately double the
%            maximum sampling rate compared to data mode 0.
%         3: 1x24 bits for the real part of each decimated sample of a single
%            mixer (discarding the imaginary part). This allows to
%            approximately quadruple the maximum sampling rate compared to
%            data mode 0.
%
%   ClampCoil
%       Structure with settings for clamping during acquisitions.
%
%   Amplitude2LnaUin
%       Receive coil efficiency at the default Larmor frequency (HW.fLarmor)
%       in V/T. If empty, a default of 1./HW.TX.PaUout2Amplitude is used.
%
%   Amplitude2LnaUinX
%       Receive coil efficiency at the Larmor frequency of the X-nucleus
%       (HW.fLarmorX) in V/T. If empty, a default of
%       1./HW.TX.PaUout2AmplitudeX is used.
%
%   AreaCoil
%       Cross section of the receive coil. This is used to calculate the
%       magnetic field density per voxel in images.
%
%   EffectiveCoilLength
%       Effective length of the receive coil. This is used to calculate voxel
%       volumes if there is no encoding in y-direction.
%
%   RiseTime
%       Rise time of the oscillation in the RX coil in seconds. This is used
%       to ensure a sufficiently large pause between excitation and
%       acquisition in some sequences.
%
%   LNAGain
%       Gain of an (optionally used) low-noise amplifier in V/V.
%
%   AmplitudeName
%       String with a (short) description that is used in plots that show the
%       the received amplitude.
%
%   AmplitudeUnit
%       String with the unit symbol for the amplitude of the received signal.
%
%   AmplitudeUnitScale
%       Scalar value with the scaling factor that is used for the amplitude of
%       the received signal.
%
%   CoilFrequencyResponse
%       Structure with fields for the frequency response of the RX coil.
%
%   FrequencyDivisor
%       Default frequency divisor for the frequency of the mixer.
%
%   Calibration
%       Structure with calibration properties of the signal amplifiers.
%
%   Amplitude2Norm
%       Read-only property with the conversion factor from the received
%       amplitude in the rf coil to the normalized output of the ADC at the
%       default Larmor frequency (HW.fLarmor). Consider using
%       "get_RX_Amplitude" or the corresponding field in the "data" structure
%       instead.
%
%   Amplitude2Uin
%       Read-only property with the conversion factor from the received
%       amplitude in the rf coil to voltage at the TRx port at the default
%       Larmor frequency (HW.fLarmor). Consider using "get_RX_Amplitude"
%       instead.
%
%   Amplitude2Adc
%       Read-only property with the conversion factor from the received
%       amplitude in the rf coil to digitized output values of the ADC at
%       the default Larmor frequency (HW.fLarmor). Consider using
%       "get_RX_Amplitude" instead.
%
%   Norm2Amplitude
%       Read-only property with the value 1./HW.RX.Amplitude2Norm.
%
%   GainMin
%       Read-only property with the minimum gain from the amplitude in the rf
%       coil to the normalized output of the ADC in 1/T at the default Larmor
%       frequency. The range is determined by the range of the VGA.
%
%   GainMax
%       Read-only property with the maximum gain from the amplitude in the rf
%       coil to the normalized output of the ADC in 1/T at the default Larmor
%       frequency. The range is determined by the range of the VGA.
%
%   GainDef
%       Read-only property with the default gain from the amplitude in the rf
%       coil to the normalized output of the ADC in 1/T at the default Larmor
%       frequency. The range is determined by the range of the VGA.
%
%
% ----------------------------------------------------------------------------
% (C) Copyright 2016-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ----------------------------------------------------------------------------

end


%#function PD.HWClass
%#function PD.Talker
%#function handleHidden

