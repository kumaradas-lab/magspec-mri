classdef GradClass < handleHidden
%% Gradient system configuration
%
%     GradObj = PD.GradClass()
%
% This object manages gradient amplifier and gradient system properties for
% connected mobile MRT device(s). A user shouldn't manually create an object
% of this type. Instead it will be created automatically as a property of "HW"
% while executing "LoadSystem".
%
%
% PROPERTIES:
%
%   TimeDelay
%       Time delay of gradient pulses in seconds (e.g., due to induction
%       induced counter voltage or eddy currents).
%
%   PaEnable
%       Boolean value to un-mute the gradient amplifiers.
%
%   TubeDiameter
%       Default value for the diameter of a used sample in meters. This value
%       is used by some sequences for default image sizes.
%
%   SliceTimeDelay
%       Default time delay for slice selective gradient pulses in seconds.
%       This value is used by some sequences as the default time delay of
%       slice selective gradient pulses.
%
%   ReadTimeDelay
%       Default time delay for readout encoding gradient pulses in seconds.
%       This value is used by some sequences as the default time delay of
%       readout encoding gradient pulses.
%
%   PhaseTimeDelay
%       Default time delay for phase encoding gradient pulses in seconds.
%       This value is used by some sequences as the default time delay of
%       phase encoding gradient pulses.
%
%   Name
%       1x4 cell array of strings that are used in some plots of gradient
%       pulses.
%
%   AmpUnit
%       1x4 cell array of strings with the amplitude units used in some plots
%       of gradient pulses.
%
%   AmpUnitScale
%       1x4 vector with the scale factors for the amplitude units used in some
%       plots of gradient pulses.
%
%   LengthUnit
%       1x4 cell array of strings with the corresponding length units for the
%       gradient encoding direction used in some imaging plots.
%
%   LengthUnitScale
%       1x4 vector with the scale factors for the length units used in some
%       imaging plots.
%
%   Inductance
%       1x4 vector with the self inductance of the connected gradient coils in
%       Henry. This is used to check if gradient rise times can be reached
%       with the power of the used gradient amplifier.
%
%   ShimGradients
%       1x4 vector with Boolean values indicating which gradients are used to
%       shim the magnet.
%
%   CombineCurrentOutputs
%       Matrix with indices of output channels where each row indicates which
%       channels are connected to the same gradient coil. By default, this is
%       empty meaning that each gradient channel is connected to a separate
%       gradient coil.
%
%   ImageVol
%       1x6 vector with the limits of the default image volume:
%           [xmin xmax ymin ymax zmin zmax]
%       This is used in some (imaging) sequences to determine the default
%       image size.
%
%   ImageVolOffset
%       1x3 vector with the x, y, and z coordinates of the image center in
%       meters. This is used in some (imaging) sequences to shift the default
%       image volume.
%
%   SystemTimeDelay
%       1x4 vector with time delays of each gradient channel in seconds.
%
%   tEC
%       Scalar value with the time between gradient ramps and subsequent
%       acquisitions or rf pulses in seconds. This time is used in some
%       (imaging) sequences to make sure potential eddy currents have decayed
%       sufficiently for the next element in the pulse program.
%
%   CoilThermalGroup
%       Before a pulse sequence executes, model calculations are applied to
%       check if the duty cycle doesn't damage the gradient system. This 1x4
%       vector contains indices that assign each gradient channel to a group
%       in these model calculations.
%
%   CoilPowerDissipation
%       1xn vector with a value of the power dissipation in Watts of each
%       thermal group.
%
%   CoilTemperatur
%       1xn vector with the temperature in degrees Celsius of the coils in
%       each group at the start of a measurement. Additionally, this is used
%       as the temperature of an "infinite reservoir" to which the heat
%       transfers thoughout the measurement.
%
%   CoilMaxTemperature
%       Maximum temperature of the coil in degrees Celsius. A model is used
%       for which no heat is transferred if the coil is at
%       HW.Grad.CoilTemperatur and HW.CoilPowerDissipation is transferred at
%       this temperature. At intermediate temperature, a linear model for the
%       heat transfer is used. If the model exceeds this temperature at some
%       point during the measurement, the measurement is not started. (See
%       also "check_CoilTemperature".)
%
%   CoilThermalCapacity
%       1xn vector with the thermal capacity of each coil group in J/K.
%
%   CoilMaxDcCurrent
%       1x4 vector with the maximum (nominal) DC current of the fuse for each
%       gradient channel in Ampere.
%
%   CoilCurrentSquareTime
%       1x4 vector with the time-lag fuse parameter (I^2 t) of the fuse for
%       each gradient channel in A^2*s.
%
%   PaUoutMax
%       1x4 vector with the maximum power of each channel of the gradient
%       amplifer in Watts.
%
%   PaUoutMin
%       1x4 vector with the minimum power of each channel of the gradient
%       amplifer in Watts. (This is useful if a power amplifier is used that
%       has assymetric power for positive and negative currents.)
%
%   PaName
%       String with an identifier for a connected external gradient amplifier.
%       This is informational only.
%
%   PaNameList
%       Cell array of strings with a list of (potentially) connected gradient
%       amplifiers.
%
%   LoadPaPath
%       String with the name of a script in the Matlab search path with
%       calibration values of the connected external gradient amplifier.
%
%   MMRTTimeOffset
%       1x4 vector with time delays of each gradient channel in seconds due to
%       the MMRT implementation.
%
%   HoldShim
%       Boolean value to indicate whether the gradient amplifier should be
%       muted after a measurement is finished (false), or whether it should be
%       kept on (true).
%
%   HoldShimNormMax
%       1x4 vector with relative values of the maximum gradient amplitude at
%       each gradient (xyzB).
%
%   PowerDown
%       Boolean value to turn off the gradient amplifiers.
%
%   PaPmaxInt
%       1x4 vector with maximum internal power dissipation in W for each
%       gradient amplifier channel. Informational only - unused.
%
%   tRamp
%       Default ramp time of gradient pulses in seconds. This is used in most
%       (imaging) sequences.
%
%   Status1
%       Boolean value to indicate whether the "status 1" signal (power supply
%       ok) of the gradient amplifier should be checked before starting a
%       measurement.
%
%   Status2
%       Boolean value to indicate whether the "status 1" signal (over-voltage
%       or over-temperature ok) of the gradient amplifier should be checked
%       before starting a measurement.
%
%   ExtGradSN
%       Serial number of a connected external gradient amplifier. This is
%       informational only.
%
%   n
%       Read-only property with the number of available gradient channels.
%
%   x
%       Index of the gradient amplifier channel that is connected to the
%       x-gradient coils.
%
%   y
%       Index of the gradient amplifier channel that is connected to the
%       y-gradient coils.
%
%   z
%       Index of the gradient amplifier channel that is connected to the
%       z-gradient coils.
%
%   B
%       1xn vector with indices of additional gradient amplifier channel(s).
%       The number of available additional channels depends on the
%       implementation. If the number of elements in this property is
%       different from 1, properties that are documented to contain 4 elements
%       will contain more elements to accomodate for all available gradient
%       channels. Changing the number of elements in this property manually
%       leads to undefined behavior.
%
%   xyzBDir
%       1x4 vector with plus or minus 1 to indicate the sign of the connected
%       gradient coils (invert the current for the coils with -1).
%
%   LoadRin
%       1x4 vector with the resistance of the gradient coils. This is used if
%       the gradient amplifiers are voltage controlled.
%
%   LoadIin2Amp
%       1x4 vector with the gradient coil efficiencies in (T/m)/A. It is only
%       necessary to set either HW.Grad.LoadIin2Amp or HW.Grad.LoadUin2Amp.
%
%   LoadUin2Amp
%       1x4 vector with the gradient coil efficiencies in (T/m)/V. It is only
%       necessary to set either HW.Grad.LoadIin2Amp or HW.Grad.LoadUin2Amp.
%
%   PaCurrentControlled
%       1x4 vector with Boolean values to indicate if the corresponding
%       gradient amplifier channel is current controlled (true) or voltage
%       controlled (false). Setting this incorrectly could lead to undefined
%       behavior.
%
%   AmpCurrentDependent
%       1x4 vector with Boolean values to indicate whether the output
%       amplitudes of the respective gradient amplifier channel are actually
%       current dependent (true, e.g., for gradient coils) or voltage
%       dependent (false, e.g., for Piezo elements).
%
%   PaRin
%       1x4 vector with input impedance of the respective gradient amplifier
%       channel in Ohm (important if smaller than 1 kOhm).
%
%   PaOffsetU
%       1x4 vector with the amplifier output offset of the respective channel
%       in Volts. This is only valid for channels that are set to
%       HW.Grad.PaCurrentControlled(iChannel)=false. Setting a non-zero value
%       for other channels might lead to undefined behavior.
%
%   PaOffsetI
%       1x4 vector with the amplifier output offset of the respective channel
%       in Amperes. This is only valid for channels that are set to
%       HW.Grad.PaCurrentControlled(iChannel)=true. Setting a non-zero value
%       for other channels might lead to undefined behavior.
%
%   PaUin2PaUout
%       1x4 vector with gain for each gradient amplifier channel in V/V. This
%       is only valid for channels that are set to
%       HW.Grad.PaCurrentControlled(iChannel)=false. Setting a non-zero value
%       for other channels might lead to undefined behavior.
%
%   PaUin2PaIout
%       1x4 vector with gain for each gradient amplifier channel in A/V. This
%       is only valid for channels that are set to
%       HW.Grad.PaCurrentControlled(iChannel)=true. Setting a non-zero value
%       for other channels might lead to undefined behavior.
%
%   PaRout
%       1x4 vector with output impedance of each gradient amplifier channel in
%       Ohm.
%
%   DacBits
%       1x4 vector with the number of bits in the input of the DAC for each
%       gradient amplifier channel. Changing this from the set value might
%       lead to undefined behavior.
%
%   ExtGradOffsetDAC
%       1x4 vector with the offset in DAC bits of the "Ext. Grad" connector of
%       the console.
%
%   Dac2ExtGradUout
%       1x4 vector with the voltage per 1 bit DAC input for each gradient
%       amplifier channel without load in Volts.
%
%   Dac2ExtGradUout500Ohm
%       1x4 vector with the voltage per 1 bit DAC input for each gradient
%       amplifier channel with 500 Ohm load in Volts.
%
%   OffsetinDACMin
%       1x4 vector with lower tolerance for HW.Grad.OffsetinDAC for each
%       gradient amplifier channel.
%
%   OffsetinDACMax
%       1x4 vector with upper tolerance for HW.Grad.OffsetinDAC for each
%       gradient amplifier channel.
%
%   MaxAmp
%       1x4 vector with maximum amplitude of each gradient coil in T/m (with
%       respect to the shimmed magnet).
%
%   MaxAmpSlice
%       Maximum amplitude for a slice encoding gradient in T/m. This is used
%       in some (imaging) sequences to set a different (lower) limit for slice
%       encoding gradients than for other gradient pulses.
%
%   MaxAmpTotal
%       1x4 vector (read-only) with maximum amplitude of each gradient coil in
%       T/m (with respect to the un-shimmed magnet).
%
%   AmpOffset
%       1x4 vector with the amplitude offset of each gradient coil in T/m.
%       This corresponds to the values in HW.MagnetShim for the respective
%       console.
%
%   AmpOffsetEstimated
%       1x4 vector with the estimated amplitude offset of each gradient coil
%       in T/m. These values are used by default as start values for the shim
%       amplitudes in "Find_Shim".
%
%   AmpOffsetExtra
%       1x4 vector with additional amplitude offsets for each gradient coil in
%       T/m. This acts the same as HW.Grad.AmpOffset. But it is handled
%       differently, e.g., when plotting pulse programs.
%
%   EddyCurrent
%       Structure with settings for compensating (uniform) Eddy currents.
%
%   CalibrationLoadIin
%       Structure with calibration values for a load dependent gain.
%
%   PaOffsetDAC
%       Read-only property. 1x4 vector with the offsets of the gradient
%       amplifier output for each gradient amplifier channel in DAC values.
%
%   Dac2Norm
%       Read-only property. 1x4 vector with the conversion factors from DAC
%       values to normalized output amplitude for each gradient amplifier
%       channel.
%
%   ExtGradRout
%       Read-only property. 1x4 vector with the output impedances of the "Ext.
%       Grad" connector for each gradient amplifier channel in Ohm.
%
%   xyzB
%       Read-only property. 1x4 vector with the indices of the gradient
%       amplifier channels used for each gradient coil.
%
%   Channel2xyzB
%       Read-only property. 1x4 vector with the gradient coil "index" that is
%       driven by each gradient amplifier channel. The gradient coil indices
%       correspond to the x, y, z, and B coils.
%
%   LoadUin2LoadIin
%       Read-only property. 1x4 vector with the conductance of each connected
%       load (i.e., gradient coil) in A/V.
%
%   Dac2PaUin
%       Read-only property. 1x4 vector with the conversion factors from DAC
%       values to input voltage at the gradient amplifier in V/bit for each
%       gradient amplifier channel.
%
%   PaCCindex
%       Read-only property. Row vector with the indices of the gradient
%       amplifier channels that are current controlled (see
%       HW.Grad.PaCurrentControlled).
%
%   PaVCindex
%       Read-only property. Row vector with the indices of the gradient
%       amplifier channels that are voltage controlled (see
%       HW.Grad.PaCurrentControlled).
%
%   AmpCDindex
%       Read-only property. Row vector with the indices of the gradient
%       amplifier channels that are current dependent (see
%       HW.Grad.AmpCurrentDependent).
%
%   AmpVDindex
%       Read-only property. Row vector with the indices of the gradient
%       amplifier channels that are voltage dependent (see
%       HW.Grad.AmpCurrentDependent).
%
%   Dac2PaUout
%       Read-only property. 1x4 vector with the conversion factors from DAC
%       input to the voltage at the gradient amplifier output in V/bit for
%       each gradient amplifier channel.
%
%   Dac2PaIout
%       Read-only property. 1x4 vector with the conversion factors from DAC
%       input to the current at the gradient amplifier output in A/bit for
%       each gradient amplifier channel.
%
%   Dac2LoadUin
%       Read-only property. 1x4 vector with the conversion factors from DAC
%       input to the voltage at the gradient coils in V/bit for each gradient
%       amplifier channel.
%
%   Dac2Amp
%       Read-only property. 1x4 vector with the conversion factors from DAC
%       input to the amplitude in the respective gradient coil in (T/m)/bit
%       for each gradient amplifier channel.
%
%   OffsetinDAC
%       Read-only property. 1x4 vector with the total offset in DAC bits (from
%       both the console and an external gradient amplifer) in bits for each
%       gradient amplifier channel.
%
%   Amp2PaUin
%       Read-only property. 1x4 vector with the conversion factors from the
%       input voltage of the gradient amplifier to the amplitude in the
%       gradient coils in V/(T/m) for each gradient amplifier channel.
%
%   PaUin2Amp
%       Read-only property. 1x4 vector with the conversion factors from the
%       amplitude in the gradient coils to the input voltage of the gradient
%       amplifier in (T/m)/V for each gradient amplifier channel.
%
%   Amp2ExtGradUout
%       Read-only property. 1x4 vector with the conversion factors from the
%       amplitude in the gradient coils to the output voltage at the "Ext.
%       Grad" port of the console in V/(T/m) for each gradient amplifier
%       channel.
%
%   ExtGradUout2Amp
%       Read-only property. 1x4 vector with the conversion factors from the
%       output voltage at the "Ext. Grad" port of the console to the amplitude
%       in the gradient coils in (T/m)/V for each gradient amplifier channel.
%
%   Amp2PaUout
%       Read-only property. 1x4 vector with the conversion factors from the
%       output voltage of the gradient amplifier to the amplitude in the
%       gradient coils in V/(T/m) for each gradient amplifier channel.
%
%   Amp2PaIout
%       Read-only property. 1x4 vector with the conversion factors from the
%       output current of the gradient amplifier to the amplitude in the
%       gradient coils in A/(T/m) for each gradient amplifier channel.
%
%   PaUout2Amp
%       Read-only property. 1x4 vector with the conversion factors from the
%       amplitude in the gradient coils to the output voltage of the gradient
%       amplifier in (T/m)/V for each gradient amplifier channel.
%
%   Amp2LoadUin
%       Read-only property. 1x4 vector with the conversion factors from the
%       amplitude in the gradient coils to the voltage at the gradient coils
%       in V/(T/m) for each gradient amplifier channel.
%
%   Amp2LoadIin
%       Read-only property. 1x4 vector with the conversion factors from the
%       amplitude in the gradient coils to the current though the gradient
%       coils in A/(T/m) for each gradient amplifier channel.
%
%
% ----------------------------------------------------------------------------
% (C) Copyright 2015-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ----------------------------------------------------------------------------

end


%#function LoadGradAmp
%#function PD.HWClass
%#function handleHidden

