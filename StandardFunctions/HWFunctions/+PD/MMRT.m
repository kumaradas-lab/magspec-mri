classdef MMRT < handleHidden
%% MMRT  Class storing properties of the MMRT
%
%       MMRTobj = PD.MMRT(HW)
%
% This object manages properties of connected mobile MRT device(s).
% A user shouldn't manually create an object of this type. Instead it will be
% created automatically as a property of "HW" while executing "LoadSystem".
%
%
% INPUT:
%
%   HW
%       Object of class PD.HWClass
%
%
% PROPERTIES:
%
%   fSystem
%       Clock frequency of system in Hertz. Don't change this property
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
%   UsbDelay
%       Delay to start the usb transfer (PC to MMRT) in seconds. It can be 0
%       because of CommandCacheSize buffer. But it must have a non-negative
%       value.
%
%   UsbSpeed
%       USB speed in Bytes/sec (up to 40 MB/s). This is used to calculate the
%       minimum duration of the CLTime between tRep.
%
%   TimeToSequenceDelay
%       Reduction of command loadtime between sequences in seconds. When
%       timing subsequent measurement sequences to each other, the CLTime
%       between the sequences is reduced by this duration. An end user
%       shouldn't need to change the default value.
%
%   CommandCacheSize
%       Number of commands that can be buffered in the FPGA. This is used to
%       calculate the minimum duration of the CLTime between tReps. An end
%       user shouldn't need to change the default value.
%
%   fPowerSync
%       Optimum frequency for DC/DC converter in Hertz.
%
%   fPowerSyncMin
%       Minimum frequency for DC/DC converter in Hertz.
%
%   fPowerSyncMax
%       Maximum frequency for DC/DC converter in Hertz.
%
%   fPowerSyncDist
%       Relative distance of AQ frequency between peaks from harmonics of
%       DC/DC converter.
%
%   PowerSyncFrequencyList
%       Vector with a list of possible frequencies for DC/DC converter in
%       Hertz.
%
%   FanDutyCycleMin
%       Minimum fan duty cycle in percent.
%
%   actEndBitWidth
%       Number of bits used to count ticks in a tRep. Used to check the
%       maximum possible duration of a tRep. Changing it from its default
%       value can lead to undefined behavior.
%
%   CLTimeBitWidth
%       Number of bits used to count ticks in the CLTime between tReps. Used
%       to check the maximum possible duration of a CLTime. Changing it from
%       its default value can lead to undefined behavior.
%
%   TriggerDebounceTime
%       Default trigger debounce time in seconds. This might require a
%       non-standard Firmware to have an effect.
%
%   Command
%       Object with a list of commands for the MMRT device. Changing it from
%       its default values can lead to undefined behavior.
%
%   ReduceCommands
%       Boolean value to indicate if commands that don't change from the
%       previous tRep can be removed. This allows for possible shorter CLTimes
%       between tReps.
%
%   ClockSyncEnable
%       Boolean value to indicate whether the internal clocks should be
%       synchronized. It must be set to 0 to change clock divisors. Instead of
%       changing this value manually, use "set_system_frequency".
%
%   useExtClock10MHz
%       Boolean value to switch between internal clocks (false) or an external
%       10 MHz clock at Clk port of the device (true).
%
%   CommunicationType
%       Type of the talker that is used for the connection to the MMRT device.
%       Currently, the following talker types are supported:
%           0: .NET assembly (default)
%           1: dummy (no actual connection to an MMRT device)
%           2: MATLAB socket connection (experimental)
%           3: .NET socket connection (experimental)
%
%   FPGA_Libs
%       Version number of the currently loaded .NET assembly libraries.
%       Changing it from its set value in places other than the LoadMySystem.m
%       file can lead to undefined behavior.
%
%   TCP_Server
%       IP address of the server that is connected to the MMRT. Changing it
%       from its set value in places other than the LoadMySystem.m file can
%       lead to undefined behavior. (Default: 'localhost')
%
%   TCP_Port
%       TCP port of the socket server that is connected to the MMRT. Changing
%       it from its set value in places other than the LoadMySystem.m file can
%       lead to undefined behavior. (Default: 8887)
%
%   FPGA_Firmware
%       Version number of the currently loaded FPGA firmware. Changing it from
%       its set value in places other than the LoadMySystem.m file can lead to
%       undefined behavior.
%
%   type
%       String with the type of the connected device. Currently this can only
%       be 'drive-l'.  Changing it from its set value can lead to undefined
%       behavior.
%
%   OptLoadMySystemStr
%       Optional string that is added to a new LoadMySystem file. This is used
%       only internally.
%
%   initializeOnLoad
%       Send initialization packet during LoadSystem. Among other things this
%       initialization package sets the currently calibrated shim values for
%       the console.
%
%   initialized
%       Boolean value to indicate whether the device was successfully
%       initialized since Matlab was started. Changing it from its set value
%       can lead to undefined behavior.
%
%   DeviceSerialStr
%       Serial number of the connected console as a string.
%
%   DeviceSerial
%       Numeric part of the serial number of the connected console.
%
%   OpenMatlabRevision
%       Revision number of the currently used version of OpenMatlab.
%
%   OpenMatlabBuildTime
%       Build time of the currently used version of OpenMatlab in Matlab's
%       "now" format.
%
%
% See also:
%   PD.HWClass
%
%
% ----------------------------------------------------------------------------
% (C) Copyright 2016-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ----------------------------------------------------------------------------

end


%#function PD.Commands
%#function PD.HWClass
%#function PD.RX
%#function PD.TXClass
%#function PD.Talker
%#function handleHidden

