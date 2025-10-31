classdef MMRT < handleHidden
%% MMRT  Class storing properties of the MMRT
%
% constructor:
%   MMRTobj = MMRT(HW)
%
% Input:
%   HW      object of class PD.HW
%
% MMRT properties:
%   fSystem                 clock frequency of system
%   UsbDelay                delay to start the usb transfer (PC to MMRT)
%                           (down to 0 because of CommandCacheSize Buffer)
%   UsbSpeed                USB speed in Bytes/sec (up to 40 MB/s)
%   TimeToSequenceDelay     reduction of command loadtime between sequences
%   CommandCacheSize        number of commands that can be buffered in FPGA
%   fPowerSync              DC/DC failure, optimum frequency and range
%   fPowerSyncMin
%   fPowerSyncMax
%   fPowerSyncDist          Ziel AQFrequenz Position zwischen Oberwellen des
%                           DC/DCs
%   PowerSyncFrequencyList  ???
%
%   Command                 List of commands for the MMRT device.
%
%   initializeOnLoad        send "initialization packet" with LoadHW
%
%   FPGA_Libs               version of FPGA libraries
%   FPGA_Firmware           version of FPGA firmware
%   OptLoadMySystemStr      optional string added to a new LoadMySystem file
%
%   DeviceSerialStr         Serial number of device as a string
%   DeviceSerial            Serial number of device
%
%
% See also:
%   PD.HW
%
%
% ----------------------------------------------------------------------------
% (C) Copyright 2016-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ----------------------------------------------------------------------------

end


%#function PD.Commands
%#function PD.HW
%#function handleHidden

