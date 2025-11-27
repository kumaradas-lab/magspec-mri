function HW = InitializeMMRT(HW, iDevice, reInit)
%% Use "empty" sequence to initialize basic MMRT settings
%
%   HW = InitializeMMRT(HW, iDevice, reInit)
%
% If HW.MMRT.initializeOnLoad is true, an "empty" sequence is started to
% initialize basic MMRT settings. This is done each time "LoadSystem" is
% executed:
% Enables or disables the following (condition in brackets):
%   * soft mute of the gradient amplifier (HW.Grad.HoldShim)
%   * "manual" mute of the gradient amplifier (HW.Grad.PaEnable)
%   * "manual" sleep of the gradient amplifier (HW.Grad.PowerDown)
%   * pre-amplifier gain (HW.RX.VGAGainDef)
% Sets the following:
%   * ResetDataCounter
%   * ResetDDS
%   * ResetSystemTime
%   * SetStartPulseProgramOn
%
% The first time after the MMRT is initialized (noticeable by fan spooling up),
% other MMRT properties are set.
%
%
% INPUT:
%
%   HW
%       HW object or structure
%
%   iDevice
%       Index for the MMRT device
%       (Default: 1)
%
%   reInit
%       Initialize MMRT even if HW.MMRT.initializeOnLoad is false.
%       (Default: false)
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2017-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function PD.Commands
%#function PD.IgnoreWarning
%#function PD.MagnetBase
%#function PD.MagnetSampleHeater
%#function PD.SampleHeaterBase
%#function PD.Talker
%#function PD.UnwindProtectGuard
%#function get_COM_PDControl
%#function get_MagnetMacFromSerial
%#function get_SampleHeaterMacFromSerial
%#function isemptyfield

