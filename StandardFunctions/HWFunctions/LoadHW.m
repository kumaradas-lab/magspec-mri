function [HW, mySave, ShimMatrix, DDS] = LoadHW(varargin)
%% Connect to MRT device and load system standards
%
%   [HW, mySave, ShimMatrix] = LoadHW('HW', HW, 'mySave', mySave, ...
%          'ShimMatrix', ShimMatrix, 'doForce', doForce, 'oop', useOop, ...
%          'initOnce', initOnce, 'userName', userName)
%
% When called, the connection to the MMRT device is established and default
% values are set. Default values can be changed by editing "LoadMySystem.m".
%
%
% INPUT:
%
%   The following optional property-value-pairs are allowed:
%   'HW'          HW-structure or object. Default: []
%   'mySave'      mySave structure. Default: []
%   'ShimMatrix'  ShimMatrix object. Default: []
%   'doForce'     logical or string. If set to true, HW (and ShimMatrix) objects
%                 are re-created even if they were already passed as valid input
%                 arguments. The strings 'true', 'force', 'reload' or 'init' are
%                 treated as true. Everything else is treated as false.
%                 Default: false.
%   'oop'         logical. If set to false, a structure is created. Otherwise,
%                 HW is an object of class PD.HWClass.
%   'initOnce'    logical. If set to true, HW.Grad.HoldShim and
%                 HW.MMRT.initializeOnLoad are set to true. If set to false,
%                 HW.MMRT.initializeOnLoad is set to false. If empty, the
%                 default values are kept. Default: []
%   'userName'    string. User name that is used for different User folders.
%                 Default: 'default'
%
%
% OUTPUT:
%
%   HW            HW-structure or object
%   mySave        mySave-stucture
%   ShimMatrix    ShimMatrix object if applicable, empty otherwise
%   DDS           DDS object if applicable, empty otherwise
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function InitializeMMRT
%#function LoadCalcHW
%#function LoadCoil
%#function LoadGradAmp
%#function LoadLibs
%#function LoadRXTX_Cal
%#function LoadShimMatrix_Standard
%#function LoadSystem
%#function LoadSystem_Specific
%#function LoadSystem_Standard
%#function PD.HWClass
%#function PD.MagnetReadyState
%#function PD.MySave
%#function PD.UnwindProtectGuard
%#function PD.helper_HW_Shim
%#function checkDeviceSerial
%#function check_FPGA_Firmware_Libs
%#function getOpenMatlabRootPath
%#function get_MRDevice
%#function initialize_MRDevice
%#function isAbsolutePath
%#function isemptyfield
%#function openMatlabRevision
%#function setOpenMatlabSearchPath
%#function set_FanDutyCycle
%#function set_system_frequency
%#function sleep
%#function use_ExternalClock

