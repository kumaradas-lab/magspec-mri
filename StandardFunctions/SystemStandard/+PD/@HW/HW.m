classdef HW < handleHidden
%% HW - Class with standard interface to Pure Devices ResearchLab
%
%   HW = PD.HW(rootPath)
%
% An instance of this class will be loaded when LoadSystem is executed.
% This is the main object holding a.o. default parameters and calibration
% settings of the MMRT. Some of these parameters are implementation specific
% and should not be changed by the user (unless they know what they are
% doing).
%
% A user should not create an object of this class. Instead, call "LoadSystem"
% to initialize openMatlab, which will also create a valid HW object.
%
% The object is singleton. That means, only one instance of this class should
% exist in one Matlab session. Call "PD.HW.GetInstance()" to get a reference
% to that instance in any scope (after executing "LoadSystem" in the base
% workspace at least once).
%
% See also:
% LOADSYSTEM
%
% ----------------------------------------------------------------------------
% (C) Copyright 2015-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ----------------------------------------------------------------------------

end


%#function LoadGradAmp_Standard
%#function LoadGradSystem_Standard
%#function LoadHFAmp_Standard
%#function LoadHWStruct_Standard
%#function LoadLift_Standard
%#function LoadMMRT_Standard
%#function LoadMagnet_Standard
%#function LoadPiezo_Standard
%#function LoadProbe_Standard
%#function LoadRXTXCoil_Standard
%#function LoadRXTX_Standard
%#function PD.DICOM
%#function PD.Grad
%#function PD.MMRT
%#function PD.RX
%#function PD.TX
%#function PD.Talker
%#function handleHidden

