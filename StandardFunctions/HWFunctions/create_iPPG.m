function iPPG = create_iPPG(HW, iDevice)
%% Create a simple pulse program that sets the MRT to a defined state
%
%   iPPG = create_iPPG(HW, iDevice)
%
% The returned pulse program does not contain any actual pulses or acquisitions
% during the tRep. It can be used to initialize the device to a defined state or
% to set a defined state after a running measurement was aborted (potentially in
% a state that could be harmful if it were kept for a long time).
%
% ------------------------------------------------------------------------------
% (C) Copyright 2017-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function PD.DigitalIO
%#function PD.Gradient
%#function PD.HFAcquisition
%#function PD.HFPulses
%#function PD.MRISequence
%#function PD.SequenceCommands
%#function PD.TXMaxDef
%#function PD.Talker

