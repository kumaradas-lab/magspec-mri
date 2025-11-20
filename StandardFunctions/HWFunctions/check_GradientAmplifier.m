function check_GradientAmplifier(HW, checkLevel)
%% Check status of gradient amplifier
%
%     check_GradientAmplifier(HW, checkLevel)
%
% If the system is configured to be using the external gradient amplifier
% DC-600 and it cannot be detected or it signals an error state, this function
% throws an error.
% Also, if a DC-600 is connected but the software is not configured to use
% it, an error is thrown.
%
% There are two pins on the "Ext. Grad." interface of the console that are used
% for status checks. High level on these pins means the status is ok.
% If HW.Grad.Status1 is 1 or checkLevel is positive, the level of the Status1
% pin is checked. If HW.Grad.Status1 is 2, the level of the Status2 pin is only
% checked if checkLevel is positive.
% The same applies for the value of HW.Grad.Status2 with respect to the Status2
% pin of the "Ext. Grad." interface of the console.
%
% This function is called with checkLevel=0 immediately before a measurement
% starts. It is repeatedly called with checkLevel=1 while polling measurement
% data (i.e., during a running measurement).
%
% The function does nothing when using a dummy MRI device.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function AbortMeasurement
%#function PD.Talker

