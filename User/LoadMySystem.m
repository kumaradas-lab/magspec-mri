% LoadMySystem - default

% This file was created for the device with serial number 103.
% Check if the currently connected device matches that serial number.
checkDeviceSerial(HW, 103, mfilename('fullpath'));

HW.fLarmor = 24610000.000; HW.B0 = HW.fLarmor/(HW.Gamma.H1/2/pi);

% DC-600 external gradient amplifier
%LoadGradAmp_DC600_SN_16;
