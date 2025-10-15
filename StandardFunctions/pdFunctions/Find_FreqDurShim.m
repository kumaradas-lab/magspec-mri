function [ HW ,mySave] = Find_FreqDurShim( HW, mySave, maxtime, doplot)
% automatically find the correct frequency, RF 90 degree pulse duration and
% shim values; values will be appended to "LoadMagnet_magnet22_NrXX.m
% aferwards

[HW,mySave]  = Find_Frequency_Sweep( HW, mySave,maxtime,doplot);
[HW,mySave]  = Find_PulseDuration( HW, mySave,maxtime,doplot);
[HW,mySave]  = Find_Frequency_Sweep( HW, mySave,maxtime,doplot);
[HW,mySave]  = Find_Shim( HW, mySave, maxtime,doplot);
[HW,mySave]  = Find_PulseDuration( HW, mySave,maxtime,doplot);
LoadSystem
[HW,mySave]  = Find_Frequency_Sweep( HW, mySave,maxtime);

    
%% ------------------------------------------------------------------------    
% (C) Copyright 2012 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%------------------------------------------------------------------------    