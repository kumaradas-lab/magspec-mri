classdef MPIClass < handle
%% Class with properties for MPI device (& communication to MPI device)
%
% A user should not create an instance of this class manually. Instead, it
% will be created automatically as HW.MPI if HW.MPI.useMPI is set to `true`,
% e.g., in the LoadMySystem.m script file.
%
% This class holds pulse program or device specific properties regarding MPI
% measurements. Additionally, it can be used to, e.g., query the temperature
% of a sensor that is integrated or connected to the MPI devices.
%
%
% ----------------------------------------------------------------------------
% (C) Copyright 2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ----------------------------------------------------------------------------

end
