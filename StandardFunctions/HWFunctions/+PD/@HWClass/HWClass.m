classdef HWClass < handleHidden
%% HWClass - Class with standard interface to Pure Devices ResearchLab
%
%   HW = PD.HWClass.GetInstance(rootPath)
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
% exist in one Matlab session. Call "PD.HWClass.GetInstance()" to get a
% reference to that instance in any scope (after executing "LoadSystem" in the
% base workspace at least once).
%
% See also:
% LoadSystem
%
%
% PROPERTIES:
%
%   fLarmor
%       Larmor frequency in Hertz. This frequency is used by default in most
%       measurement sequences.
%
%   B0
%       Amplitude of the magnetic field density in Tesla. The magnetic field
%       density of a permanent magnet can change with its temperature. In many
%       cases, its amplitude is determined by measuring the Larmor frequency
%       with NMR.
%
%   GammaDef
%       Gyromagnetic ratio of the default nucleus in rad*Hz/Tesla. The
%       gyromagnetic ratio connects the Larmor frequency and the magnetic
%       field density: HW.fLarmor = HW.GammaDef / (2*pi) * HW.B0
%
%   fLarmorX
%       Larmor frequency of an X-nucleus in Hertz. This is used by default as
%       the second frequency in experiments which use dual-frequency rf-pulses
%       and/or dual-frequency acquisitions.
%
%   GammaX
%       Gyromagnetic ration of the X-nucleus in rad*Hz/Tesla. The
%       gyromagnetic ratio connects the Larmor frequency and the magnetic
%       field density for the X-nucleus:
%       HW.fLarmorX = HW.GammaX / (2*pi) * HW.B0
%
%   RootPath
%       String with the absolute root path of the OpenMatlab installation that
%       is currently used.
%
%   UserPath
%       String with the absolute or relative path to the User folder
%       containing calibration and settings files. In case, this is a relative
%       path, it is relative to HW.RoorPath. Among others, these files can
%       include "LoadMySystem.m" with general configuration settings,
%       "MagnetShimCal.m" with calibration values for the magnet shim, or
%       "PaUout2AmplitudeCal.m" with calibration value for the rf coil
%       transmission efficiency.
%
%   MagnetShimPath
%       String with the absolute path of the file containing calibration
%       settings for the magnet shim. By default, this is a file named
%       "MagnetShimCal.m" inside the folder given by HW.UserPath.
%
%   NetworkCalPath
%       String with the absolute path of the file containing calibration
%       settings for the rf coil network. By default, this is a .mat file
%       named "NetworkCal.mat" inside the folder given by HW.UserPath.
%
%   UserName
%       String with the currently used user name. User names can be used to
%       manage different configurations within the same installation of
%       OpenMatlab. Each user has its own subfolder inside the folder "User"
%       that contains files with calibration values and settings. If
%       HW.UserName is empty or 'default', the files directly inside the
%       folder "User" are used (not any of the potential subfolders).
%
%   UserNameDefault
%       String with the default user name. The user that will be loaded for
%       the connected console if "LoadSystem" is called without specifying a
%       user name.
%
%   UserNameList
%       Cell array of strings with user names that have pre-configured
%       settings.
%
%   UserNameListComment
%       Cell array of strings the same size as HW.UserNameList. These strings
%       will be displayed as comments to each user name in the list that might
%       be displayed on "LoadSystem".
%
%   UserNameMessage
%       Boolean value that can be used to turn off displaying the list of
%       pre-configured users on each call of "LoadSystem".
%
%   Function_Prepare_Measurement
%       Function handle with the following signature:
%         [HW, Seq, AQ, TX, Grad] = @(HW, Seq, AQ, TX, Grad)
%       If not empty, this function is called by "set_sequence" on its input
%       before the pulse program is created.
%
%   Function_Before_Measurement
%       Function handle with the following signature:
%         [HW] = @(HW)
%       If not empty, this function is called immediately before a measurement
%       is started.
%
%   Function_While_Measurement
%       Function handle with the following signature:
%         [HW] = @(HW)
%       If not empty, this function is called repeatedly while waiting for the
%       measurement to finish.
%
%   Function_After_Measurement
%       Function handle with the following signature:
%         [HW] = @(HW)
%       If not empty, this function is called immediately before receiving the
%       measurement data from the device.
%
%   Function_Abort_Measurement
%       Function handle with the following signature:
%         @(HW)
%       If not empty, this function is called when a measurement is aborted.
%
%   PostProcessingFcn
%       Function handle the wth following signature:
%         [Seq, data_1D, data_S_AQs_TRs] = @(nargout_set_sequence, ...
%           Seq, raw_data_AQ, num_averages, channel, iDevice, nucleusX)
%       This function is called by set_sequence to sort convert the raw data
%       received by the device into data structures that are more convenient
%       in Matlab. (Default: @Postprocessing)
%
%   Function_Run_Drawnow_While_Measurement
%       Boolean value to select whether the Matlab function "drawnow()" is
%       called (repeatedly) while waiting for a measurement to finish. The
%       function "drawnow()" allows (among other things) graphics to update.
%       However, calling it could also lead to delays. That might not be
%       desirable in time-critical measurements.
%
%   DefSeqValues
%       Structure with fields for holding default values for some pulse
%       programs.
%
%   FindFrequencyPlot
%       Boolean value to select whether the "Find_Frequency_*" functions plot
%       their measurement results.
%
%   FindFrequencyPause
%       Scalar value with the pause time in seconds after the frequency is
%       determined by the "Find_Frequency_Sweep" function. Adjust that value
%       to be appropriate for the (longest) T1 time of the used sample.
%
%   FindFrequencyPauseFID
%       Scalar value with the pause time in seconds after the frequency is
%       determined by the "Find_Frequency_FID" function. Adjust that value
%       to be appropriate for the (longest) T1 time of the used sample.
%
%   FindFrequencySweep
%       Structure with default settings for the "Find_Frequency_Sweep"
%       function.
%
%   FindFrequencyFID
%       Structure with default settings for the "Find_Frequency_FID" function.
%
%   FindShim
%       Structure with default settings for the "Find_Shim" function.
%
%   FindPulseDuration
%       Structure with default settings for the "Find_PulseDuration" function.
%
%   Network
%       Structure with default settings for "GUI_NetworkAnalyzer".
%
%   PlotSequence
%       Structure with default settings for the "plot_Seq" function.
%
%   RecoveryCPMG
%       Structure with default settings for the "sequence_RecoveryCPMG"
%       function.
%
%   Spectro
%       Structure with default settings for spectroscopy functions.
%
%   ReInit
%       Boolean value to indicate whether a short tRep should be added to the
%       front of each experiment to make sure all amplifiers have powered up
%       before the actual measurement starts.
%
%   tRepInit
%       Duration in seconds of the tRep that is prepended to each measurement
%       if HW.ReInit is true.
%
%   DICOM
%       Structure with DICOM file settings. (Currently unused)
%
%   Gamma
%       Structure with literature values of gyromagnetic ratios of different
%       nuclei.
%
%   FindFrequencyGamma
%       Gyromagnetic ration in rad*Hz/Tesla that is used in the
%       "Find_Frequency_*" functions. That gyromagnetic ration might differ
%       from HW.Gamma (e.g., to determine the magnet frequency at 1H frequency
%       while running the actual measurement for a X-nucleus).
%
%   Constant
%       Structure with literature values for some physical constants.
%
%   RelativCableSpeed
%       Relative factor (between 0 and 1) to the speed of light to estimate
%       the phase shift due to cable length. This is used in network analyzer
%       mode.
%
%   AddCableLength
%       Cable length in meter used in network analyzer mode to compensate
%       phase shift.
%
%   TempR50
%       Temperature in Kelvin for noise considerations.
%
%   MMRT
%       Object of type PD.MMRT with properties of the used mobile MRT. In case
%       multiple devices are connected, this is a vector of these objects
%       corresponding to each connected device.
%
%   TX
%       Object of type PD.TXClass with transmission properties of the
%       connected device. In case multiple devices are connected, this is a
%       vector of these objects corresponding to each connected device.
%
%   RX
%       Object of type PD.RX with receiption properties of the connected
%       device. In case multiple devices are connected, this is a vector of
%       these objects corresponding to each connected device.
%
%   Grad
%       Object of type PD.GradClass with gradient system and amplifier
%       properties of the connected device. In case multiple devices are
%       connected, this is a vector of these objects corresponding to each
%       connected device.
%
%   DigitalIO
%       Structure with setting of the digital input and output port. In case
%       multiple devices are connected, this is a vector of these structures
%       corresponding to each connected device.
%
%   MPI
%       Structure with properties for MPI (magnetic particle imaging)
%       measurement systems (if applicable).
%
%   Magnet
%       Optionally, this is an object of type PD.MagnetBase (or of a derived
%       class) with properties and interfaces to the control in a connected
%       magnet.
%
%   SampleHeater
%       Optionally, this is an object of type PD.SampleHeaterBase (or of a
%       derived class) with properties and interfaces to the control in a
%       connected sample heater device.
%
%   Piezo
%       Structure with properties of an optionally connected Piezo element.
%
%   Lift
%       Optionally, this is an object of type PD.SampleLift with properties
%       and interfaces to a connected sample lift.
%
%   PowerSupply
%       Optionally, this is an object of a supported type with properties and
%       interfaces to a connected power supply.
%
%   GUI
%       Enumeration object to select optional graphical user interfaces for
%       special applications.
%
%   DefaultDeviceOrder
%       If multiple devices are connected, this is a vector of serial numbers
%       of these devices for a stable order in the properties HW.MMRT, HW.RX,
%       HW.TX, and HW.Grad.
%
%   tFlip180Def
%       Read-only value with the duration of an rf pulse in seconds that flips
%       the spin system by 180° at the default rf amplitude (HW.TX.AmpDef) at
%       the default frequency (HW.fLarmor).
%
%   tFlip90Def
%       Read-only value with the duration of an rf pulse in seconds that flips
%       the spin system by 90° at the default rf amplitude (HW.TX.AmpDef) at
%       the default frequency (HW.fLarmor).
%
%   tFlip1Def
%       Read-only value with the duration of an rf pulse in seconds that flips
%       the spin system by 1° at the default rf amplitude (HW.TX.AmpDef) at
%       the default frequency (HW.fLarmor).
%
%   tFlip180min
%       Read-only value with the duration of an rf pulse in seconds that flips
%       the spin system by 180° at the default rf amplitude (HW.TX.AmpMax) at
%       the default frequency (HW.fLarmor).
%
%   tFlip90min
%       Read-only value with the duration of an rf pulse in seconds that flips
%       the spin system by 90° at the default rf amplitude (HW.TX.AmpMax) at
%       the default frequency (HW.fLarmor).
%
%   tFlip1min
%       Read-only value with the duration of an rf pulse in seconds that flips
%       the spin system by 1° at the default rf amplitude (HW.TX.AmpMax) at
%       the default frequency (HW.fLarmor).
%
%   MagnetShim
%       1x4 vector with values for the magnet shim in T/m. If multiple devices
%       are connected, the size of this vector increases to describe the
%       possible shim values for the gradient channels of all connected
%       devices.
%
%   B0Estimated
%       The estimated (or expected) value of the B0 magnetic field amplitude
%       in Tesla.
%
%   B0Target
%       If the connected magnet has optional B0 shift coils, this value holds
%       the target value of the B0 field in Tesla after the B0 shift is
%       applied.
%
%   MagnetT2i
%       Approximate apparent T2* of a bulk water sample with the connected
%       magnet in seconds.
%
%   TX2RXdeadTime
%       Dead time between rf pulse and acquisition in seconds. Used by
%       default, e.g., in imaging sequences. Use "get_DeadTimeTX2RX" instead.
%
%   RX2TXdeadTime
%       "Dead time" between the end of acquisitions and the next rf pulse in
%       seconds. Used by default, e.g., in imaging sequences. Use
%       "get_DeadTimeRX2TX" instead.
%
%   NetworkCalPathStandard
%       String with the absolute path of the file containing default
%       calibration settings for the rf coil network that should be working
%       for most devices.
%
%   NetworkCal
%       Structure with calibration values for the rf coil network as loaded
%       from the file given in HW.NetworkCalPath
%
%   Dummy
%       Vector the same length as the number of "connected" devices with the
%       serial numbers of the devices for which the settings are loaded if
%       OpenMatlab is using the respective device in dummy mode  (i.e.,
%       without actually connected device) or zero otherwise. If non-zero,
%       random numbers are returned instead of actual measurement data.
%       To initiate dummy mode, call "LoadSystem dummy". To select settings
%       for a specific device, set mySave.DummySerial to the serial number of
%       the device for which the settings should be loaded, or append the
%       serial number of the respective device to the command, e.g.,
%       "LoadSystem dummy123". Changing the value of this property directly
%       does not have any actual effect.
%
%   DummyFastForwardFactor
%       If HW.Dummy is non-zero, this factor divides the time that is used for
%       dummy measurements. This value is ignored if HW.Dummy is zero.
%       (Default: 1)
%
%
% METHODS:
%
%   HW_struct = HW.ToStruct()
%       Return structure corresponding to the current HW object. Consider
%       saving that structure to a file instead of the HW object itself to
%       make sure all fields correspond to the current values when loading it.
%       (No update methods will be run automatically for a structure.)
%
%
% ----------------------------------------------------------------------------
% (C) Copyright 2015-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ----------------------------------------------------------------------------

end


%#function LoadGradAmp_Standard
%#function LoadGradSystem_Standard
%#function LoadHFAmp_Standard
%#function LoadHWStruct_Standard
%#function LoadLift_Standard
%#function LoadMMRT_Standard
%#function LoadMPI_Standard
%#function LoadMagnet_Standard
%#function LoadPiezo_Standard
%#function LoadPowerSupply_Standard
%#function LoadProbe_Standard
%#function LoadRXTXCoil_Standard
%#function LoadRXTX_Standard
%#function LoadSampleHeater_Standard
%#function PD.DICOM
%#function PD.GradClass
%#function PD.IgnoreWarning
%#function PD.MMRT
%#function PD.MPIClass
%#function PD.RX
%#function PD.TXClass
%#function PD.Talker
%#function handleHidden

