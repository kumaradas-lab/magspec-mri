function [SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave)
%% Highly customizable FLASH (or gradient echo) sequence for imaging in 1, 2, or 3 dimensions
%
%       [SeqLoop, mySave] = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave)
%
% This function executes a Fast Low-Angle Shot (FLASH) or gradient echo (GRE)
% imaging experiment and reconstructs the corresponding image. The imaging block
% is highly customizable (including EPI read-outs). Different correction
% measurements can be executed and applied automatically.
% By default, a FLASH measurement is executed. For a GRE measurement, set the
% excitation angle to 90 degrees (see Seq.FlipAngle).
%
%
% INPUT:
%
%   HW
%           HW structure or object (see LoadSystem).
%
%   Seq
%           structure containing the settings for the imaging block. All
%           settings are optional and default values are used if they are
%           omitted. Among others the following fields can be set:
%
%       tEcho
%               Echo time in seconds. For EPI measurements, this is the
%               effective echo time, i.e. the time between the center of the
%               excitation pulse and the center of the echo with zero k-space
%               encoding. See also Seq.AQSlice.nImages for further details.
%               If set to 0, run a ZTE experiment (experimental support).
%               (Default: 5e-3)
%
%       RepetitionTime
%               Repetition time (i.e. the time between excitation pulses) in
%               seconds.
%
%       FlipAngle
%               Flip angle of the excitation in degrees. Set to 90 for a
%               gradient echo measurement.
%
%       T1
%               If Seq.FlipAngle is omitted or empty, it is set to the Ernst
%               angle (alpha_E) corresponding to a sample with the given T1 time
%               in seconds:
%                   alpha_E = acos (exp(-RepetitionTime/T1))
%
%       SingletRep
%               Boolean value to indicate whether the pulse program for each
%               EPI segment should be realized as single tReps.
%               (default: Seq.AQSlice.EPIFactor == 1)
%
%       Loops
%               Scalar number of repeated measurements (default: 1). This can be
%               used to increase SNR. The measurements are run independently. So
%               a.o. the Larmor frequency can be adapted to account for
%               temperature drift (see Seq.Find_Frequency_interval).
%
%       LoopsBreak
%               Scalar with the duration of breaks between loops in seconds.
%
%       LoopsRepetitionTime
%               Scalar with the repetition time (i.e. time between excitation
%               pulses) of loops in seconds. Seq.LoopsBreak and
%               Seq.LoopRepetitionTime cannot be set at the same time.
%
%       LoopsBreakExactly
%               Boolean value whether subsequent loops start exactly on console
%               time. Otherwise, the computer clock is used for timing.
%               (default: 0, i.e. computer clock).
%
%       average
%               Scalar number of repeated measurements (default: 1). In contrast
%               to Seq.Loops, these averages are done at a lower level. On the
%               one hand, that means that the memory foodprint is lower. On the
%               other hand, all average measurements are exactly the same (i.e.
%               it is not possible to adapt to a changed Larmor frequency
%               between repetitions).
%
%       averageBreak
%               Scalar with the duration of bresks between successive average
%               measurements in seconds.
%
%       CorrectAQWindowPhase
%               Boolean value to indicate whether the phase of the acquisition
%               windows should be corrected (default: 1). This is helpful if the
%               actual Larmor frequency drifts slightly from the used excitation
%               and acquisition frequencies. See also function "get_data".
%
%       CorrectPhase
%               Boolean value. If true, a short acquisition window is placed
%               immediately after the excitation pulse. This signal is used to
%               correct the phase of the acquisition window at the encoded
%               gradient echo to compensate for the drift of the Larmor
%               frequency.
%
%       CorrectPhaseDuration
%               Duration of the acquisition window used for Seq.CorrectPhase in
%               seconds.
%
%       CorrectPhaseAQtOffset
%               Minimum offset from excitation pulse center to start of the
%               acquisition window for Seq.CorrectPhase in seconds. If
%               necessary, this offset is extented to acommodate for the
%               excitation pulse (+ dead time) or the slice rephase
%               gradient (+ eddy current time CorrectPhaseAQtEC). (Default: 0)
% 
%       CorrectPhaseAQtEC
%               eddy current time between the slice rephase and start of the
%               acquisition window for Seq.CorrectPhase in seconds.
%
%       CorrectPhaseNRead
%               Number of samples taken at the acquisition window used for
%               Seq.CorrectPhase.
%
%       CorrectPhaseSeparate
%               Boolean value. Only applies if Seq.CorrectPhase is also true. If
%               true, separate excitations are used to track the frequency
%               (instead of a short acquisition window immediately after the
%               excitation pulse).
%
%       CorrectPhaseBlockSize
%               If Seq.CorrectPhaseSeparate is true, an acquisition window for
%               frequency tracking is inserted after each CorrectPhaseBlockSize
%               excitations with image encoding.
%               (Default: 4 if Seq.CorrectPhaseSeparate, 0 otherwise)
%
%       CorrectSliceRephase
%               Boolean value. If true, a correction measurement is taken
%               immediately before the actual image measurement to determine an
%               offset for the slice rephase pulse (eddy currents).
%
%       CorrectReadRephase
%               Boolean value. If true, a correction measurement is taken
%               immediately before the actual image measurement to determine an
%               offset for the read encoding pulse (eddy currents).
%
%       Function_Prepare_Measurement
%               A function handle to a prepare function with the following
%               signature:
%                     [HW, Seq, AQ, TX, Grad] = @(HW, Seq, AQ, TX, Grad)
%               It is executed after the read, phase, and slice pulses ("logical
%               blocks") are created with get_ReadParameter, get_PhaseParameter,
%               and get_SliceParameter and before the corresponding pulses are
%               added to the final pulse program (AQ, TX, and Grad). It can be
%               used to modify the sequence before it is eventually executed.
%               Also see the programming notes below.
%               Alternatively, this can be a cell array of function handles with
%               the same signature as described before. In that case, the
%               prepare functions are executed in the order as they appear in
%               the cell array.
%               (Default: [])
%
%       plotSeqAQ
%               Plot the pulse program wrapped after each gradient echo (train)
%               and zoomed in to a reasonable region. The value selects the
%               gradient channels for which the pulse program is plotted. This
%               overrides the settings for Seq.plotSeq (see set_sequence).
%               (Default: [])
%       plotZTEEndPoints
%               Plot graphics illustrating the end points of the acquired rays
%               in k-space for ZTE measurements. (Default: false)
%
%       CorrectB0Read
%               structure with settings for B0 measurement and read correction.
%               This structure might contain the following fields (default
%               values apply if the fields are omitted or empty):
%           Use
%                   Boolean value (default: false). If it is true, the B0 map is
%                   used to compensate the offset in read direction due to B0
%                   deviations. This can also be used with a B0 map that is
%                   acquired in a separate measurement. In this case, the B0
%                   data has to be supplied in Seq.dataB0 (see description of
%                   output values below).
%           Get
%                   Boolean value to set if the B0 map should be acquired
%                   (default: Seq.CorrectB0Read.Use && isempty(Seq.dataB0)). If
%                   it is true, two additional measurements are performed with
%                   different echo times. The phase difference between the
%                   resulting images is used to determine the B0 deviation. The
%                   results are returned in SeqLoop.dataB0.
%           Plot
%                   Boolean value (default: false). If true, the results for
%                   both images acquired for the B0 map are plotted with the
%                   same settings like the final image(s).
%           tEchoIncr
%                   echo time increment between the two measurements in seconds
%                   (default: 1e-3). The first of the two measurements is
%                   performed with Seq.tEcho. For the second measurement, the
%                   echo time is increased by this value.
%           MaxFreqOffset
%                   maximum frequency offset (for selection of RoI) of the B0
%                   map in Hz (default: 2000).
%           MinRelAmp
%                   minimum amplitude (relative to maximum amplitude in image)
%                   that is part of the RoI (default: 0.2).
%           MaxRelAmpDiff
%                   maximum relative amplitude deviation (from 1.0) between the
%                   images from the two echo times (default: 0.25).
%           RoIExtension
%                   extend the region of interest in all directions by
%                   approximately the set value in meter. If this is a scalar,
%                   that value applies to all directions. Otherwise, it applies
%                   to the dimensions as sorted in data.ImageZ (see below).
%                   (Default: 0.15*image size)
%           WaitForSample
%                   Boolean value (default: false). If true, a message dialog is
%                   shown after acquisition of the B0 map and the measurement is
%                   paused until confirmation to allow for a sample exchange.
%           ZeroFillWindowSize
%                   relative size of the k-space filter used for the B0 map
%                   (default: 0.4).
%
%       AQSlice
%               scalar structure defining the orientation and other properties
%               of the acquired slice. Among others the following fields can be
%               used:
%
%           thickness
%                   Thickness of the slice in meter (default: 0.005). During the
%                   excitation pulses a gradient is applied such that the
%                   bandwidth of the pulse corresponds to the set thickness. For
%                   a definition of the bandwidth see the used pulse shape
%                   function (see: Seq.AQSlice(1).excitationPulse).
%
%           excitationPulse
%                   Function handle to the pulse shape function used for the
%                   exciation pulses (default: @Pulse_Rect if thickness > 1m,
%                   @Pulse_RaisedCos otherwise). Type "Pulse_" followed by
%                   hitting the tabulator key to see a list of available pulse
%                   shape functions (all in folder "pdFunctions"). If you want
%                   to add your own pulse shape function, see function
%                   Pulse_Rect_Composite90. Note that not all pulse shape
%                   functions are sensible exciation pulses.
%
%           sizeRead
%                   scalar with the size of the image in read direction in meter
%                   (default: HW.Grad.ImageVol(6)-HW.Grad.ImageVol(5)). For ZTE
%                   experiments, this setting isn't used. Use
%                   Seq.AQSlice.sizePhase instead.
%
%           sizePhase
%                   1x3 vector with the size of the image in phase1, phase2, and
%                   phase3 directions in meter (default according to
%                   HW.Grad.ImageVol). By default, phase1 is parallel to the
%                   slice direction, phase3 is parallel to the read direction
%                   and the three phase directions form a right handed
%                   coordinate system.
%                   For ZTE experiments (Seq.tEcho=0), these values are used to
%                   define the resulting image size. Values must be set from the
%                   start, e.g., for 2d ZTE experiments, set sizePhase(1) and
%                   sizePhase(2). Additionally, these values are used to create
%                   a regular grid on the surface of an n-dimensional cuboid. A
%                   k-line is acquired in the direction pointing to each of
%                   these points. At the same time, these values are used when
%                   reconstructing the image.
%
%           nRead
%                   scalar with the integer number of pixels in read direction.
%                   For ZTE experiments, this is the number of samples read in
%                   each k-line.
%
%           nPhase
%                   1x3 vector with integer numbers of pixels in phase1, phase2,
%                   and phase3 direction in meter. For a 2d image, set nRead and
%                   nPhase(2) to values > 1. For a 3d image, set nRead,
%                   nPhase(1), and nPhase(2) to values > 1. For a 3d image
%                   (CSI), set nPhase(1), nPhase(2), and nPhase(3) to values >
%                   1. For ZTE experiments (Seq.tEcho=0), these values are used
%                   to create a regular grid on the surface of an n-dimensional
%                   cuboid. A k-line is acquired in the direction pointing to
%                   each of these points. At the same time, these values are
%                   used when reconstructing the image.
%
%           Resolution
%                   1x3 vector with the resolution of a 3d image in voxels per
%                   meter (phase1 x phase2 x phase3/read).
%                   Or 1x2 vector with the resolution of a 2d image in pixels
%                   per meter (phase2 x phase3/read).
%                   Or scalar with the resolution of a 1d profile in pixels per
%                   meter (phase3/read).
%                   If both nRead and sizeRead or nPhase and sizePhase are set,
%                   they take precedence.
%
%           ReadOS
%                   scalar with the oversampling factor in read direction
%                   (default: [] for an automatic oversampling factor).
%                   Seq.AQSlice(1).nRead * Seq.AQSlice(1).ReadOS must be an
%                   integer.
%
%           PhaseOS
%                   1x3 vector with the oversampling factors in phase1, phase2,
%                   and phase3 direction (default: [1 1 1]).
%                   Seq.AQSlice(1).nPhase .* Seq.AQSlice(1).PhaseOS must all be
%                   integers.
%                   This can be used to effectively remove image artifact
%                   "entering" from the opposite side of the sample and to
%                   increase SNR.
%
%           SamplingFactor
%                   Integer sampling factor that is used to (down-)sample the
%                   received signal. Both "ReadOS" and "SamplingFactor"
%                   effectively increase the sampling frequency of the
%                   acquisition. However, samples acquired for "ReadOS" are used
%                   in image reconstruction, thus potentially reducing ghosting
%                   and damping effects. In contrast, the acquired signal is
%                   downsampled by "SamplingFactor" before image reconstruction,
%                   thus reducing the total amount of memory needed for k-space
%                   and image. (Default: 1 or however large necessary to reach a
%                   possible sampling rate.)
%
%           HzPerPixMin
%                   Minimum bandwidth per pixel in Hz. That means the duration
%                   of the acquisition windows is approximately
%                   1/Seq.AQSlice(1).HzPerPixMin
%                   If 0, the acquisition time is maximized.
%
%           AcquisitionTime
%                   Approximate acquisition time in s (default: SeqtEcho/1.5).
%                   If Seq.AQSlice(1).HzPerPixMin is also set, it takes
%                   precedence.
%
%           alfa, phi, theta
%                   Euler angles in radians defining the orientation of the
%                   slice in the magnet (default: all 0). By default, phase1 is
%                   oriented towards positive x, phase2 towards positive y, and
%                   phase3/read towards positive z. "theta" corresponds to a
%                   rotation around the z-axis. "phi" corresponds to a rotation
%                   around the y'-axis (the prior y-axis after the previous
%                   rotation). "alfa" corresponds to a rotation around the
%                   z''-axis (the prior z-axis after the previous rotations). By
%                   default "alfa" rotated inside the plane of a 2d image.
%
%           nImages
%                   scalar integer with number of images for each excitation
%                   pulse (Multi-Gradient-Echo). By default, the corresponding
%                   gradient echoes are spaced linearly. The time between
%                   equivalent k-line acquisitions is given by Seq.tEcho. It is
%                   also possible to specify the effective echo time of each
%                   image individually. In this case, Seq.tEcho must be a vector
%                   of length nImages with increasing values.
%                   (Default: numel(Seq.tEcho) )
%
%           EPIFactor
%                   Number of k-lines (per image) acquired in one EPI gradient
%                   echo segment. (Default: 1)
%
%           EPISegments
%                   Number of segments for EPI. If Seq.AQSlice.EPIFactor and
%                   Seq.AQSlice.EPISegments are defined simultaneously, their
%                   product must match the total number of k-lines to be
%                   acquired.
%
%           EPIEchoSpacing
%                   Time between subsequent EPI echoes in seconds (also known as
%                   interecho time).
%                   (Default: Seq.tEcho / Seq.AQSlice.EPIFactor)
%
%           EPIReadBipolar
%                   If true, the read encoders in an EPI segment are emitted in
%                   alternating polarity. Otherwise, all read out gradients are
%                   emitted with equal polarity. (Default: 0)
%
%           EPIGradient
%                   If false, the spins are rephased after each k-line
%                   acquisition. Otherwise, the spins are dephased incrementally
%                   in each EPI segment. See also EPIGradientPhase and
%                   EPIGradientRead. (Default: Seq.AQSlice.EPIFactor > 1)
%
%           EPIGradientPhase
%                   If false, rephase encoders are used after each read out in
%                   an EPI segment to rephase the phase encoding the spins.
%                   Otherwise, the spins are dephased incrementally in phase
%                   direction within each EPI segment.
%                   (Default: Seq.AQSlice.EPIGradient)
%
%           EPIGradientRead
%                   If false, rephase encoders are used after each read out in
%                   an EPI segment to rephase the read encoding of the spins.
%                   Otherwise, the read rephase encoder is moved to the read
%                   dephase encoder of the following readout. This effectively
%                   allows the two encoders to cancel each other out if
%                   Seq.AQSlice.EPIReadBipolar is true.
%                   (Default: Seq.AQSlice.EPIGradient & Seq.AQSlice.EPIReadBipolar)
%
%           oddEvenEchoes
%                   Boolean value to acquire separate images consisting of only
%                   odd or even echoes, respectively. For this, the echoes
%                   are split into two echoes (effectively increasing the
%                   bandwidth and doubling the number of echoes acquired per
%                   EPI segment!). This is helpful to reduce ghost images when
%                   acquiring segmented Multi-Gradient-Echo images.
%
%           kLineOrderType
%                   String that selects in which order the k-lines are acquired.
%                   Possible (case insensitive) values are:
%                   * 'increasing': Acquire the k-space from the lowest k-space
%                     line to the highest (default).
%                   * 'decreasing': Acquire the k-space from the highest k-space
%                     line to the lowest.
%                   * 'centerOut': Start acquisition at the center of the
%                     k-space and proceed alternatingly to lower and higher
%                     k-space lines.
%                   * 'random': Acquire the k-space lines in a random order.
%
%           kLineOrder
%                   Permutation or selection vector for the k-space lines (1
%                   being the lowest k-space line index). If this field exists
%                   and is not empty, it takes precedence over
%                   Seq.AQSlice.kLineOrderType (Default: []).
%
%           SliceRephaseImmediately
%                   Logical value to indicate that the slice gradient should be
%                   rephased immediately after the slice selection gradient. If
%                   false and Seq.CorrectPhase is false, the slice rephase
%                   gradient is simultaneous to the read out dephase gradient.
%                   If Seq.CorrectPhase is true, the slice gradient is always
%                   rephased immediately after the slice selection gradient
%                   independent from this setting. (Default: Seq.CorrectPhase)
%
%           DephaseLengthFactor
%                   scalar value with a factor for the duration of the dephase
%                   pulses (default: 1). A factor of 1 means that the dephase
%                   pulse duration is approximately have the acquisition time.
%
%           DephaseTimeOffset
%                   scalar value for the time offset of the dephase pulses in
%                   seconds (default: 0). A value of 0 means that the dephase
%                   pulses are immediately before the read out gradient pulse.
%
%           SpoilFactor
%                   1x3 vector with spoil factors in slice/phase1, phase2, and
%                   read/phase3 direction. A factor of 1 means that the spoilers
%                   dephase the spins over the thickness/pixels/voxels by 360
%                   degrees.
%
%           sizePhaseSpoil
%                   1x3 vector with spoil sizes in slice/phase1, phase2, and
%                   read/phase3 direction in mm. Spins dephase by 180 degrees
%                   over this distance. If Seq.AQSlice.SpoilFactor is set at the
%                   same time, this value takes precedence.
%
%           SpoilLengthFactor
%                   scalar value with a factor for the duration of the spoiler
%                   pulses (default: 1). A factor of 1 means that the spoiler
%                   has the same duration as the phase gradient pulses.
%
%           SpoilTimeOffset
%                   scalar value for the time offset of the spoiler pulses in
%                   seconds (default: 0). A value of 0 means that the spoiler
%                   pulses coincide with the phase gradient pulses.
%
%
%   AQ
%           structure with settings for the acquisition windows. It is not
%           necessary to set anything here. The default AQ structure that is
%           created when calling "LoadSystem" can be passed. Any acquisition
%           windows defined in this structure are added to the ones from the
%           main imaging block.
%
%   TX
%           1x2 structure with settings for the rf pulses. It is not necessary
%           to set anything here. The default TX structure that is created when
%           calling "LoadSystem" can be passed. Any rf pulses defined in this
%           structure are added to the ones from the main imaging block.
%
%   Grad
%           1x4 structure with settings for the gradients. It is not necessary
%           to set anything here. The default Grad structure that is created
%           when calling "LoadSystem" can be passed. Any gradient pulses defined
%           in this structure are added to the ones from the main imaging block.
%
%   mySave
%           [optional] mySave structure as generated when calling "LoadSystem".
%           It is used for tracking the current Larmor frequency and timing of
%           Larmor frequency sweeps between different functions (or function
%           calls).
%
%
% OUTPUT:
%
%   SeqLoop
%           Same as the input structure Seq with added fields containing the
%           actually used field values, the pulse program (AQ, TX, Grad,
%           DigitalIO) and the results (data).
%           Amongst others the added fields are:
%
%     TX
%             Structure containing matrices describing the rf pulses in the
%             pulse program. See "manual openMatlab" for further details.
%
%     AQ
%             Structure containing matrices describing the acquisition windows
%             in the pulse program. See "manual openMatlab" for further details.
%
%     Grad
%             1x4 structure containing matrices describing the gradient pulses
%             in the pulse program. See "manual openMatlab" for further details.
%
%     DigitalIO
%             Structure containing matrices describing digital output signals in
%             the pulse program. See "manual openMatlab" for further details.
%
%     HW
%             Structure holding all content of HW as used for the measurement.
%
%     data
%             Structure containing the measured data. Amongst others these
%             fields are:
%
%       time_of_tRep, time_all, data, fft1_data, f_fft1_data:
%               See function "get_data".
%
%       kSpaceOsRaw, ImageOsRaw, kSpaceOsZ, ImageOsZ, kSpaceZ, ImageZ, kSpace,
%       Image:
%               See function "get_kSpaceAndImage".
%
%       kTicks, Ticks:
%               See function "get_kSpaceAndImageTicks".
%
%       tImageZ
%               Time (in s) measured between center of excitation pulse and echo
%               center for each image in the echo train.
%
%     dataB0
%             In case Seq.CorrectB0Read.Get is set to true, this structure
%             contains data with the measured B0 map.
%
%   mySave
%           mySave structure containing (most notably) the current Larmor
%           frequency and timing information for Larmor frequency sweeps.
%
%
% PROGRAMMING NOTES:
%
%   The pulse program is constructed using logical blocks that correspond to
%   slice, read, and phase encoders.
%   The following functions translate the logical blocks to the actual pulse
%   program:
%     For slice units:      get_SliceParameter
%     For read units:       get_ReadParameter
%     For phase units:      get_PhaseParameter
%   Please, see the documentation of these functions for further details.
%
%   The pulse program consists of multiple tReps. It is implemented by using
%   separate tReps for each (slice-selective) excitation pulse and for each
%   read-out.
%
%   When manipulating the pulse program with Seq.Function_Prepare_Measurement,
%   the settings for these units can be changed. Currently, the following units
%   are used:
%     Seq.Read(1)
%       Read out at tEcho
%     Seq.Read(2)
%       Dummy read out at tEcho of pre and post shots
%     Seq.Read(3)
%       Dummy read out to compensate gradients at last tRep of each EPI segment
%       if Seq.AQSlice(1).EPIGradientRead
%     Seq.Read(4)
%       Dummy read out to compensate gradients at last tRep of each EPI segment
%        of pre and post shots if Seq.AQSlice(1).EPIGradientRead
%     Seq.Read(5)
%       (Un-encoded) read out immediately after pulse to track the frequency if
%       Seq.CorrectPhase > 0
%     Seq.Slice(1)
%       excitation pulse (with optional slice selection)
%     Seq.Phase(1)
%       Phase encoding in phase(1) direction. By default, the slice direction is
%       parallel to phase(1).
%     Seq.Phase(2)
%       Phase encoding in phase(2) direction. By default, the read and slice
%       directions are perpendicular to phase(2).
%     Seq.Phase(3)
%       Phase encoding in phase(3) direction. By default, the read direction is
%       parallel to phase(3).
%     Seq.Phase(4)
%       Spoiler/crusher in phase(1)/slice direction aligned with read rephase at
%       last echo.
%     Seq.Phase(5)
%       Spoiler/crusher in phase(2) direction aligned with read rephase at last
%       echo.
%     Seq.Phase(6)
%       Spoiler/crusher in phase(3)/read direction aligned with read rephase at
%       last echo.
%     Seq.Phase(7)
%       Dummy spoiler that moves the rephase gradient in phase(1)/slice
%       direction to the dephase gradient of the next EPI echo within the EPI
%       echo train if Seq.AQSlice(1).EPIGradientPhase
%     Seq.Phase(8)
%       Dummy spoiler that moves the rephase gradient in phase(2) direction to
%       the dephase gradient of the next EPI echo within the EPI echo train if
%       Seq.AQSlice(1).EPIGradientPhase
%     Seq.Phase(9)
%       Dummy spoiler that moves the rephase gradient in phase(3)/read
%       direction to the dephase gradient of the next EPI echo within the EPI
%       echo train if Seq.AQSlice(1).EPIGradientPhase
%
%
%   See the documentation of "set_sequence" for information about exact timing
%   between subsequent measurements.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2013-2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
%
% See also: sequence_Spin_Echo, get_ReadParameter, get_PhaseParameter,
%           get_SliceParameter, set_sequence, get_data, get_kSpaceAndImage,
%           get_kSpaceAndImageTicks, plot_kSpaceAndImage
%
% Note:

%% Flash V2.0           11.04.2014
%                      Toni, Driessle - Pure Devices GmbH
% Flash 1 2 and 3 dimension
%
% LoadSystem                                              % Load system parameters (Reset to default: HW Seq AQ TX)
% Seq.Loops=1;                                            % Number of loop averages    1...
% Seq.LoopsBreak=[];                                      % Pause between two loop averages in seconds ([]= fast as possible)
% Seq.LoopSaveAllData=0;                                  % Save all data of the Loops
% Seq.LoopPlot = 1;                                       % Plot image after each loop
% Seq.Find_Frequency_interval = 100;                      % Time between two frequency search
%
% Seq.plotSeqTR=1:3;                                      % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
% Seq.plotSeq=1:3;                                        % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
% Seq.plotSeqStart=8;                                     % Plot sequence start with tRep x
% Seq.plotSeqEnd=12;                                      % Plot sequence stop with tRep x
%
% Seq.T1=100e-3;                                          % T1 of probe excitation is acos(exp(-Seq.tRep/Seq.T1))/pi*180
% Seq.tEcho=5e-3;                                         % Echo time in seconds eg. 5e-3
% Seq.RepetitionTime=10e-3;                               % Repetition time in seconds (default is Seq.tEcho*2)
%
%
% % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.AQSlice(1).nRead=16;                                % Number of Pixels in read, if nRead>1 nPhase(3)=1
% Seq.AQSlice(1).nPhase(1)=16;                            % Number of Pixels in phase(1)
% Seq.AQSlice(1).nPhase(2)=16;                            % Number of Pixels in phase(2)
% Seq.AQSlice(1).nPhase(3)=1;                             % Number of Pixels in CSI, if nPhase(3)>1 nRead=1 sizeRead=1e12
% Seq.AQSlice(1).HzPerPixMin=500;                         % Bandwidth per pixel in Hz (1/HzPerPixMin= duration of AQ) eg. 0 vor max tAQ
% Seq.AQSlice(1).sizeRead=0.016;                          % Image size in read in meter (for CSI set to 1e12)
% Seq.AQSlice(1).sizePhase(1)=0.016;                      % Image size in phase in meter
% Seq.AQSlice(1).sizePhase(2)=0.016;                      % Image size in phase(2) in meter
% Seq.AQSlice(1).sizePhase(3)=0.016;                      % Image size in phase(3) in meter
% Seq.AQSlice(1).thickness=0.005;                         % Image thickness in slice direction  used for 2D and 3D! ([] for no Slice) in meter
% Seq.AQSlice(1).excitationPulse=@Pulse_Rect;             % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
%
% % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.AQSlice(1).ReadOS=16;                               % Oversampling read ([] for automatic, recommended >=2) 1...
% Seq.AQSlice(1).PhaseOS(1)=1;                            % Oversampling phase(1)  1...
% Seq.AQSlice(1).PhaseOS(2)=2;                            % Oversampling phase(2)  1...
% Seq.AQSlice(1).PhaseOS(3)=1;                            % Oversampling phase(3)  1...% % rotate all image coordinates around... rotate
%
% % Plot        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.LoopPlot=0;                                         % Plot every Loop
% Seq.AQSlice(1).plotkSpace=0;                            % Plot k-Space
% Seq.AQSlice(1).plotImage=1;                             % Plot Image
% Seq.AQSlice(1).plotPhase=0;                             % Plot phase of k-Space or Image
% Seq.AQSlice(1).plotB0ppm=0;                             % Plot B0 ppm
% Seq.AQSlice(1).plotB0Hz=0;                              % Plot B0 Hz
% Seq.AQSlice(1).ZeroFillWindowSize=1.4;                  % Zero Fill Window Size (k-Space)
% Seq.AQSlice(1).ZeroFillFactor=2;                        % Zero fill resolution factor
% Seq.AQSlice(1).plotImageHandle=111+Loop;                % figure handle of the 2D plot
%
% % Orientation in Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.AQSlice(1).alfa=0.0*pi;                             % first rotation around x axis in RAD
% Seq.AQSlice(1).phi=0.0*pi;                              % second rotation around y axis in RAD
% Seq.AQSlice(1).theta=0.0*pi;                            % third rotation around z axis in RAD
% Seq.AQSlice(1).Center2OriginImage = [0,0,0];            % coordinate of the origin (of the gradient system) in the image coordinate system (slice or phase(1) x phase(2) x read or phase(3)) in m
% Seq.AQSlice(1).ReadCoordinate=1;                        % direction of Read:  x = 1, y = 2, z = 3
% Seq.AQSlice(1).PhaseCoordinate(1)=1;                    % direction of Phase(1):  x = 1, y = 2, z = 3
% Seq.AQSlice(1).PhaseCoordinate(2)=2;                    % direction of Phase(2):  x = 1, y = 2, z = 3
% Seq.AQSlice(1).PhaseCoordinate(3)=3;                    % direction of Phase(3):  x = 1, y = 2, z = 3
% Seq.AQSlice(1).SliceCoordinate=3;                       % direction of Slice:  x = 1, y = 2, z = 3
%
% % Some corrections    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.SteadyState_PreShots=10;                            % number of tReps before AQ starts(used for CorrectPhase filter and steady state)
% Seq.SteadyState_PostShots=10;                           % number of tReps after AQ stops (used for CorrectPhase filter)
% Seq.CorrectPhase=1;                                     % correct linear phase and offset of each tRep
% Seq.CorrectPlotFrequency=0;                             % plot frequency vs time
% HW.Grad.SliceTimeDelay=[0,0,0]*1e-6;                    % delay only for the Slice gradient
% HW.Grad.ReadTimeDelay=[0,0,0]*1e-6;                     % delay only for the read gradient
% HW.Grad.PhaseTimeDelay=[0,0,0]*1e-6;                    % delay only for the phase gradient
% Seq.average=1;                                          % Number of averages 1... (Seq.average >= 2 is not working properly with Seq.CorrectPhase=1, please use Seq.Loops)
% Seq.averageBreak=1;                                     % Pause between two averages in seconds
% Seq.LoopsBreakExactly=0;                                % Lock Time between two loops
% Seq.AQSlice(1).ReadOSUsedForImage=17;                   % Number of samples used in CSI
% Seq.SweepPhase=1;                                       % cumsum phase increment each tRep 1 on 0 off
% Seq.PhaseIncrement=0;                                   % phase increment each tRep in deg
% Seq.MaxGradAmpSlice=0.1;                                % limit slice gradient strength
% Seq.FlipAngle=90;                                       % set excitation flip angle
% Seq.CorrectAQWindowPhase=1                              % Correct phase of AQ Window if the AQ.Frequency is not HW.fLarmor
% Seq.AQSlice(1).ReadGradRephaseSign=1;                   % 1 dephase, -1 rephase on read direction
% Seq.AQSlice(1).sizePhaseSpoil=0.001;                    % Spoiler size in Seq.AQSlice(1).PhaseCoordinate(1) direction
% Seq.AQSlice(1).SpoilLengthFactor=1;                     % Double Spoil Length to decrease Amp
% Seq.AQSlice(1).RephasePhaseEncoding=0;                  % Rephase Phase Encoding
% Seq.AQSlice(1).SliceDephase=0;                          % Dephase Slice Encoding
% Seq.PhaseIncrementLoop=0;                               % Phase Increment Loop to Loop
%
% % Seq.CorrectSliceRephase=0;                            % Correct SliceGradTimeIntegralRephaseOffset
% % Seq.CorrectReadRephase=0;                             % Correct ReadGradTimeIntegralOffset
% % Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset=SeqLoop.AQSlice(1).SliceGradTimeIntegralRephaseOffset;
% % Seq.AQSlice(1).ReadGradTimeIntegralOffset=SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset;
% Seq.CorrectSliceRephasePlot=0;
% Seq.CorrectReadRephasePlot=0;
%
% SeqLoop = sequence_Flash(HW, Seq, AQ, TX, Grad, mySave);


%% default parameters
if nargin < 3
  AQ.Start = [];        % matrix containing the starting times (center of the first Sample – 0.5/AQ.fSample) of each acquisition window for each TR; up to 510 AQ windows are possible per TR
  AQ.nSamples = [];     % matrix containing the number of samples to be acquired for each AQ window for each TR
  AQ.fSample = [];      % matrix containing the sampling frequency of each acquisition window for each TR
  AQ.Frequency = [];    % matrix containing the mixing frequency of each acquisition windows for each TR
  AQ.Phase = [];        % matrix containing the phase offset for each acquisition window
  AQ.ResetPhases = [];  % if set to 1 the TX and RX phase is matched at the beginning of the TR
  AQ.Gain = [];         % relative gain of the acquisition path (1/AQ.Gain = max input Amplitude). Use HW.RX.GainDef for the best noise figure.
  AQ.Repeat = [];       % if there is a ‘1’ in the array, the settings for the window of the prior TR are used. (reduces data traffic between the PC and the MRI device)
end

if nargin < 4
  TX.Channel = HW.TX.ChannelDef;  % use transmitting channel 1 or 2
  TX.Start = [];        % matrix containing the starting time of the RF pulses for each TR column by column; up to 510 RF pulses are possible per TR (column)
  TX.Frequency = [];    % matrix containing transmitting frequency for each pulse of each TR
  TX.Duration = [];     % matrix containing the duration of the RF pulses for each TR column by column; up to 510 RF pulses are possible per TR
  TX.Amplitude = [];    % matrix containing the amplitude for each pulse of each TR
  TX.Phase = [];        % matrix containing the phase offset for each pulse of each TR

  TX.BlankOffset = [];  % array of times to define, when the blanking signal is applied prior to the RF pulse; can be adjusted every TR
  TX.BlankPostset = []; % array of times to define how long blanking stays active after a RF pulse; can be adjusted every TR
  TX.Repeat = [];       % if there is a ‘1’ in the array, the pulses of the prior TR are used. (reduces data traffic between the PC and the MRI device)
end

if nargin < 5
  [Grad(1:HW.Grad.n).Time] = deal([]);   % time the values in Grad(t).Amp with the corresponding index are set
  [Grad(1:HW.Grad.n).Amp] = deal([]);    % amplitude of the gradient at the time Grad(t).Time (linearly interpolated)
  [Grad(1:HW.Grad.n).Shim] = deal([]);   % additional shim; magnet shim is already considered in HW.MagnetShim. Caution: Using high values over a long time will damage the gradient coils and amplifiers!
  [Grad(1:HW.Grad.n).Repeat] = deal([]); % if there is a one in the array the gradients of the prior TR are used. (reduces data traffic between the PC and the MRI device)
end

if nargin < 6, mySave = struct(); end

Seq = set_EmptyField(Seq, 'PreProcessSequence', 1);
Seq = set_EmptyField(Seq, 'StartSequence', 1);
Seq = set_EmptyField(Seq, 'PollPPGfast', 1);
Seq = set_EmptyField(Seq, 'GetRawData', 1);
Seq = set_EmptyField(Seq, 'PostProcessSequence', 1);

% add a few empty fields to make them appear early on in the Seq structure
if ~isfield(Seq, 'tEcho'), Seq.tEcho = []; end
if ~isfield(Seq, 'RepetitionTime'), Seq.RepetitionTime = []; end
if ~isfield(Seq, 'T1'), Seq.T1 = []; end
if ~isfield(Seq, 'FlipAngle'), Seq.FlipAngle = []; end

if ~isfield(Seq, 'AQSlice'), Seq.AQSlice = struct(); end
if ~isfield(Seq.AQSlice, 'iDevice'), Seq.AQSlice(1).iDevice = 1; end  % make sure that AQSlice is not empty

Seq = set_EmptyField(Seq, 'Loops', 1);
if ~isfield(Seq, 'LoopsBreak'), Seq.LoopsBreak = []; end
if ~isfield(Seq, 'LoopsRepetitionTime'), Seq.LoopsRepetitionTime = []; end
Seq = set_EmptyField(Seq, 'LoopsBreakExactly', 0);
if ~isfield(Seq, 'StartSequenceTime'), Seq.StartSequenceTime = []; end
Seq = set_EmptyField(Seq, 'CorrectReadRephasePlot', 0);
Seq = set_EmptyField(Seq, 'CorrectSliceRephasePlot', 0);
Seq = set_EmptyField(Seq, 'CLTime', 10e-6);
Seq = set_EmptyField(Seq, 'LoopPlot', 1);
Seq = set_EmptyField(Seq, 'LoopPlotAverages', 1);
Seq = set_EmptyField(Seq, 'LoopSeqPlot', 0);
Seq = set_EmptyField(Seq, 'CorrectAQWindowPhase', 1);
Seq = set_EmptyField(Seq, 'average', 1);
if isemptyfield(Seq, 'Reinitialize'), Seq.Reinitialize = 1; end
if Seq.tEcho == 0
  isZTE = true;
else
  isZTE = false;
end
if isemptyfield(Seq, 'CorrectPhase')
  % use (un-encoded) acquisitions to track magnet frequency
  Seq.CorrectPhase = 1;
end
if isemptyfield(Seq, 'CorrectPhaseSeparate')
  % use separate excitation for frequency tracking windows
  Seq.CorrectPhaseSeparate = isZTE;
end
if Seq.CorrectPhase && isZTE && ~Seq.CorrectPhaseSeparate
  error('PD:sequence_Flash:ZTENoCorrectPhase', ...
    'ZTE experiments require the frequency tracking windows to be at separate excitations.');
end
if isemptyfield(Seq, 'CorrectPhaseBlockSize')
  % number of encoded echoes before acquiring a non-encoded echo
  Seq.CorrectPhaseBlockSize = 4;
end
if ~Seq.CorrectPhaseSeparate || ~Seq.CorrectPhase
  Seq.CorrectPhaseBlockSize = 0;
end
Seq = set_EmptyField(Seq, 'CorrectPhaseDuration', 0.5e-3);
Seq = set_EmptyField(Seq, 'CorrectPhaseAQtOffset', 0);
Seq = set_EmptyField(Seq, 'CorrectPhaseAQtEC', HW.Grad(Seq.AQSlice(1).iDevice).tEC);
Seq = set_EmptyField(Seq, 'CorrectPhaseNRead', 16);
Seq = set_EmptyField(Seq, 'CorrectSliceRephase', 0);
Seq = set_EmptyField(Seq, 'CorrectReadRephase', 0);
Seq = set_EmptyField(Seq, 'CorrectPlotFrequency', 0);
Seq = set_EmptyField(Seq, 'SweepPhaseQuadratically', 0);
Seq = set_EmptyField(Seq, 'PhaseIncrement', 0);  % e.g.; HW.Constant.GoldenAngle138
Seq = set_EmptyField(Seq, 'PhaseIncrementLoop', 0);

if ~isfield(Seq, 'Function_Prepare_Measurement')
  Seq.Function_Prepare_Measurement = [];
end

Seq = set_EmptyField(Seq, 'MaxGradAmpSlice', HW.Grad(Seq.AQSlice(1).iDevice).MaxAmpSlice);
if isemptyfield(Seq.AQSlice(1), 'AmplitudeUnit'), Seq.AQSlice(1).AmplitudeUnit = HW.RX(Seq.AQSlice(1).iDevice).AmplitudeUnit; end
if isemptyfield(Seq.AQSlice(1), 'AmplitudeUnitScale'), Seq.AQSlice(1).AmplitudeUnitScale = HW.RX(Seq.AQSlice(1).iDevice).AmplitudeUnitScale; end
if isemptyfield(Seq.AQSlice(1), 'LengthUnit'), Seq.AQSlice(1).LengthUnit = HW.Grad(Seq.AQSlice(1).iDevice).LengthUnit; end
if isemptyfield(Seq.AQSlice(1), 'LengthUnitScale'), Seq.AQSlice(1).LengthUnitScale = HW.Grad(Seq.AQSlice(1).iDevice).LengthUnitScale; end
if isemptyfield(Seq.AQSlice(1), 'ReadCoordinate'), Seq.AQSlice(1).ReadCoordinate = 3; end
if isemptyfield(Seq.AQSlice(1), 'PhaseCoordinate'), Seq.AQSlice(1).PhaseCoordinate = [1 2 3]; end
if isemptyfield(Seq.AQSlice(1), 'SliceCoordinate'), Seq.AQSlice(1).SliceCoordinate = 1; end
if isemptyfield(Seq.AQSlice(1), 'plotImageHandle'), Seq.AQSlice(1).plotImageHandle = 110; end
if isemptyfield(Seq.AQSlice(1), 'UseAQWindow'), Seq.AQSlice(1).UseAQWindow = 1; end
if isemptyfield(Seq.AQSlice(1), 'SliceGradTimeIntegralRephaseOffset'), Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset = 0; end
if isemptyfield(Seq.AQSlice(1), 'ReadGradTimeIntegralOffset'), Seq.AQSlice(1).ReadGradTimeIntegralOffset = 0; end
if isemptyfield(Seq.AQSlice(1), 'ReadGradTimeDelayOffset'), Seq.AQSlice(1).ReadGradTimeDelayOffset = 0; end
if isemptyfield(Seq.AQSlice(1), 'SliceRephaseTimeOffset')
  % move slice rephase gradient (only used if Seq.CorrectPhase is true)
  Seq.AQSlice(1).SliceRephaseTimeOffset = 0;
end
if ~isfield(Seq.AQSlice(1), 'thickness') || isnan(Seq.AQSlice(1).thickness), Seq.AQSlice(1).thickness = []; end
if isemptyfield(Seq.AQSlice(1), 'SliceRephaseImmediately')
  % slice rephase gradient immediately after slice gradient.
  % Otherwise, simultaneously with read dephase.
  Seq.AQSlice(1).SliceRephaseImmediately = Seq.CorrectPhase;
end
if Seq.CorrectPhase && ~Seq.CorrectPhaseSeparate
  % FIXME: In principle, rephasing the slice immediately is only needed for the
  % tReps with tracking windows for Seq.CorrectPhaseSeparate.
  Seq.AQSlice(1).SliceRephaseImmediately = 1;
end

% structure with data for B0 correction
if ~isfield(Seq, 'dataB0'), Seq.dataB0 = []; end
% settings for B0 measurement and correction
if ~isfield(Seq, 'CorrectB0Read'), Seq.CorrectB0Read = struct(); end
% correct read offset with B0 map
if isemptyfield(Seq.CorrectB0Read, 'Use'), Seq.CorrectB0Read.Use = false; end
% get data for B0 correction
if isemptyfield(Seq.CorrectB0Read, 'Get'), Seq.CorrectB0Read.Get = Seq.CorrectB0Read.Use && isempty(Seq.dataB0); end
% plot results for both images acquired for B0 map
if isemptyfield(Seq.CorrectB0Read, 'Plot'), Seq.CorrectB0Read.Plot = false; end
% echo time increment for second measurement in s
if isemptyfield(Seq.CorrectB0Read, 'tEchoIncr'), Seq.CorrectB0Read.tEchoIncr = 1e-3; end
% maximum frequency offset of B0 map in Hz
if isemptyfield(Seq.CorrectB0Read, 'MaxFreqOffset'), Seq.CorrectB0Read.MaxFreqOffset = 2000; end
% minimum amplitude relative to maximum amplitude
if isemptyfield(Seq.CorrectB0Read, 'MinRelAmp'), Seq.CorrectB0Read.MinRelAmp = 0.2; end
% maximum relative amplitude deviation from 1.0 between measurements 1 and 2
if isemptyfield(Seq.CorrectB0Read, 'MaxRelAmpDiff'), Seq.CorrectB0Read.MaxRelAmpDiff = 0.25; end
% after measuring B0 map, wait for sample exchange
if isemptyfield(Seq.CorrectB0Read, 'WaitForSample'), Seq.CorrectB0Read.WaitForSample = false; end
% ZeroFillWindowSize used for B0map
if isemptyfield(Seq.CorrectB0Read, 'ZeroFillWindowSize'), Seq.CorrectB0Read.ZeroFillWindowSize = 0.4; end
% if Seq.CorrectB0Read.Get, Seq.LoopSaveAllData = 1; end  % Force LoopSaveAllData to on
if Seq.CorrectB0Read.Use && ~Seq.CorrectB0Read.Get && isempty(Seq.dataB0)
  warning('PD:sequence_Flash:CorrectB0ReadImpossible', ...
    ['To be able to correct the read error caused by B0 deviations either the ', ...
    'B0 map must be measured (Seq.CorrectB0Read.Get) or a previously measured ', ...
    'B0 map must be provided (Seq.dataB0).\nB0 read correction will be skipped.']);
  Seq.CorrectB0Read.Use = false;
end

% settings for plot_kSpaceAndImage
if ~isfield(Seq.AQSlice(1), 'iSlice'),              Seq.AQSlice(1).iSlice         = [];   end
if isemptyfield(Seq.AQSlice(1), 'plotPhase'),       Seq.AQSlice(1).plotPhase      = 1;    end
if isemptyfield(Seq.AQSlice(1), 'plotkSpace'),      Seq.AQSlice(1).plotkSpace     = 1;    end
if ~isfield(Seq.AQSlice(1), 'plotkSpacePhase'),     Seq.AQSlice(1).plotkSpacePhase = [];  end
if isemptyfield(Seq.AQSlice(1), 'plotImage'),       Seq.AQSlice(1).plotImage      = 1;    end
if ~isfield(Seq.AQSlice(1), 'plotImagePhase'),      Seq.AQSlice(1).plotImagePhase = [];   end
if ~isfield(Seq.AQSlice(1), 'plotImageOs'),         Seq.AQSlice(1).plotImageOs    = [];   end
if ~isfield(Seq.AQSlice(1), 'plotImageOsPhase'),    Seq.AQSlice(1).plotImageOsPhase = []; end
if isemptyfield(Seq.AQSlice(1), 'plotImagehAxes'),  Seq.AQSlice(1).plotImagehAxes = cell(1, 4); end
if isemptyfield(Seq.AQSlice(1), 'plotImageOshAxes'),Seq.AQSlice(1).plotImageOshAxes = cell(1, 4); end
  if isemptyfield(Seq.AQSlice(1), 'plotDatahAxes'), Seq.AQSlice(1).plotDatahAxes = cell(1, 4); end
if ~isfield(Seq.AQSlice(1), 'plotB0ppm'),           Seq.AQSlice(1).plotB0ppm      = [];   end
if ~isfield(Seq.AQSlice(1), 'plotB0PpmPhase'),      Seq.AQSlice(1).plotB0PpmPhase = [];   end
if ~isfield(Seq.AQSlice(1), 'plotB0ppmGradient'),   Seq.AQSlice(1).plotB0ppmGradient = []; end
if ~isfield(Seq.AQSlice(1), 'plotB0Gradient'),      Seq.AQSlice(1).plotB0Gradient = []; end
if isemptyfield(Seq.AQSlice(1), 'plotB0ppmhAxes'),  Seq.AQSlice(1).plotB0ppmhAxes = cell(1, 3); end
if isemptyfield(Seq.AQSlice(1), 'plotB0GradientppmhAxes'), Seq.AQSlice(1).plotB0GradienthAxes = cell(1, 3); end
if ~isfield(Seq.AQSlice(1), 'plotB0Hz'),            Seq.AQSlice(1).plotB0Hz       = [];   end
if ~isfield(Seq.AQSlice(1), 'plotB0HzPhase'),       Seq.AQSlice(1).plotB0HzPhase  = [];   end
if ~isfield(Seq.AQSlice(1), 'plotB0HzGradient'),    Seq.AQSlice(1).plotB0HzGradient = []; end
if ~isfield(Seq.AQSlice(1), 'plotB0GradientPhase'), Seq.AQSlice(1).plotB0GradientPhase = []; end
if isemptyfield(Seq.AQSlice(1), 'plotFft1_data'),   Seq.AQSlice(1).plotFft1_data    = []; end
if isemptyfield(Seq.AQSlice(1), 'plotData'),        Seq.AQSlice(1).plotData    = []; end
if isemptyfield(Seq.AQSlice(1), 'plotB0HzhAxes'),   Seq.AQSlice(1).plotB0HzhAxes  = cell(1, 3); end
if ~isfield(Seq.AQSlice(1), 'plotB1percent'),       Seq.AQSlice(1).plotB1percent  = [];   end
if isemptyfield(Seq.AQSlice(1), 'plotB1hAxes'),     Seq.AQSlice(1).plotB1hAxes    = cell(1); end
if ~isfield(Seq.AQSlice(1), 'RoiRelativeValue'),    Seq.AQSlice(1).RoiRelativeValue = []; end
if ~isfield(Seq.AQSlice(1), 'RoiCutOffPercentile'), Seq.AQSlice(1).RoiCutOffPercentile = []; end
if ~isfield(Seq.AQSlice(1), 'ZeroFillFactor'),      Seq.AQSlice(1).ZeroFillFactor = [];   end
if ~isfield(Seq.AQSlice(1), 'ZeroFillWindowSize'),  Seq.AQSlice(1).ZeroFillWindowSize = []; end
if ~isfield(Seq.AQSlice(1), 'sliceomaticProps'),    Seq.AQSlice(1).sliceomaticProps = []; end
if isemptyfield(Seq.AQSlice(1), 'raiseFigures'),    Seq.AQSlice(1).raiseFigures   = 0;    end
if isemptyfield(Seq.AQSlice(1), 'SliceCartesianAxis'), Seq.AQSlice(1).SliceCartesianAxis = cell(1); end
if isemptyfield(Seq.AQSlice(1), 'ReadCartesianAxis'), Seq.AQSlice(1).ReadCartesianAxis = cell(1); end
if isemptyfield(Seq.AQSlice(1), 'PhaseCartesianAxis'), Seq.AQSlice(1).PhaseCartesianAxis = cell(3,1); end

% slice orientation
if isemptyfield(Seq.AQSlice(1), 'Gamma'), Seq.AQSlice(1).Gamma = HW.GammaDef; end % Gamma in RAD
if isemptyfield(Seq.AQSlice(1), 'alfa'), Seq.AQSlice(1).alfa = 0; end
if isemptyfield(Seq.AQSlice(1), 'phi'), Seq.AQSlice(1).phi = 0; end
if isemptyfield(Seq.AQSlice(1), 'theta'), Seq.AQSlice(1).theta = 0; end
if isemptyfield(Seq.AQSlice(1), 'angle2Turns'), Seq.AQSlice(1).angle2Turns = 1/(2*pi); end
if isemptyfield(Seq.AQSlice(1), 'SliceGradSign'), Seq.AQSlice(1).SliceGradSign = 1; end
if isemptyfield(Seq.AQSlice(1), 'ReadGradSign'), Seq.AQSlice(1).ReadGradSign = 1; end
if isemptyfield(Seq.AQSlice(1), 'ReadGradRephaseSign'), Seq.AQSlice(1).ReadGradRephaseSign = 1; end
if isemptyfield(Seq.AQSlice(1), 'RephaseReadEncoding'), Seq.AQSlice(1).RephaseReadEncoding = 1; end
if isemptyfield(Seq.AQSlice(1), 'RephasePhaseEncoding'), Seq.AQSlice(1).RephasePhaseEncoding = 1; end

% image size
if isemptyfield(Seq.AQSlice(1), 'sizeRead'),  Seq.AQSlice(1).sizeRead = 0; end
if isemptyfield(Seq.AQSlice(1), 'sizePhase'), Seq.AQSlice(1).sizePhase = [0,0,0]; end
if Seq.AQSlice(1).sizePhase(1) == 0
  Seq.AQSlice(1).sizePhase(1) = HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(2) ...
    - HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(1);
end
if numel(Seq.AQSlice(1).sizePhase) < 2 || Seq.AQSlice(1).sizePhase(2) == 0
  Seq.AQSlice(1).sizePhase(2) = HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(4) ...
    - HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(3);
end
if numel(Seq.AQSlice(1).sizePhase) < 3 || Seq.AQSlice(1).sizePhase(3) == 0
  Seq.AQSlice(1).sizePhase(3) = HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(6) ...
    - HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(5);
end
if Seq.AQSlice(1).sizeRead == 0
  Seq.AQSlice(1).sizeRead = HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(6) ...
    - HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(5);
end

% Resolution
if isemptyfield(Seq.AQSlice(1), 'Resolution')
  Seq.AQSlice(1).Resolution = [];
else
  switch numel(Seq.AQSlice(1).Resolution)
    case 1
      if ~isemptyfield(Seq.AQSlice(1), 'thickness') && Seq.AQSlice(1).thickness < 1000
        % 2D iso
        Seq.AQSlice(1).Resolution = [Inf, repmat(Seq.AQSlice(1).Resolution(1), 1, 2)];
      else
        % 1D 2D 3D CSI iso
        Seq.AQSlice(1).Resolution = repmat(Seq.AQSlice(1).Resolution(1), 1, 3);
      end
    case 2
      % 2D
      Seq.AQSlice(1).Resolution = [Inf, Seq.AQSlice(1).Resolution(1), Seq.AQSlice(1).Resolution(2)];
  end
  if Seq.AQSlice(1).sizePhase(1) < 1000
    Seq.AQSlice(1).nPhase(1) = round(Seq.AQSlice(1).sizePhase(1)/Seq.AQSlice(1).Resolution(1));
    Seq.AQSlice(1).sizePhase(1) = Seq.AQSlice(1).nPhase(1) * Seq.AQSlice(1).Resolution(1);
  else
    if isemptyfield(Seq.AQSlice(1), 'thickness')
      Seq.AQSlice(1).thickness = Seq.AQSlice(1).Resolution(1);
    end
    Seq.AQSlice(1).nPhase(1) = 1;
  end
  Seq.AQSlice(1).nPhase(2) = round(Seq.AQSlice(1).sizePhase(2) / Seq.AQSlice(1).Resolution(2));
  Seq.AQSlice(1).sizePhase(2) = Seq.AQSlice(1).nPhase(2) * Seq.AQSlice(1).Resolution(2);
  if Seq.AQSlice(1).sizePhase(3) < Seq.AQSlice(1).sizeRead
    Seq.AQSlice(1).nPhase(3) = round(Seq.AQSlice(1).sizePhase(3) / Seq.AQSlice(1).Resolution(3));
    Seq.AQSlice(1).sizePhase(3) = Seq.AQSlice(1).nPhase(3) * Seq.AQSlice(1).Resolution(3);
    Seq.AQSlice(1).nRead = 1;
  else
    Seq.AQSlice(1).nRead = round(Seq.AQSlice(1).sizeRead / Seq.AQSlice(1).Resolution(3));
    Seq.AQSlice(1).sizeRead = Seq.AQSlice(1).nRead * Seq.AQSlice(1).Resolution(3);
    Seq.AQSlice(1).nPhase(3) = 1;
  end
  Seq.AQSlice(1).nPhase(isinf(Seq.AQSlice(1).nPhase) | ...
                        isnan(Seq.AQSlice(1).nPhase) | ...
                        (Seq.AQSlice(1).nPhase==0)) = 1;
  Seq.AQSlice(1).nRead(isinf(Seq.AQSlice(1).nRead) | ...
                       isnan(Seq.AQSlice(1).nRead) | ...
                       (Seq.AQSlice(1).nRead)==0) = 1;
  Seq.AQSlice(1).sizePhase(isinf(Seq.AQSlice(1).sizePhase) | ...
                           isnan(Seq.AQSlice(1).sizePhase) | ...
                           (Seq.AQSlice(1).sizePhase)==0) = Inf;
  Seq.AQSlice(1).sizeRead(isinf(Seq.AQSlice(1).sizeRead) | ...
                          isnan(Seq.AQSlice(1).sizeRead) | ...
                          (Seq.AQSlice(1).sizeRead)==0) = Inf;
end

if isemptyfield(Seq.AQSlice(1), 'nRead')
  Seq.AQSlice(1).nRead = 32;
end
% dimensions of AQSlice.nPhase
if isemptyfield(Seq.AQSlice(1), 'nPhase')
  if isZTE
    Seq.AQSlice(1).nPhase = [32 32 1];
  else
    Seq.AQSlice(1).nPhase = [1 32 1];
  end
else
  if Seq.AQSlice(1).nPhase(1)==0
    Seq.AQSlice(1).nPhase(1) = 1;
  end
  if numel(Seq.AQSlice(1).nPhase)==1 || Seq.AQSlice(1).nPhase(2)==0
    Seq.AQSlice(1).nPhase(2) = 1;
  end
  if numel(Seq.AQSlice(1).nPhase)==2 || Seq.AQSlice(1).nPhase(3)==0
    Seq.AQSlice(1).nPhase(3) = 1;
  end
end
if isZTE && sum(Seq.AQSlice(1).nPhase > 1) < 2
  error('PD:sequence_Flash:No1dZTE', ...
    'ZTE measurements are only implemented for 2d and 3d.');
end

% dimensions of AQSlice.PhaseOS
if isemptyfield(Seq.AQSlice(1), 'PhaseOS')
  Seq.AQSlice(1).PhaseOS = [1 1 1];
else
  if Seq.AQSlice(1).PhaseOS(1)==0
    Seq.AQSlice(1).PhaseOS(1) = 1;
  end
  if numel(Seq.AQSlice(1).PhaseOS)==1 || Seq.AQSlice(1).PhaseOS(2)==0
    Seq.AQSlice(1).PhaseOS(2) = 1;
  end
  if numel(Seq.AQSlice(1).PhaseOS)==2 || Seq.AQSlice(1).PhaseOS(3)==0
    Seq.AQSlice(1).PhaseOS(3) = 1;
  end
end
Seq.AQSlice(1).PhaseOS = round(Seq.AQSlice(1).nPhase.*Seq.AQSlice(1).PhaseOS)./Seq.AQSlice(1).nPhase; % integer number of phase steps

Seq.AQSlice(1).sizePhase(Seq.AQSlice(1).nPhase==1) = Inf;
Seq.AQSlice(1).sizeRead(Seq.AQSlice(1).nRead==1) = Inf;

% (approximate) extension of the RoI for B0 map in meters
if isemptyfield(Seq.CorrectB0Read, 'RoIExtension')
  Seq.CorrectB0Read.RoIExtension = 0.15 * [Seq.AQSlice(1).sizeRead, Seq.AQSlice(1).sizePhase(1:2)];
end
if isscalar(Seq.CorrectB0Read.RoIExtension)
  Seq.CorrectB0Read.RoIExtension(1:3) = Seq.CorrectB0Read.RoIExtension;
end

if Seq.CorrectB0Read.Get && sum(Seq.AQSlice(1).nPhase~=1)~=2
  error('PD:sequence_Flash:CorrectB0ReadNot3d', ...
    'Measuring the B0 read error only works for 3d images (read-phase(1)-phase(2)).');
end

if isemptyfield(Seq.AQSlice(1), 'thickness'), Seq.AQSlice(1).thickness = Inf; end
if isemptyfield(Seq.AQSlice(1), 'SliceDephase'), Seq.AQSlice(1).SliceDephase = 0; end
Seq.AQSlice(1).SliceDephase = Seq.AQSlice(1).SliceDephase & (Seq.AQSlice(1).thickness>=1000);

% number of images in "gradient echo train" (Multi-Gradient-Echo)
if isemptyfield(Seq.AQSlice(1), 'nImages')
  if isemptyfield(Seq, 'tEcho')
    Seq.AQSlice(1).nImages = 1;
  else
    Seq.AQSlice(1).nImages = numel(Seq.tEcho);
  end
end
if Seq.AQSlice(1).nImages ~= round(Seq.AQSlice(1).nImages) || Seq.AQSlice(1).nImages < 1
  error('PD:sequence_Flash:NonIntegerImages', ...
    'Seq.AQSlice.nImages must be a positive integer.');
end

% EPI factor and segments
% number of k-lines per image
if isZTE
  % one k-line pointing towards each grid point on a N-dimensional cuboid surface
  numKLines = ...
    2*((Seq.AQSlice(1).nPhase(1)*Seq.AQSlice(1).PhaseOS(1)-1) + ...
       (Seq.AQSlice(1).nPhase(2)*Seq.AQSlice(1).PhaseOS(2)-1)) * ...
       Seq.AQSlice(1).nPhase(3)*Seq.AQSlice(1).PhaseOS(3) + ...
    (Seq.AQSlice(1).nPhase(3)>1) * ...
    2 * prod(Seq.AQSlice(1).nPhase(1:2).*Seq.AQSlice(1).PhaseOS(1:2)-2);
else
  numKLines = prod(Seq.AQSlice(1).nPhase) * prod(Seq.AQSlice(1).PhaseOS);
end
if isemptyfield(Seq.AQSlice(1), 'EPISegments'), Seq.AQSlice(1).EPISegments = numKLines; end
% limit to total number of k-lines
Seq.AQSlice(1).EPISegments = min([numKLines*Seq.AQSlice(1).nImages, Seq.AQSlice(1).EPISegments]);
% find next integer divisor
Seq.AQSlice(1).EPISegments = Seq.AQSlice(1).EPISegments - 1 + ...
  find(~rem(numKLines./(Seq.AQSlice(1).EPISegments:numKLines),1), 1, 'first');

% number of image k-lines per excitation
if isemptyfield(Seq.AQSlice(1), 'EPIFactor')
  Seq.AQSlice(1).EPIFactor = numKLines / Seq.AQSlice(1).EPISegments;
end
Seq.AQSlice(1).EPISegments = numKLines / Seq.AQSlice(1).EPIFactor;
if (Seq.AQSlice(1).EPISegments ~= round(Seq.AQSlice(1).EPISegments)) || ...
    (Seq.AQSlice(1).EPIFactor ~= round(Seq.AQSlice(1).EPIFactor))
  error('PD:sequence_Flash:IncompatibleEPIFactor', ...
    ['"AQSlice.EPIFactor" (%f) or "AQSlice.EPISegments" (%f) cannot be used ', ...
    'with that number of kLines (%d) to acquire.'], ...
    Seq.AQSlice(1).EPIFactor, Seq.AQSlice(1).EPISegments, numKLines);
end
if isemptyfield(Seq, 'SingletRep')
  Seq.SingletRep = (Seq.AQSlice(1).EPIFactor==1);
end

% use bipolar read gradients
if isemptyfield(Seq.AQSlice(1), 'EPIReadBipolar')
  Seq.AQSlice(1).EPIReadBipolar = false;
end

% rephase gradient at dephase of next tRep
if isemptyfield(Seq.AQSlice(1), 'EPIGradient')
  Seq.AQSlice(1).EPIGradient = Seq.AQSlice(1).EPIFactor > 1;
end
% rephase gradient in phase direction at dephase of next tRep
if isemptyfield(Seq.AQSlice(1), 'EPIGradientPhase')
  Seq.AQSlice(1).EPIGradientPhase = Seq.AQSlice(1).EPIGradient;
end
% rephase gradient in read direction at dephase of next tRep
if isemptyfield(Seq.AQSlice(1), 'EPIGradientRead')
  Seq.AQSlice(1).EPIGradientRead = Seq.AQSlice(1).EPIGradient && Seq.AQSlice(1).EPIReadBipolar;
end

if Seq.AQSlice(1).EPIFactor > 1 && Seq.CorrectPhase && Seq.CorrectPhaseSeparate
  error('PD:sequence_Flash:EPIWithCorrectPhaseSeparate', ...
    'Separate excitations for frequency tracking with EPI are not yet implemented.');
end

% acquire separate images consisting of only odd or even Echoes respectively
if isemptyfield(Seq.AQSlice(1), 'oddEvenEchoes')
  Seq.AQSlice(1).oddEvenEchoes = false;
end
Seq.AQSlice(1).oddEvenEchoes = Seq.AQSlice(1).oddEvenEchoes > 0;


if isemptyfield(Seq.AQSlice(1), 'kLineOrderType'), Seq.AQSlice(1).kLineOrderType = 'increasing';  end
if isemptyfield(Seq.AQSlice(1), 'kLineOrder')
  switch lower(Seq.AQSlice(1).kLineOrderType)
    case 'increasing'
      Seq.AQSlice(1).kLineOrder = 1:numKLines;

    case 'centerout'
      nKLinesPhase = (Seq.AQSlice(1).nPhase .* Seq.AQSlice(1).PhaseOS);
      kLineOrder1 = reshape(cumsum((0:nKLinesPhase(1)-1) .* (-1).^(0:nKLinesPhase(1)-1)) + floor(nKLinesPhase(1)/2)+1, [],1,1);
      kLineOrder2 = reshape((cumsum((0:nKLinesPhase(2)-1) .* (-1).^(0:nKLinesPhase(2)-1)) + floor(nKLinesPhase(2)/2))*numel(kLineOrder1), 1,[],1);
      kLineOrder3 = reshape((cumsum((0:nKLinesPhase(3)-1) .* (-1).^(0:nKLinesPhase(3)-1)) + floor(nKLinesPhase(3)/2))*numel(kLineOrder1)*numel(kLineOrder2), 1,1,[]);
      Seq.AQSlice(1).kLineOrder = reshape(bsxfun(@plus, bsxfun(@plus, kLineOrder1, kLineOrder2), kLineOrder3), 1,[]);

    case 'decreasing'
      Seq.AQSlice(1).kLineOrder = numKLines:-1:1;

    case 'random'
      Seq.AQSlice(1).kLineOrder = randperm(numKLines);

    otherwise
      error('PD:sequence_Flash:incompatiblekLineOrderType', ...
      'AQSlice.kLineOrder must be empty, "increasing", "centerOut", "decreasing", or "random".');

  end
else
  Seq.AQSlice(1).kLineOrder = unique(reshape(Seq.AQSlice(1).kLineOrder, 1, []), 'stable');
  Seq.AQSlice(1).kLineOrderType = 'manual';
  if numel(Seq.AQSlice(1).kLineOrder) > numKLines || ...
      ~isreal(Seq.AQSlice(1).kLineOrder) || ...
      any(Seq.AQSlice(1).kLineOrder < 1) || any(Seq.AQSlice(1).kLineOrder > numKLines) || ...
      any(abs(Seq.AQSlice(1).kLineOrder - round(Seq.AQSlice(1).kLineOrder)) > eps(numKLines))
    % FIXME: Is there a better way to check that "kLineOrder" is a
    % permutation or selection vector?
    error('PD:sequence_Flash:incompatiblekLineOrder', ...
      '"AQSlice.kLineOrder" must be a permutation or selection vector with the numbers of kLines to be acquired.');
  end
end


if isemptyfield(Seq, 'tEcho'),
  if isemptyfield(Seq.AQSlice(1), 'EPIEchoSpacing')
    Seq.tEcho = 5e-3;
  else
    if Seq.AQSlice(1).nImages > 1
      Seq.tEcho = Seq.AQSlice(1).EPIEchoSpacing * ...
        Seq.AQSlice(1).EPIFactor;
    else
      Seq.tEcho = Seq.AQSlice(1).EPIEchoSpacing * ...
        (floor(Seq.AQSlice(1).EPIFactor/2) + 1);
    end
  end
end

if Seq.AQSlice(1).nImages > 1 && isscalar(Seq.tEcho)
  Seq.tEcho = (1:Seq.AQSlice(1).nImages) * Seq.tEcho;
end

if numel(Seq.tEcho) ~= Seq.AQSlice(1).nImages
  error('PD:sequence_Flash:nImages', ...
    'Seq.AQSlice(1).nImages does not match number of echo times in Seq.tEcho.');
end

% spacing between echoes in each EPI segment
if isemptyfield(Seq.AQSlice(1), 'EPIEchoSpacing')
  if Seq.AQSlice(1).nImages > 1
    Seq.AQSlice(1).EPIEchoSpacing = min(diff([0; Seq.tEcho(:)])) ...
      / Seq.AQSlice(1).EPIFactor;
  else
    Seq.AQSlice(1).EPIEchoSpacing = Seq.tEcho ...
      / (floor(Seq.AQSlice(1).EPIFactor/2) + 1);
  end
end

if Seq.AQSlice(1).EPIFactor == 1 && Seq.AQSlice(1).nImages == 1 && ~Seq.AQSlice(1).oddEvenEchoes
  % set to zero if it is not used
  Seq.AQSlice(1).EPIEchoSpacing = 0;
end


% if ~isfield(Seq.AQSlice(1), 'tEcho'), Seq.AQSlice(1).tEcho = []; end
if isemptyfield(Seq.AQSlice(1), 'tEchoOffset')
  if isfield(Seq, 'tEchoOffset')  % backwards compatibility
    Seq.AQSlice(1).tEchoOffset = Seq.tEchoOffset;
  else
    Seq.AQSlice(1).tEchoOffset = 0;
  end
end

if isemptyfield(Seq, 'RepetitionTime')  % repetition time of excitations in seconds
  if isemptyfield(Seq, 'tRep')
    Seq.RepetitionTime = 2*Seq.tEcho;
  else
    Seq.RepetitionTime = Seq.tRep(1);
  end
end
if isemptyfield(Seq, 'tRep')
  Seq.tRep = Seq.RepetitionTime;
end

if (Seq.AQSlice(1).EPIFactor > 1 || Seq.AQSlice(1).nImages > 1) && ...
    Seq.RepetitionTime < ...
    Seq.tEcho(end) + ...
    ceil(Seq.AQSlice(1).EPIFactor/2) * Seq.AQSlice(1).EPIEchoSpacing
  error('PD:sequence_Flash:RepetitionTimeTooShort', ...
    'Repetition time is too short for selected gradient echo train duration.');
end

if (Seq.AQSlice(1).EPIFactor > 1 || Seq.AQSlice(1).nImages > 1) && ...
    any(diff([0; Seq.tEcho(:)]) + eps(Seq.AQSlice(1).EPIEchoSpacing) < ...
    (floor(Seq.AQSlice(1).EPIFactor/2) + 1) * Seq.AQSlice(1).EPIEchoSpacing)
  error('PD:sequence_Flash:EchoTimeTooShort', ...
    'Echo time (spacing) is too short for selected EPI echo spacing and EPI factor.');
end

if isemptyfield(Seq, 'T1'), Seq.T1 = 100e-3; end
if isemptyfield(Seq, 'FlipAngle'), Seq.FlipAngle = acosd(exp(-Seq.RepetitionTime/Seq.T1)); end

if ~isfield(Seq.AQSlice(1), 'SpoilFactor') || isempty(Seq.AQSlice(1).SpoilFactor);
  if Seq.AQSlice(1).nRead > 1
    Seq.AQSlice(1).SpoilFactor = [0, 0, Inf]; % Spoil read direction
  else
    Seq.AQSlice(1).SpoilFactor = [1, 1, 1]; % Spoil all phase directions
  end
  if Seq.AQSlice(1).nPhase(1)==1 && Seq.AQSlice(1).thickness<1000
    Seq.AQSlice(1).SpoilFactor(1) = 8;
  end
end
if ~isZTE
  % FIXME: Spoilers for ZTE?
  if isemptyfield(Seq.AQSlice(1), 'sizePhaseSpoil'),
    if Seq.AQSlice(1).nRead > 1
      % Spoil read direction
      Seq.AQSlice(1).sizePhaseSpoil = ...
        [Seq.AQSlice(1).sizePhase(1)/Seq.AQSlice(1).nPhase(1), ...
        Seq.AQSlice(1).sizePhase(2)/Seq.AQSlice(1).nPhase(2), ...
        Seq.AQSlice(1).sizeRead/Seq.AQSlice(1).nRead] ./ Seq.AQSlice(1).SpoilFactor./2;
      Seq.AQSlice(1).sizePhaseSpoil(isnan(Seq.AQSlice(1).sizePhaseSpoil)) = Inf;
    else
      % Spoil all phase directions
      Seq.AQSlice(1).sizePhaseSpoil = ...
        [Seq.AQSlice(1).sizePhase(1)/Seq.AQSlice(1).nPhase(1), ...
        Seq.AQSlice(1).sizePhase(2)/Seq.AQSlice(1).nPhase(2), ...
        Seq.AQSlice(1).sizePhase(3)/Seq.AQSlice(1).nPhase(3)] ./ Seq.AQSlice(1).SpoilFactor./2;
    end
    if Seq.AQSlice(1).nPhase(1)==1 && Seq.AQSlice(1).thickness<1000
      Seq.AQSlice(1).sizePhaseSpoil(1) = Seq.AQSlice(1).thickness./Seq.AQSlice(1).SpoilFactor(1);
    end
  end

  if size(Seq.AQSlice(1).sizePhaseSpoil)==1
    Seq.AQSlice(1).sizePhaseSpoil = repmat(Seq.AQSlice(1).sizePhaseSpoil(1), 1, 3);
  end
end

if isemptyfield(Seq.AQSlice(1), 'DephaseLengthFactor')
  if isZTE
    % There are no dephase gradient pulses in ZTE sequences. Set their (default)
    % length factor to Inf to avoid (bogus) error messages about gradient ramps
    % not matching inside the pulse duration.
    Seq.AQSlice(1).DephaseLengthFactor = Inf;
  else
    Seq.AQSlice(1).DephaseLengthFactor = 1;
  end
end
if isemptyfield(Seq.AQSlice(1), 'SpoilLengthFactor'), Seq.AQSlice(1).SpoilLengthFactor = Seq.AQSlice(1).DephaseLengthFactor; end
if isemptyfield(Seq.AQSlice(1), 'DephaseTimeOffset'), Seq.AQSlice(1).DephaseTimeOffset = 0; end
if isemptyfield(Seq.AQSlice(1), 'SpoilTimeOffset'), Seq.AQSlice(1).SpoilTimeOffset = Seq.AQSlice(1).DephaseTimeOffset; end
if isemptyfield(Seq.AQSlice(1), 'excitationPulse')
  if Seq.AQSlice(1).thickness >= 1000
    Seq.AQSlice(1).excitationPulse = @Pulse_Rect;
  else
    Seq.AQSlice(1).excitationPulse = @Pulse_RaisedCos;
  end
end
if ~isfield(Seq.AQSlice(1), 'excitationFlipAngleComposite'), Seq.AQSlice(1).excitationFlipAngleComposite = []; end
if ~isfield(Seq.AQSlice(1), 'excitationFlipPhaseComposite'), Seq.AQSlice(1).excitationFlipPhaseComposite = []; end

if isemptyfield(Seq.AQSlice(1), 'Center2OriginImage')
  Seq.AQSlice(1).Center2OriginImage = [0, 0, 0];
end

if any(Seq.AQSlice(1).Center2OriginImage) && isZTE
  % FIXME: Add support for image shift?
  error('PD:sequence_Flash:ZTENoImageOffset', ...
    'For ZTE experiments, the image center cannot be shifted from the origin.');
end

if isemptyfield(Seq, 'SteadyState_PreShots')
  Seq.SteadyState_PreShots = 3 * double((prod(Seq.AQSlice(1).nPhase.*Seq.AQSlice(1).PhaseOS)>1) | (Seq.CorrectPhase~=0));
end
if isemptyfield(Seq, 'SteadyState_PostShots')
  Seq.SteadyState_PostShots = 4 * double(Seq.CorrectPhase~=0);
end

if Seq.CorrectPhase && Seq.CorrectPhaseSeparate && ...
    (Seq.SteadyState_PreShots == 0 || Seq.SteadyState_PostShots == 0)
  error('PD:sequence_Flash:CorrectPhaseSeparateNoPrePost', ...
    ['There must be at least one pre-shot (Seq.SteadyState_PreShots = %d) and ', ...
    'post-shot (Seq.SteadyState_PostShots = %d) when using Seq.CorrectPhaseSeparate'], ...
    Seq.SteadyState_PreShots, Seq.SteadyState_PostShots);
end

% settings for sequence plot
if ~isfield(Seq, 'plotSeqAQ'), Seq.plotSeqAQ = []; end
if ~isempty(Seq.plotSeqAQ)
  Seq.plotSeq = Seq.plotSeqAQ;
  if Seq.CorrectPhaseSeparate
    nCorrectPhaseSeparate = floor((Seq.AQSlice(1).EPISegments-1) / Seq.CorrectPhaseBlockSize);
  else
    nCorrectPhaseSeparate = 0;
  end
  if Seq.SingletRep
    if isemptyfield(Seq, 'plotSeqStart')
      Seq.plotSeqStart = Seq.SteadyState_PreShots + 1;
    end
    if isemptyfield(Seq, 'plotSeqEnd')
      Seq.plotSeqEnd = Seq.plotSeqStart + Seq.AQSlice(1).EPISegments ...
        + nCorrectPhaseSeparate - 1;
    end
  else
    if isemptyfield(Seq, 'plotSeqStart')
      Seq.plotSeqStart = Seq.SteadyState_PreShots * ...
        (Seq.AQSlice(1).EPIFactor * (Seq.AQSlice(1).oddEvenEchoes+1) * Seq.AQSlice(1).nImages+1) ...
        + 1;
    end
    if isemptyfield(Seq, 'plotSeqEnd')
      Seq.plotSeqEnd = Seq.plotSeqStart ...
        + (Seq.AQSlice(1).EPISegments + nCorrectPhaseSeparate) ...
          * (Seq.AQSlice(1).EPIFactor * (Seq.AQSlice(1).oddEvenEchoes+1) * Seq.AQSlice(1).nImages + 1) ...
        - 1;
    end
  end
  if ~isfield(Seq, 'plotSequence'), Seq.plotSequence = struct(); end
  if isemptyfield(Seq.plotSequence, 'wraps')
    Seq.plotSequence.wraps = Seq.AQSlice(1).EPISegments + nCorrectPhaseSeparate;
  end
  if isemptyfield(Seq.plotSequence, 'xLim')
    if isZTE
      % FIXME: Calculate better auto-limits
      Seq.plotSequence.xLim = [-Inf, 2/Seq.AQSlice(1).HzPerPixMin];
    else
      Seq.plotSequence.xLim = [-Inf, ...
        (min([Seq.tEcho(end) ...
                + (Seq.AQSlice(1).EPIEchoSpacing + Seq.tEcho(1)*(Seq.AQSlice(1).EPIEchoSpacing==0)) * ...
                  ((floor(Seq.AQSlice(1).EPIFactor/2) + 1) ...
                   + Seq.AQSlice(1).SpoilLengthFactor - Seq.AQSlice(1).DephaseLengthFactor) ...
                + Seq.AQSlice(1).SpoilTimeOffset, ...
              Seq.RepetitionTime]))];
    end
  end
end

if isemptyfield(Seq, 'plotZTEEndPoints')
  % show plots with end points of rays in k-space
  Seq.plotZTEEndPoints = false;
end

Seq = set_EmptyField(Seq, 'LoopSaveAllData', 0);
Seq = set_EmptyField(Seq, 'Find_Frequency_interval', HW.FindFrequencySweep.maxTime);

%% setup types of loops
CSRLoops = repmat([1,2], 1, Seq.CorrectSliceRephase);
Seq.LoopName = repmat({'CSRLoopPlus', 'CSRLoopMinus'}, 1, Seq.CorrectSliceRephase);

CRRLoops = ones(1, Seq.CorrectReadRephase);
Seq.LoopName = [Seq.LoopName, repmat({'CRRLoop'}, 1, Seq.CorrectReadRephase)];

CB0Loops = repmat([1,2], 1, Seq.CorrectB0Read.Get);
Seq.LoopName = [Seq.LoopName, repmat({'B0map_tEcho1', 'B0map_tEcho2'}, 1, Seq.CorrectB0Read.Get)];

Seq.LoopName = [Seq.LoopName, repmat({'normal'}, 1, Seq.Loops)];

if isempty(Grad)
  %% Gradient parameters
  % Grad(1) x
  % Grad(2) y
  % Grad(3) z
  % Grad(4) B0

  % time the values in Grad(t).Amp with the corresponding index are set
  [Grad(1:HW.Grad(Seq.AQSlice(1).iDevice).n).Time] = deal([]);
  % amplitude of the gradient at the time Grad(t).Time (linearly interpolated)
  [Grad(1:HW.Grad(Seq.AQSlice(1).iDevice).n).Amp] = deal([]);
  % additional shim; magnet shim is already considered in HW.MagnetShim.
  % Caution: Using high values over a long time will damage the gradient
  % coils and amplifiers!
  [Grad(1:HW.Grad(Seq.AQSlice(1).iDevice).n).Shim] = deal([]);
  % if there is a one in the array the gradients of the prior TR are used.
  % (reduces data traffic between the PC and the MRI device)
  [Grad(1:HW.Grad(Seq.AQSlice(1).iDevice).n).Repeat] = deal([]);
end

%% initial values for all loops
init.AQ = AQ;
init.TX = TX;
init.Seq = Seq;
init.Grad = Grad;
LoopCount = 0;
dataB0 = Seq.dataB0;

%% loops
for Loop = [CSRLoops, CRRLoops, CB0Loops, 1:Seq.Loops]
  LoopCount = LoopCount+1;
  AQ = init.AQ;
  TX = init.TX;
  Grad = init.Grad;
  Seq = init.Seq;
  Seq.Loop = Loop;
  Seq.LoopNameCount = LoopCount;

  if exist('SeqLoop', 'var')
    if isfield(SeqLoop.data, 'SliceGradTimeIntegralRephaseOffset')
      Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset = ...
        SeqLoop.AQSlice(1).SliceGradTimeIntegralRephaseOffset + SeqLoop.data.SliceGradTimeIntegralRephaseOffset;
    else
      Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset = ...
        SeqLoop.AQSlice(1).SliceGradTimeIntegralRephaseOffset;
    end
  end

  switch Seq.LoopName{Seq.LoopNameCount}
    case 'CSRLoopPlus'
      Seq.AQSlice(1).SliceGradSign = 1;

      Seq.AQSlice(1).ReadGradTimeIntegralOffset = 0;
      Seq.AQSlice(1).ReadGradTimeDelayOffset = 0;

      Seq.AQSlice(1).sizePhase(1:3) = Inf;      % size in phase(1:3)
      Seq.AQSlice(1).nPhase(1:3) = 1;           % Number of Pixels in phase(1:3)
      Seq.AQSlice(1).PhaseOS(1) = 4;            % Oversampling phase(1)  1...
      Seq.AQSlice(1).PhaseOS(2:3) = 1;          % Oversampling phase(2:3)  1...
      Seq.AQSlice(1).ReadCoordinate = Seq.AQSlice(1).SliceCoordinate;  % direction of read set to direction of slice
      Seq.AQSlice(1).EPIFactor = 1;             % number of echoes per excitation
      Seq.AQSlice(1).EPISegments = prod(Seq.AQSlice(1).nPhase) * prod(Seq.AQSlice(1).PhaseOS) / Seq.AQSlice(1).EPIFactor;

      if Seq.LoopSeqPlot
        % The pulse program layout is different from the "normal" loop.
        % Delete fields that might prevent displaying the pulse program with
        % plotSeq.
        Seq.plotSequence = [];
        Seq.plotSeqStart = [];
        Seq.plotSeqEnd = [];
      else
        Seq.CorrectPlotFrequency = 0;
        Seq.plotSeqTR = [];                     % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeq = [];                       % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
      end

    case 'CSRLoopMinus'
      Seq.AQSlice(1).SliceGradSign = -1;
      Seq.AQSlice(1).sizePhase(1:3) = Inf;      % size in phase(1:3)
      Seq.AQSlice(1).nPhase(1:3) = 1;           % Number of Pixels in phase(1:3)
      Seq.AQSlice(1).PhaseOS(1) = 4;            % Oversampling phase(1)  1...
      Seq.AQSlice(1).PhaseOS(2:3) = 1;          % Oversampling phase(2:3)  1...
      Seq.AQSlice(1).ReadCoordinate = Seq.AQSlice(1).SliceCoordinate;  % direction of read set to direction of slice
      Seq.AQSlice(1).EPIFactor = 1;             % number of echoes per excitation
      Seq.AQSlice(1).EPISegments = prod(Seq.AQSlice(1).nPhase) * prod(Seq.AQSlice(1).PhaseOS) / Seq.AQSlice(1).EPIFactor;

      if Seq.LoopSeqPlot
        % The pulse program layout is different from the "normal" loop.
        % Delete fields that might prevent displaying the pulse program with
        % plotSeq.
        Seq.plotSequence = [];
        Seq.plotSeqStart = [];
        Seq.plotSeqEnd = [];
      else
        Seq.CorrectPlotFrequency = 0;
        Seq.plotSeqTR = [];                     % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeq = [];                       % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
      end

    case 'CRRLoop'
      if exist('SeqLoop', 'var')
        if isfield(SeqLoop.data, 'ReadGradTimeIntegralOffset')
          Seq.AQSlice(1).ReadGradTimeIntegralOffset = ...
            SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset + SeqLoop.data.ReadGradTimeIntegralOffset;
          Seq.AQSlice(1).ReadGradTimeDelayOffset = ...
            SeqLoop.AQSlice(1).ReadGradTimeDelayOffset + SeqLoop.data.ReadGradTimeDelayOffset;
        else
          Seq.AQSlice(1).ReadGradTimeIntegralOffset = ...
            SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset;
          Seq.AQSlice(1).ReadGradTimeDelayOffset = ...
            SeqLoop.AQSlice(1).ReadGradTimeDelayOffset;
        end
      end

      Seq.AQSlice(1).nPhase(1:3) = 1;           % Number of Pixels in phase(1:3)
      Seq.AQSlice(1).PhaseOS(1) = 4;            % Oversampling phase(1)  1...
      Seq.AQSlice(1).PhaseOS(2:3) = 1;          % Oversampling phase(2:3)  1...

      % convert an EPI echo train to an images train
      Seq.AQSlice(1).nImages = Seq.AQSlice(1).nImages * Seq.AQSlice(1).EPIFactor;
      Seq.tEcho = Seq.tEcho(1) + ((1:Seq.AQSlice(1).EPIFactor) - Seq.AQSlice(1).EPIFactor/2 - 1) * Seq.AQSlice(1).EPIEchoSpacing;
      Seq.AQSlice(1).EPIFactor = 1;             % number of echoes per excitation
      Seq.AQSlice(1).EPISegments = prod(Seq.AQSlice(1).nPhase) * prod(Seq.AQSlice(1).PhaseOS);
      Seq.AQSlice(1).partsAverage = 0;

      if Seq.LoopSeqPlot
        % The pulse program layout is different from the "normal" loop.
        % Delete fields that might prevent displaying the pulse program with
        % plotSeq.
        Seq.plotSequence = [];
        Seq.plotSeqStart = [];
        Seq.plotSeqEnd = [];
      else
        Seq.CorrectPlotFrequency = 0;
        Seq.plotSeqTR = [];                     % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeq = [];                       % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
      end

    case 'B0map_tEcho1'
      if Seq.AQSlice(1).nRead < 2
        continue;
      end
      if Seq.CorrectB0Read.tEchoIncr == 0
        warning('PD:sequence_Flash:B0ZeroEchoIncr', ...
          'Seq.CorrectB0Read.tEchoIncr must not equal 0 for B0map_tEcho. Skipping B0 measurement...')
        continue;
      end

      Seq.tEcho = Seq.tEcho(1);
      % Seq.AQSlice(1).EPIFactor = 1;             % number of echoes per excitation
      % Seq.AQSlice(1).EPISegments = prod(Seq.AQSlice(1).nPhase) * prod(Seq.AQSlice(1).PhaseOS) / Seq.AQSlice(1).EPIFactor;
      % Seq.AQSlice(1).EPIEchoSpacing = 0;
      Seq.AQSlice(1).nImages = 1;               % number of images (k-lines) acquired in one gradient echo train

      if ~isinf(HW.FindFrequencySweep.maxTime)
        Seq.Find_Frequency_interval = 0;
      end
      Seq.AQSlice(1).ZeroFillWindowSize = Seq.CorrectB0Read.ZeroFillWindowSize;
      Seq.CorrectB0Read.Use = false;

    case 'B0map_tEcho2'
      if Seq.AQSlice(1).nRead < 2 || Seq.CorrectB0Read.tEchoIncr == 0
        continue;
      end

      Seq.tEcho = Seq.tEcho(1);
      Seq.AQSlice(1).tEchoOffset = Seq.AQSlice(1).tEchoOffset + Seq.CorrectB0Read.tEchoIncr;
      % Seq.AQSlice(1).EPIFactor = 1;             % number of echoes per excitation
      % Seq.AQSlice(1).EPISegments = prod(Seq.AQSlice(1).nPhase) * prod(Seq.AQSlice(1).PhaseOS) / Seq.AQSlice(1).EPIFactor;
      % Seq.AQSlice(1).EPIEchoSpacing = 0;
      Seq.AQSlice(1).nImages = 1;               % number of images (k-lines) acquired in one gradient echo train

      Seq.Find_Frequency_interval = Inf;  % FIXME: Find frequency again for second measurement?
      Seq.AQSlice(1).ZeroFillWindowSize = Seq.CorrectB0Read.ZeroFillWindowSize;
      Seq.CorrectB0Read.Use = false;

    case 'normal'
      if exist('SeqLoop', 'var')
        if isfield(SeqLoop.data, 'ReadGradTimeIntegralOffset')
          Seq.AQSlice(1).ReadGradTimeIntegralOffset = ...
            SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset + SeqLoop.data.ReadGradTimeIntegralOffset;
          Seq.AQSlice(1).ReadGradTimeDelayOffset = ...
            SeqLoop.AQSlice(1).ReadGradTimeDelayOffset + SeqLoop.data.ReadGradTimeDelayOffset;
        else
          Seq.AQSlice(1).ReadGradTimeIntegralOffset = ...
            SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset;
          Seq.AQSlice(1).ReadGradTimeDelayOffset = ...
            SeqLoop.AQSlice(1).ReadGradTimeDelayOffset;
        end
      end

    otherwise
      % FIXME: Issue warning?

  end

  if isemptyfield(Seq, 'SweepPhase')
    % foldout
    Seq.SweepPhase = ...
      round(prod(Seq.AQSlice(1).nPhase.*Seq.AQSlice(1).PhaseOS) / ...
            max([2, Seq.AQSlice(1).PhaseOS(find(Seq.AQSlice(1).nPhase>1, 1, 'first'))]));
  end

  %% construct pulse sequence
  % Programming notes:
  % The excitation pulses (with slice selection) are in a separate tRep. The
  % pulse center is at t=0 of these tReps.
  % The acquisition windows and corresponding encoding gradient pulses share a
  % tRep. The gradient echo center is at t=0 of these tReps.
  % All tReps have the same length of Seq.AQSlice(1).EPIEchoSpacing. But the
  % tRep before each EPI segment is adjusted such the the effective echo time is
  % Seq.tEcho. Additionally, the last tRep of the last EPI segment has a length
  % such that the time between excitation pulses is Seq.RepetitionTime.

  % assign acquisition windows of each tRep to kLine of all images
  % number of images (real and for corrections) acquired per echo train
  numImages = Seq.AQSlice(1).nImages * (Seq.AQSlice(1).oddEvenEchoes+1);
  % number of k-lines per image (has to be repeated from above for the
  % correction loops)
  if isZTE
    % one k-line pointing towards each grid point on a N-dimensional cuboid surface
    numKLines = ...
      2*((Seq.AQSlice(1).nPhase(1)*Seq.AQSlice(1).PhaseOS(1)-1) + ...
         (Seq.AQSlice(1).nPhase(2)*Seq.AQSlice(1).PhaseOS(2)-1)) * ...
         Seq.AQSlice(1).nPhase(3)*Seq.AQSlice(1).PhaseOS(3) + ...
      (Seq.AQSlice(1).nPhase(3)>1) * ...
      2 * prod(Seq.AQSlice(1).nPhase(1:2).*Seq.AQSlice(1).PhaseOS(1:2)-2);
  else
    numKLines = prod(Seq.AQSlice(1).nPhase) * prod(Seq.AQSlice(1).PhaseOS);
  end
  % assign number (in phase increment order) to each k-line that will be
  % acquired in this experiment and break dimensions such that
  % oddEven x nImages x EPI segments x EPI factor
  kLinesDim = reshape(1:numImages*numKLines, ...
    Seq.AQSlice(1).oddEvenEchoes+1, Seq.AQSlice(1).nImages, ...
    [], Seq.AQSlice(1).EPIFactor);
  if any(strcmp(Seq.LoopName{Seq.LoopNameCount}, {'normal', 'B0map_tEcho1', 'B0map_tEcho2'}))
    % reshape such that EPI segments and EPI trains share the same dimension
    kLinesDimTemp = reshape(kLinesDim, ...
      (Seq.AQSlice(1).oddEvenEchoes+1), Seq.AQSlice(1).nImages, ...
      numKLines);
    % permute the order of the k-lines and reshape to the original dimensions
    kLinesDimTemp = kLinesDimTemp(:,:,[Seq.AQSlice(1).kLineOrder, setdiff(1:numKLines, Seq.AQSlice(1).kLineOrder)]);
    kLinesDimTemp = kLinesDimTemp(:,:,1:numel(Seq.AQSlice(1).kLineOrder),:);
    % consider reduced size if some k-lines are skipped
    szRed = [size(kLinesDim, 1), size(kLinesDim, 2), size(kLinesDim, 3), size(kLinesDim, 4)];
    szRed(4) = szRed(4) / (numel(kLinesDim)/numel(kLinesDimTemp));
    if szRed(4) < 1
      szRed(3) = szRed(3) * szRed(4);
      szRed(4) = 1;
    end
    if any(mod(szRed,1) ~= 0)
      error('PD:sequence_Echo_Standard:incompatibleReduction', ...
        'Number of acquired k-lines incompatible with number of EPI segments.');
    end
    kLinesDimTemp = reshape(kLinesDimTemp, szRed);
  else
    kLinesDimTemp = kLinesDim;

    % k-line order only applies for imaging loops
    Seq.AQSlice(1).kLineOrder = 1:(numel(kLinesDimTemp)/Seq.AQSlice(1).nImages);
  end

  % chronological order in which k-lines are to be acquired:
  % k-lines in EPI segments x EPI segments
  Seq.kLines = reshape(permute(kLinesDimTemp, [1, 4, 2, 3]), ...
    Seq.AQSlice(1).EPIFactor*Seq.AQSlice(1).nImages*(Seq.AQSlice(1).oddEvenEchoes+1), []);
  % order of those k-lines for image reconstruction:
  % kLines x nImages x corrections
  Seq.kLinesImages = reshape(permute(kLinesDim, [3, 4, 2, 1]), ...
    numKLines, Seq.AQSlice(1).nImages, []);

  if Seq.CorrectPhase && Seq.CorrectPhaseSeparate
    % calculate number of additional excitations for frequency tracking
    nCorrectPhaseSeparate = floor((size(Seq.kLines, 2)-1) / Seq.CorrectPhaseBlockSize);
  else
    nCorrectPhaseSeparate = 0;
  end

  % tReps with excitation pulse
  Seq.rfPulsetReps = cumsum(...
    ones(1, size(Seq.kLines, 2) + Seq.SteadyState_PreShots + Seq.SteadyState_PostShots + nCorrectPhaseSeparate) + ...
    size(Seq.kLines, 1)) - ...
    size(Seq.kLines, 1);
  % tReps with excitation pulses for image
  rfPulseImagetReps = Seq.rfPulsetReps(Seq.SteadyState_PreShots+1:end-Seq.SteadyState_PostShots);
  rfPulseFrequencytReps = [];
  if Seq.CorrectPhase && Seq.CorrectPhaseSeparate
    % select pulses that are used for frequency tracking
    rfPulseFrequencytReps = rfPulseImagetReps((Seq.CorrectPhaseBlockSize+1):(Seq.CorrectPhaseBlockSize+1):end);
    % remove pulses that are used for frequency tracking
    rfPulseImagetReps((Seq.CorrectPhaseBlockSize+1):(Seq.CorrectPhaseBlockSize+1):end) = [];
  end
  % tReps with acquisition window
  Seq.kLinestReps = repmat(rfPulseImagetReps, size(Seq.kLines, 1), 1) + ...
    repmat((1:size(Seq.kLines, 1)).', 1, numel(rfPulseImagetReps));
  % Matrix with the tReps grouped into EPI segments:
  % tReps in echo train x EPI segments
  Seq.tRepEPISegment = [rfPulseImagetReps; Seq.kLinestReps];

  % generate tRep times
  Seq.tRep = Seq.AQSlice(1).EPIEchoSpacing / (Seq.AQSlice(1).oddEvenEchoes+1) * ...
    ones(1, (Seq.AQSlice(1).EPIFactor*Seq.AQSlice(1).nImages*(Seq.AQSlice(1).oddEvenEchoes+1) + 1));
  % extend tRep before first EPI segment
  Seq.tRep(1) = Seq.tEcho(1) ...
    - (floor(Seq.AQSlice(1).EPIFactor/2) + (Seq.AQSlice(1).oddEvenEchoes)/4) * Seq.AQSlice(1).EPIEchoSpacing;
  % extend tReps between EPI segments
  nEPISegment = Seq.AQSlice(1).EPIFactor * (Seq.AQSlice(1).oddEvenEchoes+1);
  Seq.tRep(nEPISegment+1:nEPISegment:end-1) = ...
    diff(Seq.tEcho) - (Seq.AQSlice(1).EPIFactor-1+Seq.AQSlice(1).oddEvenEchoes/2) * Seq.AQSlice(1).EPIEchoSpacing;
  % extend last tRep in EPI segment to match RepetitionTime between pulses
  Seq.tRep(end) = Seq.RepetitionTime - ...
    (Seq.tEcho(end) ...  % "center" of EPI segment of last image
     + (ceil(Seq.AQSlice(1).EPIFactor/2)-1+Seq.AQSlice(1).oddEvenEchoes/4)*Seq.AQSlice(1).EPIEchoSpacing);
  % repeat for all EPI segments (as well as pre and post shots, and separate frequency tracking)
  Seq.tRep = repmat(Seq.tRep, 1, ...
    Seq.AQSlice(1).EPISegments + Seq.SteadyState_PreShots + Seq.SteadyState_PostShots + nCorrectPhaseSeparate);

  if isemptyfield(Seq.AQSlice(1), 'UsetRep')
    % k-lines x nImages x corrections
    [~, b] = sort(Seq.kLines(:));
    Seq.AQSlice(1).UsetRep = Seq.kLinestReps(b);
  end


  if Seq.StartSequence && (Loop==1  || Seq.LoopsBreakExactly==0)
    if isemptyfield(mySave, 'lastTime'), mySave.lastTime = 0; end
    if isinf(Seq.Find_Frequency_interval) || ...
        (now*24*3600-mySave.lastTime < Seq.Find_Frequency_interval-0.01)
      % Use last results of frequency sweep (i.e. mySave.HW.B0)
      [HW, mySave] = Find_Frequency_Sweep(HW, mySave, Inf);
    else
      if ~isempty(Seq.StartSequenceTime)
        sleep(max(0, Seq.StartSequenceTime - now*24*3600));
      end

      oldFindFrequencyPause = HW.FindFrequencyPause;
      hwGuard = onCleanup(@() setfield(HW, 'FindFrequencyPause', oldFindFrequencyPause));
      HW.FindFrequencyPause = 0;
      [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);  % Find magnet frequency
      delete(hwGuard);

      Seq.StartSequenceTime = now*24*3600 + max([Seq.LoopsBreak, HW.FindFrequencyPause]);
    end
  end


  if Loop == 1
    % % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad
    % % (1==x, 2==y, 3==z, 0 no Grad)
    % Seq.plotSeqTR = [1:3];
    % % Plot sequence on real timeline, plot RF, AQ and Grad
    % Seq.plotSeq = 1:3;
  elseif ~Seq.LoopSeqPlot
    % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad
    % (1==x, 2==y, 3==z, 0 no Grad)
    Seq.plotSeqTR = [];
    % Plot sequence on real timeline, plot RF, AQ and Grad
    Seq.plotSeq = [];
  end

  if ~isZTE
    % slice excitation by Seq.FlipAngle
    Seq.Slice(1).Pulse.Function = Seq.AQSlice(1).excitationPulse;
    Seq.Slice(1).Pulse.FlipAngleComposite = Seq.AQSlice(1).excitationFlipAngleComposite;
    Seq.Slice(1).Pulse.FlipPhaseComposite = Seq.AQSlice(1).excitationFlipPhaseComposite;
    Seq.Slice(1).Pulse.FlipAngle = Seq.FlipAngle;
    Seq.Slice(1).Thickness = Seq.AQSlice(1).thickness;
    % Seq.Slice(1).Pulse.Phase = 360*rand(1,numel(Seq.tRep));
    nPulses = Seq.SteadyState_PreShots + Seq.SteadyState_PostShots + Seq.AQSlice(1).EPISegments + nCorrectPhaseSeparate;
    Seq.Slice(1).Pulse.Phase = ...
      360*rem(cumsum(((1:nPulses) - Seq.SteadyState_PreShots)) ./ nPulses*Seq.SweepPhase, 1) ...
      + (0.5*Seq.SweepPhaseQuadratically *((1:nPulses).^2+(1:nPulses)+2));
    Seq.Slice(1).Pulse.PhaseIncrement = Seq.PhaseIncrement+Seq.PhaseIncrementLoop*(Loop-1);
    Seq.Slice(1).CenterOfPulse = 0;
    Seq.Slice(1).UseCoordinate = Seq.AQSlice(1).SliceCoordinate;
    Seq.Slice(1).GradTimeDelay = HW.Grad(Seq.AQSlice(1).iDevice).SliceTimeDelay;
    Seq.Slice(1).UseAtRepetitionTime = ((1:nPulses).'-1) * ...
      (Seq.AQSlice(1).EPIFactor * Seq.AQSlice(1).nImages * (Seq.AQSlice(1).oddEvenEchoes+1) + 1) + 1;
    Seq.Slice(1).distance = Seq.AQSlice(1).Center2OriginImage(1);
    Seq.Slice(1).MaxGradAmp = Seq.MaxGradAmpSlice;
    Seq.Slice(1).GradTimeIntegralRephaseOffset = Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset;
    Seq.Slice(1).GradSign = Seq.AQSlice(1).SliceGradSign(mod(Loop-1,numel(Seq.AQSlice(1).SliceGradSign))+1);
    Seq.Slice(1).GradDephaseSign = -Seq.AQSlice(1).SliceGradSign(mod(Loop-1,numel(Seq.AQSlice(1).SliceGradSign))+1);
    Seq.Slice(1).GradRephaseSign = -Seq.AQSlice(1).SliceGradSign(mod(Loop-1,numel(Seq.AQSlice(1).SliceGradSign))+1);
    if ~Seq.AQSlice(1).SliceDephase
      % In this case, this gradient block won't be used. Set duration to Inf to
      % avoid errors with gradient ramp checks.
      % FIXME: Would be nice to have a separate SliceDephaseLengthFactor
      Seq.Slice(1).GradDephaseLength = Inf;
    end
    if ~Seq.AQSlice(1).SliceRephaseImmediately && ...
        (Seq.AQSlice(1).thickness < 1000) && (Seq.CorrectPhase == 0 || Seq.CorrectPhaseSeparate)
      % In this case, this gradient block won't be used. Set duration to Inf to
      % avoid errors with gradient ramp checks.
      % FIXME: Would be nice to have a separate SliceRephaseLengthFactor
      Seq.Slice(1).GradRephaseLength = Inf;
    end

    % generate slice parameters
    Seq = get_SliceParameter(Seq, HW);

    if (Seq.AQSlice(1).EPIFactor > 1 || Seq.AQSlice(1).nImages > 1) && ...
        (Seq.RepetitionTime + Seq.Slice(1).GradCenter - Seq.Slice(1).GradLength - Seq.CLTime - 10e-6) ...
        < Seq.tEcho(end) + (floor(Seq.AQSlice(1).EPIFactor/2) + 0.5) * Seq.AQSlice(1).EPIEchoSpacing
      error('PD:sequence_Flash:TooManyImages', ...
        'Repetition time is too short to acquire selected number of images or longest (effective) echo time is too large.');
    end
  end

  CorrectPhaseDeadTimeTX2RX = 0;

  if Seq.CorrectPhase > 0
    Seq.AQSlice(2).plotkSpace = 1;                            % plot k-space
    Seq.AQSlice(2).plotImage = 1;                             % plot image
    Seq.AQSlice(2).plotPhase = 2;                             % plot phase of k-space and image
    Seq.AQSlice(2).Gamma = Seq.AQSlice(1).Gamma;              % gamma in radians
    Seq.AQSlice(2).alfa = 0;                                  % first rotation around x-axis in radians
    Seq.AQSlice(2).phi = 0;                                   % second rotation around y-axis in radians
    Seq.AQSlice(2).theta = 0;                                 % third rotation around z-axis in radians
    Seq.AQSlice(2).angle2Turns = Seq.AQSlice(1).angle2Turns;
    if isemptyfield(Seq.AQSlice(2), 'HzPerPixMin')
      Seq.AQSlice(2).HzPerPixMin = 1/Seq.CorrectPhaseDuration;  % bandwith per pixel in Hz
    end
    if isemptyfield(Seq.AQSlice(2), 'nRead')
      Seq.AQSlice(2).nRead = Seq.CorrectPhaseNRead;  % number of pixels in read direction
    end
    % Readout at FID
    Seq.AQSlice(2).nPhase(1) = 1;                             % number of pixels in phase(1) direction
    if Seq.CorrectPhaseSeparate
      Seq.AQSlice(2).nPhase(2) = 1;                       % number of pixels in phase(2) direction
    else
      Seq.AQSlice(2).nPhase(2) = nPulses;                       % number of pixels in phase(2) direction
    end
    Seq.AQSlice(2).nPhase(3) = 1;                             % number of pixels in phase(3) direction
    Seq.AQSlice(2).nPhase3D = 1;                              % number of pixels in phase3D (????)
    Seq.AQSlice(2).sizeRead = Inf;                            % image size in read in meter
    Seq.AQSlice(2).sizePhase(1) = Inf;                        % image size in phase(1) in meter
    Seq.AQSlice(2).sizePhase(2) = Inf;                        % image size in phase(2) in meter
    Seq.AQSlice(2).sizePhase(3) = Inf;                        % image size in phase(3) in meter
    Seq.AQSlice(2).thickness = Inf;                           % slice thickness in meter (used for 2D and 3D!)
    if isempty(Seq.AQSlice(2).iDevice), Seq.AQSlice(2).iDevice = Seq.AQSlice(1).iDevice; end
    autoReadOS = false;
    if isemptyfield(Seq.AQSlice(2), 'ReadOS')
      % integer oversampling factor in read direction (recommended >=2)
      Seq.AQSlice(2).ReadOS = max(2, ...
        ceil(HW.RX(Seq.AQSlice(2).iDevice).fSample / ...
          min(8000, HW.RX(Seq.AQSlice(2).iDevice).CIC_Decimation_Max) / ...
          (Seq.AQSlice(2).HzPerPixMin*Seq.AQSlice(2).nRead)));
      autoReadOS = true;
    end
    if isemptyfield(Seq.AQSlice(2), 'SamplingFactor')
      % integer sampling factor that is used to (down-)sample the received signal
      minReadOS = ...
        ceil(HW.RX(Seq.AQSlice(2).iDevice).fSample / ...
          min(8000, HW.RX(Seq.AQSlice(2).iDevice).CIC_Decimation_Max) / ...
          (Seq.AQSlice(2).HzPerPixMin*Seq.AQSlice(2).nRead)/Seq.AQSlice(2).ReadOS)*Seq.AQSlice(2).ReadOS;
      if minReadOS > Seq.AQSlice(2).ReadOS
        Seq.AQSlice(2).SamplingFactor = ceil(minReadOS/Seq.AQSlice(2).ReadOS);
      else
        Seq.AQSlice(2).SamplingFactor = 1;
      end
    end
    if autoReadOS
      Seq.AQSlice(2).ReadOS = ceil(Seq.AQSlice(2).ReadOS / Seq.AQSlice(2).SamplingFactor);
    end
    if ~Seq.CorrectPhaseSeparate
      CorrectPhaseDeadTimeTX2RX = get_DeadTimeTX2RX(HW, Seq.AQSlice(2).HzPerPixMin*Seq.AQSlice(2).nRead*Seq.AQSlice(2).ReadOS, Seq.AQSlice(2).iDevice);
    end
    Seq.AQSlice(2).PhaseOS(1) = 1;                            % oversampling factor in phase(1) direction: 1...
    Seq.AQSlice(2).PhaseOS(2) = 1;                            % oversampling factor in phase(2) direction: 1...
    Seq.AQSlice(2).PhaseOS(3) = 1;                            % oversampling factor in phase(3) direction: 1...
    if Seq.CorrectPhaseSeparate
      Seq.AQSlice(2).UsetRep = ...
        [Seq.rfPulsetReps(1:Seq.SteadyState_PreShots), rfPulseFrequencytReps, ...
         Seq.rfPulsetReps((end-Seq.SteadyState_PostShots+1):end)]+1;
    else
      Seq.AQSlice(2).UsetRep = Seq.Slice(1).UseAtRepetitionTime;  % directly after each excitation pulse
    end
    Seq.AQSlice(2).UseAQWindow = 1;
    Seq.AQSlice(2).ReadCoordinate = Seq.AQSlice(1).ReadCoordinate;  % direction of Read:  x = 1,  y = 2, z = 3
    Seq.AQSlice(2).PhaseCoordinate = Seq.AQSlice(1).PhaseCoordinate;  % direction of Phase:   x = 1,  y = 2, z = 3
    Seq.AQSlice(2).SliceCoordinate = Seq.AQSlice(1).SliceCoordinate;  % direction of Slice:   x = 1,  y = 2, z = 3
    Seq.AQSlice(2).SliceCartesianAxis = cell(1);
    Seq.AQSlice(2).ReadCartesianAxis = cell(1);
    Seq.AQSlice(2).PhaseCartesianAxis = cell(3,1);
    Seq.AQSlice(2).plotImageHandle = 109;
    if isempty(Seq.AQSlice(2).AmplitudeUnit), Seq.AQSlice(2).AmplitudeUnit = Seq.AQSlice(1).AmplitudeUnit; end
    if isempty(Seq.AQSlice(2).AmplitudeUnitScale), Seq.AQSlice(2).AmplitudeUnitScale = Seq.AQSlice(1).AmplitudeUnitScale; end
    if isempty(Seq.AQSlice(2).LengthUnit), Seq.AQSlice(2).LengthUnit = Seq.AQSlice(1).LengthUnit; end
    if isempty(Seq.AQSlice(2).LengthUnitScale), Seq.AQSlice(2).LengthUnitScale = Seq.AQSlice(1).LengthUnitScale; end
        
    if Seq.AQSlice(1).thickness < 1000
      Seq.CorrectPhaseAQtOffset = ...
        max(Seq.Slice(1).CenterOfRephase + Seq.Slice(1).GradRephaseLength/2 ...
            + Seq.CorrectPhaseAQtEC, ...
            Seq.CorrectPhaseAQtOffset);
    else
      Seq.CorrectPhaseAQtOffset = ...
        max(get_DeadTimeTX2RX(HW, Seq.AQSlice(2).HzPerPixMin*Seq.AQSlice(2).nRead*Seq.AQSlice(2).ReadOS, Seq.AQSlice(1).iDevice), ...
             Seq.CorrectPhaseAQtOffset);
    end

  end

  % acquisition at tEcho
  if isemptyfield(Seq.AQSlice(1), 'SpoilEddyCurrentTime')
    Seq.AQSlice(1).SpoilEddyCurrentTime = HW.Grad(Seq.AQSlice(1).iDevice).tEC;
  end
  if isemptyfield(Seq.AQSlice(1), 'DephasePreEddyCurrentTime')
    Seq.AQSlice(1).DephasePreEddyCurrentTime = max(HW.Grad(Seq.AQSlice(1).iDevice).TimeDelay(1:3)) * (1 + 0.5*(Seq.CorrectReadRephase>0));
  end
  if isemptyfield(Seq.AQSlice(1), 'AcquisitionTime')
    Seq.AQSlice(1).AcquisitionTime = Inf;
  end
  if isemptyfield(Seq.AQSlice(1), 'HzPerPixMin')
    Seq.AQSlice(1).HzPerPixMin = 1/Seq.AQSlice(1).AcquisitionTime;
  end
  if ~Seq.AQSlice(1).HzPerPixMin
    % calculate maximum possible acquisition time

    if isZTE
      error('PD:sequence_Flash:ZTEneedsHzPerPixMin', ...
        'For ZTE experiments, either HzPerPixMin or AcquisitionTime must be specified.');
    end

    % timing between excitation and start of dephase pulse
    tpre = Seq.tEcho(1) - (floor(Seq.AQSlice(1).EPIFactor / 2) + (Seq.AQSlice(1).oddEvenEchoes)/2) * Seq.AQSlice(1).EPIEchoSpacing/(Seq.AQSlice(1).oddEvenEchoes+1) ... % duration from center of rf pulse to center of first echo
      - abs(Seq.AQSlice(1).tEchoOffset) ...
      - (Seq.AQSlice(1).DephasePreEddyCurrentTime*1) * ((Seq.AQSlice(1).thickness>=1000) || (Seq.CorrectPhase>0 && ~Seq.CorrectPhaseSeparate)) ...
      - Seq.AQSlice(1).ReadGradTimeDelayOffset ...
      - (Seq.CorrectPhaseAQtOffset * (Seq.AQSlice(1).thickness<1000) + Seq.CorrectPhaseDuration * (1+1./Seq.CorrectPhaseNRead)) .* (Seq.CorrectPhase>0 && ~Seq.CorrectPhaseSeparate) ...  % Readout at FID (phase correction)
      - Seq.AQSlice(1).DephaseTimeOffset ...
      - max(Seq.Slice(1).Pulse.MaxLength*0.5 ...  % rf pulse
            + CorrectPhaseDeadTimeTX2RX .* (Seq.CorrectPhase>0 && ~Seq.CorrectPhaseSeparate) ...  % CorrectPhase acquisition DeadTimeTX2RX
            + (HW.TX(Seq.AQSlice(1).iDevice).DampCoil.DigitalOutputDuration) * HW.TX((Seq.AQSlice(1).iDevice)).DampCoil.Enable * (Seq.CorrectPhase==0 || Seq.CorrectPhaseSeparate) ...  % damp coil signal
            , Seq.CorrectPhaseAQtOffset * (Seq.CorrectPhase>0 && ~Seq.CorrectPhaseSeparate)) ...
        * (Seq.AQSlice(1).thickness>=1000) ...  % no slice selection
      - (Seq.Slice(1).CenterOfPulse + Seq.Slice(1).GradLength/2 + HW.Grad(Seq.AQSlice(1).iDevice).tEC) ...
        * (Seq.AQSlice(1).thickness<1000) * (~Seq.AQSlice(1).SliceRephaseImmediately) ...  % with slice selection, no phase correction
      - Seq.CLTime(min(2, numel(Seq.CLTime))) - 1e-6;


    % timing between end of rephase pulse and excitation of next repetition (potentially also slice gradient)
    if Seq.AQSlice(1).SliceDephase
      SliceDephaseStart=(-Seq.Slice(1).CenterOfDephase + Seq.Slice(1).GradDephaseLength/2 ...
          + max(HW.Grad(Seq.AQSlice(1).iDevice).TimeDelay(1:3)) );
    else
      SliceDephaseStart=0;
    end
    tpost = Seq.tRep(end) ...  % duration from center of last AQ window to center of next rf pulse (ignoring offsets)
      - abs(Seq.AQSlice(1).tEchoOffset) ...  % offset of acquisition window
      + Seq.AQSlice(1).ReadGradTimeDelayOffset ...  % delay of readout gradient
      - Seq.AQSlice(1).SpoilEddyCurrentTime ...
      - Seq.AQSlice(1).SpoilTimeOffset ...
      - max(Seq.Slice(1).Pulse.MaxLength*0.5, max(HW.Grad(Seq.AQSlice(1).iDevice).TimeDelay(1:3))) ...
        * (Seq.AQSlice(1).thickness>=1000) ... % no Slice Grad and no SliceDephase
      - SliceDephaseStart ...
        * (Seq.AQSlice(1).thickness<1000)*(Seq.AQSlice(1).SliceDephase~=0) ... % slice grad and SliceDephase
      - (-Seq.Slice(1).CenterOfPulse + Seq.Slice(1).GradLength/2 ...
         + max(HW.Grad(Seq.AQSlice(1).iDevice).TimeDelay(1:3)) ) ...
        * (Seq.AQSlice(1).thickness<1000)*(Seq.AQSlice(1).SliceDephase==0) ... % slice grad and no SliceDephase
      - Seq.CLTime(1) - 10e-6;

    % reduce timing by duration of read dephase and rephase pulse (for tpre) or
    % duration of spoiler (for tpost)
    tAQmax = min(...
      (tpre - ...
        (HW.Grad(Seq.AQSlice(1).iDevice).tEC + HW.Grad(Seq.AQSlice(1).iDevice).tRamp) - ...
        (HW.Grad(Seq.AQSlice(1).iDevice).tEC + HW.Grad(Seq.AQSlice(1).iDevice).tRamp/2) ...
        * Seq.AQSlice(1).DephaseLengthFactor)*2 / (1+Seq.AQSlice(1).DephaseLengthFactor), ...
      (tpost - ...
        (HW.Grad(Seq.AQSlice(1).iDevice).tEC + HW.Grad(Seq.AQSlice(1).iDevice).tRamp) - ...
        (HW.Grad(Seq.AQSlice(1).iDevice).tEC + HW.Grad(Seq.AQSlice(1).iDevice).tRamp/2) ...
        * Seq.AQSlice(1).SpoilLengthFactor)*2 / (1+Seq.AQSlice(1).SpoilLengthFactor));

    if numImages*Seq.AQSlice(1).EPIFactor > 1
      % Check timing between read rephase (optional) and next dephase (optional)
      if (~Seq.AQSlice(1).EPIGradientPhase || ~Seq.AQSlice(1).EPIGradientRead)
        % account for both dephase and rephase pulse
        % Because of the symmetry, it is enough to look at half the echo time
        tInter = Seq.AQSlice(1).EPIEchoSpacing/2 / (Seq.AQSlice(1).oddEvenEchoes+1) ...
          - Seq.CLTime - 10e-6;

        % reduce timing by duration of read dephase pulse (and by symmetry for the
        % read rephase pulse as well)
        tAQmax = min(tAQmax, ...
          (tInter - ...
            (HW.Grad(Seq.AQSlice(1).iDevice).tEC + HW.Grad(Seq.AQSlice(1).iDevice).tRamp ...
             + Seq.AQSlice(1).DephaseTimeOffset) - ...
            (HW.Grad(Seq.AQSlice(1).iDevice).tEC + HW.Grad(Seq.AQSlice(1).iDevice).tRamp/2) ...
            * Seq.AQSlice(1).DephaseLengthFactor)*2 / (1+Seq.AQSlice(1).DephaseLengthFactor) );
      elseif Seq.AQSlice(1).EPIFactor > 1
        % rephase coincides with dephase of following tRep
        % So don't include the rephase pulse in the timing calculation
        tInter = Seq.AQSlice(1).EPIEchoSpacing / (Seq.AQSlice(1).oddEvenEchoes+1) ...
          - Seq.CLTime - 10e-6;

        % reduce timing by duration of read dephase pulse only
        tAQmax = min(tAQmax, ...
          (tInter - Seq.AQSlice(1).DephaseTimeOffset + ...
           HW.Grad(Seq.AQSlice(1).iDevice).tRamp*Seq.AQSlice(1).DephaseLengthFactor/2) / ...
          (1+Seq.AQSlice(1).DephaseLengthFactor/2) ...
          - 2*(HW.Grad(Seq.AQSlice(1).iDevice).tEC + HW.Grad(Seq.AQSlice(1).iDevice).tRamp + Seq.CLTime));
      else
        % Rephase and dephase in read and phase direction cancel exactly.
        % So use maximum AQ duration under readout gradient.
        tInter = min([Seq.AQSlice(1).EPIEchoSpacing / (Seq.AQSlice(1).oddEvenEchoes+1); diff(Seq.tEcho(:))]) ...
          - 2 * (HW.Grad(Seq.AQSlice(1).iDevice).tEC + HW.Grad(Seq.AQSlice(1).iDevice).tRamp) ...
          - Seq.CLTime - 10e-6;
        tAQmax = min(tAQmax, tInter);
      end
      tAQmax = min(tAQmax(:));
    end

    Seq.AQSlice(1).HzPerPixMin = 1./tAQmax;  % bandwidth per pixel in Hz (1/HzPerPixMin = duration of AQ)
  end
  clear CorrectPhaseDeadTimeTX2RX

  Seq.AQSlice(1).AcquisitionTime = 1/Seq.AQSlice(1).HzPerPixMin;
  autoReadOS = false;
  if isemptyfield(Seq.AQSlice(1), 'ReadOS')
    Seq.AQSlice(1).ReadOS = max(2, ...
      ceil(HW.RX(Seq.AQSlice(1).iDevice).fSample / ...
        min(8000, HW.RX(Seq.AQSlice(1).iDevice).CIC_Decimation_Max) / ...
        (Seq.AQSlice(1).HzPerPixMin * Seq.AQSlice(1).nRead)));
    autoReadOS = true;
  end
  if isemptyfield(Seq.AQSlice(1), 'SamplingFactor')
    % integer sampling factor that is used to (down-)sample the received signal
    % (before image reconstruction)
    minReadOS = ...
      ceil(HW.RX(Seq.AQSlice(1).iDevice).fSample / ...
        min(8000, HW.RX(Seq.AQSlice(1).iDevice).CIC_Decimation_Max) / ...
        (Seq.AQSlice(1).HzPerPixMin*Seq.AQSlice(1).nRead)/Seq.AQSlice(1).ReadOS)*Seq.AQSlice(1).ReadOS;
    if minReadOS > Seq.AQSlice(1).ReadOS
      Seq.AQSlice(1).SamplingFactor = ceil(minReadOS/Seq.AQSlice(1).ReadOS);
    else
      Seq.AQSlice(1).SamplingFactor = 1;
    end
  end
  if autoReadOS
    Seq.AQSlice(1).ReadOS = ceil(Seq.AQSlice(1).ReadOS / Seq.AQSlice(1).SamplingFactor);
  end
  if isemptyfield(Seq.AQSlice(1), 'ReadOSUsedForImage')
    % Optionally, reduce read over-sampling factor for CSI image reconstruction
    Seq.AQSlice(1).ReadOSUsedForImage = min(Seq.AQSlice(1).ReadOS, 17);
  end
  Seq.tEcho = Seq.tEcho + Seq.AQSlice(1).tEchoOffset;
  Seq.AQSlice(1).tEcho = Seq.tEcho;

  % Readout at echoes
  Seq.Read(1).HzPerPixelMin = Seq.AQSlice(1).HzPerPixMin;
  Seq.Read(1).UseCoordinate = Seq.AQSlice(1).ReadCoordinate;
  if isZTE
    % For ZTE, we need an "approximate" read out block in a first step.
    Seq.Read(1).CenterOfReadout = 1/Seq.AQSlice(1).HzPerPixMin/2;
    Seq.Read(1).StartWithKSpaceCenter = 1;
    Seq.Read(1).StartWithKSpaceCenterNoSampleShift = 1;

    % resolution in "phase(1)" direction is the reference for encoding
    Seq.Read(1).sizeRead = Seq.AQSlice(1).sizePhase(1);
    Seq.Read(1).Resolution = Seq.Read(1).sizeRead / (Seq.AQSlice(1).nPhase(1)/2);  % only half of k-space

    % sanity check of input parameters
    usedDim = Seq.AQSlice(1).nPhase > 1;
    if any(diff(usedDim) > 0)
      % FIXME: We don't allow something like nPhase(2)=1 and nPhase(3)=16. This
      %        should catch this. How to make the error message clearer?
      error('PD:sequence_Flash:invalidPhaseDim', ...
        ['For ZTE measurements, nPhase must be used from start to end dimension. ', ...
        'Singleton dimensions must be trailing.']);
    end

    if sum(usedDim) > 2
      Seq.AQSlice(1).VoxelVolume = prod(Seq.AQSlice(1).sizePhase./Seq.AQSlice(1).nPhase);
    elseif sum(usedDim) > 1
      Seq.AQSlice(1).VoxelVolume = prod(Seq.AQSlice(1).sizePhase(1:2)./Seq.AQSlice(1).nPhase(1:2)) .* ...
        HW.RX(Seq.AQSlice(1).iDevice).EffectiveCoilLength;
    else
      Seq.AQSlice(1).VoxelVolume = 1;
    end

    if ~isfinite(Seq.AQSlice(1).VoxelVolume)
      Seq.AQSlice(1).VoxelVolume = 1;
    end

    % Construct a grid on the surface of an N-dimensional cuboid from phase
    % settings and direct the readouts to each of these grid points.
    % end k-space coordinates in "phase" directions
    % If sizePhase is Inf (no encoding in that direction), use an arbitrary
    % large size to be able to aim for a finite end point of the "ray" in
    % k-space.
    maxSizePhase = 1e6;
    nSteps = Seq.AQSlice(1).nPhase.*Seq.AQSlice(1).PhaseOS;
    kmin = -floor(nSteps/2)./(min(Seq.AQSlice(1).sizePhase, maxSizePhase).*Seq.AQSlice(1).PhaseOS);
    kmax = kmin + (nSteps-1)./(min(Seq.AQSlice(1).sizePhase, maxSizePhase).*Seq.AQSlice(1).PhaseOS);
    for iPhase = 3:-1:1
      ksteps{iPhase} = get_FFTGrid(Seq.AQSlice(1).nPhase(iPhase)./min(Seq.AQSlice(1).sizePhase(iPhase), maxSizePhase), nSteps(iPhase));
    end

    endpoints_k1 = zeros(3, []);
    % end planes in "phase(1)" directions
    [k2, k3] = ndgrid(ksteps{2}, ksteps{3});
    endpoints_k1 = [endpoints_k1, ...
      [repmat(kmin(1), 1, numel(k2(2:end,:))); ...
      reshape(k2(2:end,:), 1, []); ...
      reshape(k3(2:end,:), 1, [])]];  %#ok<AGROW>
    endpoints_k1 = [endpoints_k1, ...
      [repmat(kmax(1), 1, numel(k2(1:end-1,:))); ...
      reshape(k2(1:end-1,:), 1, []); ...
      reshape(k3(1:end-1,:), 1, [])]];  %#ok<AGROW>

    % end planes in "phase(2)" directions
    endpoints_k2 = zeros(3, []);
    [k1, k3] = ndgrid(ksteps{1}, ksteps{3});
    endpoints_k2 = [endpoints_k2, ...
      [reshape(k1(1:end-1,:), 1, []); ...
      repmat(kmin(2), 1, numel(k1(1:end-1,:))); ...
      reshape(k3(1:end-1,:), 1, [])]];  %#ok<AGROW>
    endpoints_k2 = [endpoints_k2, ...
      [reshape(k1(2:end,:), 1, []); ...
      repmat(kmax(2), 1, numel(k1(2:end,:))); ...
      reshape(k3(2:end,:), 1, [])]];  %#ok<AGROW>

    % end planes in "phase(3)" directions
    endpoints_k3 = zeros(3, []);
    if Seq.AQSlice(1).nPhase(3) > 1
      [k1, k2] = ndgrid(ksteps{1}, ksteps{2});
      endpoints_k3 = [endpoints_k3, ...
        [reshape(k1(2:end-1,2:end-1), 1, []); ...
        reshape(k2(2:end-1,2:end-1), 1, []); ...
        repmat(kmax(3), 1, numel(k1(2:end-1,2:end-1)))]];  %#ok<AGROW>
      endpoints_k3 = [endpoints_k3, ...
        [reshape(k1(2:end-1,2:end-1), 1, []); ...
        reshape(k2(2:end-1,2:end-1), 1, []); ...
        repmat(kmin(3), 1, numel(k1(2:end-1,2:end-1)))]];  %#ok<AGROW>
    end

    if Seq.plotZTEEndPoints
      hf = figure(1401); clf(hf);
      hax = axes(hf);  %#ok<LAXES>
      scatter3(hax, endpoints_k1(1,:), endpoints_k1(2,:), endpoints_k1(3,:), 'o');
      hold(hax, 'all');
      scatter3(hax, endpoints_k2(1,:), endpoints_k2(2,:), endpoints_k2(3,:), 'x');
      scatter3(hax, endpoints_k3(1,:), endpoints_k3(2,:), endpoints_k3(3,:), '+');
      axis(hax, 'equal');
      if Seq.AQSlice(1).nPhase(3) <= 1
        % set reasonable limits for z-axis
        set(hax, 'ZLim', [-1, 1].*mean([diff(get(hax, 'XLim')), diff(get(hax, 'YLim'))])/2)
      end
      xlabel(hax, 'k1 in 1/m'); ylabel(hax, 'k2 in 1/m'); zlabel(hax, 'k3 in 1/m');
      title('k-space aim points on surface of cuboid');
    end

    % normalize end points
    norm_kmax = reshape(kmax/kmax(1), 3, 1);
    endpoints_k1_norm = zeros(size(endpoints_k1));
    endpoints_k1_norm(usedDim,:) = bsxfun(@rdivide, endpoints_k1(usedDim,:), ...
      sqrt(sum(bsxfun(@rdivide, endpoints_k1(usedDim,:), norm_kmax(usedDim,:)).^2, 1)));
    endpoints_k2_norm = zeros(size(endpoints_k2));
    endpoints_k2_norm(usedDim,:) = bsxfun(@rdivide, endpoints_k2(usedDim,:), ...
      sqrt(sum(bsxfun(@rdivide, endpoints_k2(usedDim,:), norm_kmax(usedDim,:)).^2, 1)));
    endpoints_k3_norm = zeros(size(endpoints_k3));
    endpoints_k3_norm(usedDim,:) = bsxfun(@rdivide, endpoints_k3(usedDim,:), ...
      sqrt(sum(bsxfun(@rdivide, endpoints_k3(usedDim,:), norm_kmax(usedDim,:)).^2, 1)));

    if Seq.plotZTEEndPoints
      hf = figure(1402); clf(hf);
      hax = axes(hf);  %#ok<LAXES>
      scatter3(hax, endpoints_k1_norm(1,:), endpoints_k1_norm(2,:), endpoints_k1_norm(3,:), 'o');
      hold(hax, 'all');
      scatter3(hax, endpoints_k2_norm(1,:), endpoints_k2_norm(2,:), endpoints_k2_norm(3,:), 'x');
      scatter3(hax, endpoints_k3_norm(1,:), endpoints_k3_norm(2,:), endpoints_k3_norm(3,:), '+');
      axis(hax, 'equal');
      xlabel(hax, 'k1'); ylabel(hax, 'k2'); zlabel(hax, 'k3');
      title({'k-space end points on surface of ellipsoid', 'normalized to max(k1)'});
    end


    % rays might be "denser" in different encoding directions
    if sum(Seq.AQSlice(1).nPhase>1) > 2
      % 3d
      kspace_step = min(Seq.AQSlice(1).sizePhase, maxSizePhase).*Seq.AQSlice(1).PhaseOS;
      resolution_factor = 1./[prod(kspace_step([2,3])), prod(kspace_step([1,3])), prod(kspace_step([1,2]))];
      resolution_factor = resolution_factor./resolution_factor(1);
    elseif sum(Seq.AQSlice(1).nPhase>1) > 1
      resolution_factor = min(Seq.AQSlice(1).sizePhase, maxSizePhase).*Seq.AQSlice(1).PhaseOS;
      resolution_factor = resolution_factor./resolution_factor(1);
    else
      resolution_factor = 1;
    end

    % correction for "denser" rays towards the diagonals
    % assumption: points are equidistant on the cuboid
    % FIXME: Is this is correct if the encoding is not the same in k1, k2, k3?

    Seq.AQSlice(1).ReadDensityCorrection = ...
      [resolution_factor(1) ./ (1 + ...
      sqrt(1 - (endpoints_k1_norm(1,:) ./ sqrt(sum(endpoints_k1_norm.^2, 1))).^2)) .^ (sum(Seq.AQSlice(1).nPhase>1)-1), ...
      resolution_factor(2) ./ (1 + ...
      sqrt(1 - (endpoints_k2_norm(2,:) ./ sqrt(sum(endpoints_k2_norm.^2, 1))).^2)) .^ (sum(Seq.AQSlice(1).nPhase>1)-1), ...
      resolution_factor(3) ./ (1 + ...
      sqrt(1 - (endpoints_k3_norm(3,:) ./ sqrt(sum(endpoints_k3_norm.^2, 1))).^2)) .^ (sum(Seq.AQSlice(1).nPhase>1)-1)];

    Seq.Read(1).kLineDir = [endpoints_k1_norm, endpoints_k2_norm, endpoints_k3_norm];

    % There is no read dephase for ZTE measurements.
    % Disable any checks for this gradient block.
    Seq.Read(1).GradDephaseLengthFactor = NaN;
  else
    Seq.Read(1).CenterOfReadout = Seq.AQSlice(1).tEchoOffset;
    % match the phase of the corresponding excitation pulse
    rfPulseImageIdx = (Seq.SteadyState_PreShots+1):(numel(Seq.Slice(1).Pulse.Phase)-Seq.SteadyState_PostShots);
    if Seq.CorrectPhase && Seq.CorrectPhaseSeparate
      % select pulses that are used for frequency tracking
      rfPulseFrequencyIdx = rfPulseImageIdx((Seq.CorrectPhaseBlockSize+1):(Seq.CorrectPhaseBlockSize+1):end);
      % remove pulses that are used for frequency tracking
      rfPulseImageIdx((Seq.CorrectPhaseBlockSize+1):(Seq.CorrectPhaseBlockSize+1):end) = [];
    end
    Seq.Read(1).Phase = repmat(Seq.Slice(1).Pulse.Phase(rfPulseImageIdx), ...
      numImages*Seq.AQSlice(1).EPIFactor, 1);
    Seq.Read(1).PhaseIncrement = Seq.Slice(1).Pulse.PhaseIncrement;

    Seq.Read(1).GradDephaseLengthFactor = Seq.AQSlice(1).DephaseLengthFactor;
  end
  Seq.Read(1).UseAtRepetitionTime = Seq.kLinestReps;
  Seq.Read(1).GradTimeDelay = HW.Grad(Seq.AQSlice(1).iDevice).ReadTimeDelay;
  Seq.Read(1).distance = Seq.AQSlice(1).Center2OriginImage(3);
  Seq.Read(1).GradRephaseLengthFactor = Seq.AQSlice(1).SpoilLengthFactor;
  if Seq.AQSlice(1).EPIReadBipolar
    % read direction polarity in order of "UseAtRepetitionTime"
    ReadDirtRep = 1 - 2*mod(bsxfun(@minus, Seq.kLinestReps, Seq.kLinestReps(1,:)), 2);
    % read direction polarity in order of k-lines per image
    Seq.AQSlice(1).ReadDir = ReadDirtRep(kLinesDim(Seq.kLinesImages));
  else
    ReadDirtRep = 1;
    Seq.AQSlice(1).ReadDir = 1;
  end
  Seq.Read(1).GradSign = ReadDirtRep(:).' * ...
    Seq.AQSlice(1).ReadGradSign(mod(Loop-1, numel(Seq.AQSlice(1).ReadGradSign)) + 1);
  Seq.Read(1).GradRephaseSign = -Seq.AQSlice(1).ReadGradRephaseSign .* ...
    Seq.AQSlice(1).ReadGradSign(mod(Loop-1, numel(Seq.AQSlice(1).ReadGradSign)) + 1);
  Seq.Read(1).GradDephaseSign = ...
    -Seq.AQSlice(1).ReadGradSign(mod(Loop-1, numel(Seq.AQSlice(1).ReadGradSign)) + 1);
  Seq.Read(1).GradTimeIntegralOffset = Seq.AQSlice(1).ReadGradTimeIntegralOffset;

  Seq = get_ReadParameter(Seq, HW);  % generate readout timings

  % multiplicate the length of the spoil and rephase to achieve shorter tEcho with reduced gradient amplitudes
  if Seq.AQSlice(1).SpoilTimeOffset
    Seq.Read(1).CenterOfRephase = Seq.Read(1).CenterOfRephase + Seq.AQSlice(1).SpoilTimeOffset;
    Seq = get_ReadParameter(Seq, HW);  % generate readout timings
  end

  if ~isZTE && Seq.AQSlice(1).sizePhaseSpoil(3) == 0  % set spoiler amp equal to read gradient amplitude
    if Seq.Read(1).GradTimeIntegral == 0  % no readout gradient
      Seq.AQSlice(1).SpoilFactor(3) = 0;
    else
      Seq.AQSlice(1).SpoilFactor(3) = Seq.Read(1).GradTimeIntegral / Seq.Read(1).GradTimeIntegralAQ * ...
        (1 + (Seq.AQSlice(1).SpoilLengthFactor-1)/2);
    end
    % Spoil read direction
    Seq.AQSlice(1).sizePhaseSpoil(3) = Seq.AQSlice(1).sizeRead / Seq.AQSlice(1).nRead ./ ...
      Seq.AQSlice(1).SpoilFactor(3)./2 .* Seq.Read(1).GradSign(end);  % FIXME: Is "end" always ok?
  end

  if Seq.AQSlice(1).DephaseTimeOffset
    Seq.Read(1).CenterOfDephase = Seq.Read(1).CenterOfDephase - Seq.AQSlice(1).DephaseTimeOffset;
    Seq = get_ReadParameter(Seq, HW);  % generate readout timings
  end

  if isZTE
    % slice excitation by Seq.FlipAngle
    Seq.Slice(1).Pulse.Function = Seq.AQSlice(1).excitationPulse;
    Seq.Slice(1).Pulse.FlipAngleComposite = Seq.AQSlice(1).excitationFlipAngleComposite;
    Seq.Slice(1).Pulse.FlipPhaseComposite = Seq.AQSlice(1).excitationFlipPhaseComposite;
    Seq.Slice(1).Pulse.FlipAngle = Seq.FlipAngle;
    % get slice pulse matching readout gradient amplitude
    Seq.Slice(1).GradAmp = sqrt(sum(Seq.Read(1).GradAmpNorm.^2, 1));
    Seq.Slice(1).GradAmp = [repmat(Seq.Slice(1).GradAmp(1), 1, Seq.SteadyState_PreShots), ...
      Seq.Slice(1).GradAmp, repmat(Seq.Slice(1).GradAmp(1), 1, nCorrectPhaseSeparate), ...
      repmat(Seq.Slice(1).GradAmp(end), 1, Seq.SteadyState_PostShots)];
    nPulses = Seq.SteadyState_PreShots + Seq.SteadyState_PostShots + Seq.AQSlice(1).EPISegments + nCorrectPhaseSeparate;
    Seq.Slice(1).Pulse.Phase = ...
      360*rem(cumsum(((1:nPulses) - Seq.SteadyState_PreShots)) ./ nPulses*Seq.SweepPhase, 1) ...
      + (0.5*Seq.SweepPhaseQuadratically *((1:nPulses).^2+(1:nPulses)+2));
    Seq.Slice(1).Pulse.PhaseIncrement = Seq.PhaseIncrement+Seq.PhaseIncrementLoop*(Loop-1);
    Seq.Slice(1).CenterOfPulse = 0;
    Seq.Slice(1).UseCoordinate = Seq.AQSlice(1).SliceCoordinate;
    Seq.Slice(1).GradTimeDelay = HW.Grad(Seq.AQSlice(1).iDevice).SliceTimeDelay;
    Seq.Slice(1).UseAtRepetitionTime = ((1:nPulses).'-1) * ...
      (Seq.AQSlice(1).EPIFactor * Seq.AQSlice(1).nImages * (Seq.AQSlice(1).oddEvenEchoes+1) + 1) + 1;
    Seq.Slice(1).Thickness = Seq.AQSlice(1).thickness;
    Seq.Slice(1).distance = Seq.AQSlice(1).Center2OriginImage(1);  % FIXME: Does that work for ZTE? Probably not...
    Seq.Slice(1).MaxGradAmp = Inf;  % FIXME: Is this ok?
    Seq.Slice(1).GradTimeIntegralRephaseOffset = Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset;
    Seq.Slice(1).GradSign = Seq.AQSlice(1).SliceGradSign(mod(Loop-1,numel(Seq.AQSlice(1).SliceGradSign))+1);
    Seq.Slice(1).GradDephaseSign = -Seq.AQSlice(1).SliceGradSign(mod(Loop-1,numel(Seq.AQSlice(1).SliceGradSign))+1);
    Seq.Slice(1).GradRephaseSign = -Seq.AQSlice(1).SliceGradSign(mod(Loop-1,numel(Seq.AQSlice(1).SliceGradSign))+1);
    % For ZTE, the slice gradient blocks won't be used. Set durations to Inf to
    % avoid errors with gradient ramp checks.
    % FIXME: Would be nice to have a SliceDephaseLengthFactor and
    % SliceRephaseLengthFactor
    Seq.Slice(1).GradDephaseLength = Inf;
    Seq.Slice(1).GradRephaseLength = Inf;

    % generate slice parameters
    Seq = get_SliceParameter(Seq, HW);

    % adjust read-out
    % FIXME: Is it possible to extend the readout gradient only to the front?
    % FIXME: Make the deadtime after the excitation pulse adjustable.
    % Over-estimate minimum dead time between TX pulse and acquisition start
    % using the HzPerPixel of the long acquisition window starting at the pulse
    % center.
    deadTimeTX2RX = get_DeadTimeTX2RX(HW, Seq.Read(1).HzPerPixel*Seq.Read(1).nRead*Seq.Read(1).ReadOS, Seq.Read(1).iDevice);
    Seq.Read(1).tEC = HW.Grad(Seq.Read(1).iDevice).tEC + Seq.Slice(1).Pulse.MaxLength + deadTimeTX2RX;
    Seq.Read(1).CenterOfReadout = Seq.Read(1).CenterOfReadout + (Seq.Slice(1).Pulse.MaxLength/2 + deadTimeTX2RX)/2;
    Seq.Read(1).HzPerPixelMin = 1/(Seq.Read(1).AcquisitionTime - (Seq.Slice(1).Pulse.MaxLength/2 + deadTimeTX2RX));
    Seq.Read(1).Resolution = Seq.Read(1).Resolution * Seq.Read(1).HzPerPixelMin / Seq.Read(1).HzPerPixel;
    Seq.Read(1).GradLength = [];
    % match the phase of the corresponding excitation pulse
    rfPulseImageIdx = (Seq.SteadyState_PreShots+1):(numel(Seq.Slice(1).Pulse.Phase)-Seq.SteadyState_PostShots);
    if Seq.CorrectPhase && Seq.CorrectPhaseSeparate
      % select pulses that are used for frequency tracking
      rfPulseFrequencyIdx = rfPulseImageIdx((Seq.CorrectPhaseBlockSize+1):(Seq.CorrectPhaseBlockSize+1):end);
      % remove pulses that are used for frequency tracking
      rfPulseImageIdx((Seq.CorrectPhaseBlockSize+1):(Seq.CorrectPhaseBlockSize+1):end) = [];
    end
    Seq.Read(1).Phase = repmat(Seq.Slice(1).Pulse.Phase(rfPulseImageIdx), ...
      numImages*Seq.AQSlice(1).EPIFactor, 1);
    Seq.Read(1).PhaseIncrement = Seq.Slice(1).Pulse.PhaseIncrement;

    Seq = get_ReadParameter(Seq, HW);  % generate readout timings

    % Now that we have the final timings, calculate the k-space trajectory,
    % i.e., the start and end points of the rays in k-space
    % in m^(-1) if encoded; or in seconds if not encoded
    encodingStrength = Seq.Read(1).Gamma/(2*pi) * Seq.Read(1).GradAmpNorm;
    nonEncodedDirs = isinf(Seq.AQSlice(1).sizePhase) & (Seq.AQSlice(1).nPhase > 1);
    encodingStrength(nonEncodedDirs,:) = Seq.Read(1).kLineDir(nonEncodedDirs,:);
    Seq.AQSlice(1).kLineStart = encodingStrength .* ...
      (Seq.Read(1).CenterOfReadout - Seq.Read(1).AcquisitionTime/2 + 1./(2*Seq.Read(1).fSample));
    Seq.AQSlice(1).kLineEnd = encodingStrength .* ...
      (Seq.Read(1).CenterOfReadout + Seq.Read(1).AcquisitionTime/2 - 1./(2*Seq.Read(1).fSample));

  else

    if ~Seq.AQSlice(1).SliceRephaseImmediately && ...
        (Seq.AQSlice(1).thickness < 1000) && (Seq.CorrectPhase == 0 || Seq.CorrectPhaseSeparate)
      % Slice rephase aligned with read dephase of first echo
      Seq.Slice(1).GradRephaseLength = Seq.Read(1).GradDephaseLength;
      Seq.Slice(1).CenterOfRephase = Seq.Read(1).CenterOfDephase - Seq.Read(1).GradTimeDelayOffset(1);
      Seq.Slice(1).UseAtRepetitionTimeRephase = Seq.Slice(1).UseAtRepetitionTimeRephase + 1;
      Seq = get_SliceParameter(Seq, HW);
    end

    % check if slice rephase interferes with readout window
    if Seq.AQSlice(1).thickness < 1e3 && ...
        Seq.Slice(1).CenterOfRephase + Seq.Slice(1).GradRephaseLength/2 > ...
        Seq.tRep(Seq.Slice(1).UseAtRepetitionTime(1)) + Seq.Read(1).CenterOfReadout - Seq.Read(1).GradLength/2  + Seq.Read(1).tRamp + 1e-12
      error('PD:sequence_Flash:SliceRephaseInReadOut', ...
        'Slice rephase gradient too close to read out.');
    end
  end

  Seq.AQSlice(1).ReadPrePostIdx = [];
  if Seq.SteadyState_PreShots || Seq.SteadyState_PostShots
    Seq.AQSlice(1).ReadPrePostIdx = 2;
    Seq.Read(2) = Seq.Read(1);
    Seq.Read(2).Phase = 0;
    Seq.Read(2).PhaseIncrement = 0;
    Seq.Read(2).kLineDir = Seq.Read(1).kLineDir(:,1);
    if Seq.AQSlice(1).EPIReadBipolar
      readGradBipolar = 2*mod(repmat(1:numImages*Seq.AQSlice(1).EPIFactor, ...
        Seq.SteadyState_PreShots+Seq.SteadyState_PostShots, 1), 2).' - 1;
    else
      readGradBipolar = 1;
    end
    Seq.Read(2).GradSign = readGradBipolar(:).' * Seq.AQSlice(1).ReadGradSign(mod(Loop-1,numel(Seq.AQSlice(1).ReadGradSign))+1);
    Seq.Read(2).GradRephaseSign = -Seq.AQSlice(1).ReadGradSign(mod(Loop-1,numel(Seq.AQSlice(1).ReadGradSign))+1).*Seq.AQSlice(1).ReadGradRephaseSign;
    Seq.Read(2).GradDephaseSign = -Seq.AQSlice(1).ReadGradSign(mod(Loop-1,numel(Seq.AQSlice(1).ReadGradSign))+1);
    Seq.Read(2).UseAtRepetitionTime = ...
      bsxfun(@plus, [(1:Seq.SteadyState_PreShots)-1, ...
                     (nPulses - (Seq.SteadyState_PostShots:-1:1))].' * (numImages*Seq.AQSlice(1).EPIFactor+1), ...
                    1:numImages*Seq.AQSlice(1).EPIFactor).' + 1;
    Seq.Read(2).UseAtRepetitionTimeDephase = Seq.Read(2).UseAtRepetitionTime;
    Seq.Read(2).UseAtRepetitionTimeRephase = Seq.Read(2).UseAtRepetitionTime;

    Seq = get_ReadParameter(Seq, HW);  % generate readout timings
  end

  Seq.AQSlice(1).AcquisitionTime = Seq.Read(1).AcquisitionTime;
  Seq.AQSlice(1).AcquisitionFrequency = Seq.Read(1).AcquisitionFrequency;


  if ~isZTE
    % no phase encoding or spoilers for ZTE
    % FIXME: Add support for spoilers to ZTE measurements.

    % phase encoding of last gradient echo in EPI segment in phase(1) direction
    % aligned with read dephase and rephase
    Seq.Phase(1).sizePhase = Seq.AQSlice(1).sizePhase(1);
    Seq.Phase(1).nPhase = Seq.AQSlice(1).nPhase(1);
    Seq.Phase(1).PhaseOS = Seq.AQSlice(1).PhaseOS(1);
    Seq.Phase(1).StepIncrement = Seq.AQSlice(1).nImages * (Seq.AQSlice(1).oddEvenEchoes+1);
    Seq.Phase(1).CenterOfDephase = Seq.Read(1).CenterOfDephase;
    Seq.Phase(1).CenterOfRephase = Seq.Read(1).CenterOfRephase;
    Seq.Phase(1).GradDephaseLength = Seq.Read(1).GradDephaseLength;
    Seq.Phase(1).GradRephaseLength = Seq.Read(1).GradRephaseLength;
    Seq.Phase(1).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(1);
    Seq.Phase(1).GradTimeDelay = HW.Grad(Seq.AQSlice(1).iDevice).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(end);
    Seq.Phase(1).UseAtRepetitionTime = Seq.AQSlice(1).UsetRep;
    Seq.Phase(1).usedkLines = sort(Seq.AQSlice(1).kLineOrder);
    Seq.Phase(1).numKLines = numKLines;
    Seq.Phase(1).kLineIncrement = Seq.Phase(1).StepIncrement;
    Seq.Phase(1).distance = Seq.AQSlice(1).Center2OriginImage(1);


    % phase encoding of last gradient echo in EPI segment in phase(2) direction
    % aligned with read dephase and rephase
    Seq.Phase(2).sizePhase = Seq.AQSlice(1).sizePhase(2);
    Seq.Phase(2).nPhase = Seq.AQSlice(1).nPhase(2);
    Seq.Phase(2).PhaseOS = Seq.AQSlice(1).PhaseOS(2);
    Seq.Phase(2).StepIncrement = Seq.AQSlice(1).nPhase(1) * Seq.AQSlice(1).PhaseOS(1) * ...
      Seq.Phase(1).StepIncrement;
    Seq.Phase(2).CenterOfDephase = Seq.Read(1).CenterOfDephase;
    Seq.Phase(2).CenterOfRephase = Seq.Read(1).CenterOfRephase;
    Seq.Phase(2).GradDephaseLength = Seq.Read(1).GradDephaseLength;
    Seq.Phase(2).GradRephaseLength = Seq.Read(1).GradRephaseLength;
    Seq.Phase(2).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(2);
    Seq.Phase(2).GradTimeDelay = HW.Grad(Seq.AQSlice(1).iDevice).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(end);
    Seq.Phase(2).UseAtRepetitionTime = Seq.Phase(1).UseAtRepetitionTime;
    Seq.Phase(2).usedkLines = Seq.Phase(1).usedkLines;
    Seq.Phase(2).numKLines = Seq.Phase(1).numKLines;
    Seq.Phase(2).kLineIncrement = Seq.Phase(1).kLineIncrement;
    Seq.Phase(2).distance = Seq.AQSlice(1).Center2OriginImage(2);


    % phase encoding of last gradient echo in EPI segment in phase(3) direction
    % aligned with read dephase and rephase
    Seq.Phase(3).sizePhase = Seq.AQSlice(1).sizePhase(3);
    Seq.Phase(3).nPhase = Seq.AQSlice(1).nPhase(3);
    Seq.Phase(3).PhaseOS = Seq.AQSlice(1).PhaseOS(3);
    Seq.Phase(3).StepIncrement = Seq.AQSlice(1).nPhase(2) * Seq.AQSlice(1).PhaseOS(2) * ...
      Seq.Phase(2).StepIncrement;
    Seq.Phase(3).CenterOfDephase = Seq.Read(1).CenterOfDephase;
    Seq.Phase(3).CenterOfRephase = Seq.Read(1).CenterOfRephase;
    Seq.Phase(3).GradDephaseLength = Seq.Read(1).GradDephaseLength;
    Seq.Phase(3).GradRephaseLength = Seq.Read(1).GradRephaseLength;
    Seq.Phase(3).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(3);
    Seq.Phase(3).GradTimeDelay = HW.Grad(Seq.AQSlice(1).iDevice).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(end);
    Seq.Phase(3).UseAtRepetitionTime = Seq.Phase(1).UseAtRepetitionTime;
    Seq.Phase(3).usedkLines = Seq.Phase(1).usedkLines;
    Seq.Phase(3).numKLines = Seq.Phase(1).numKLines;
    Seq.Phase(3).kLineIncrement = Seq.Phase(1).kLineIncrement;
    Seq.Phase(3).distance = Seq.AQSlice(1).Center2OriginImage(3);


    % spoiler encoding in phase(1) or slice direction aligned with read rephase at last echo
    Seq.Phase(4).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(1);
    Seq.Phase(4).nPhase = 1;
    Seq.Phase(4).PhaseOS = 2;
    Seq.Phase(4).StepOrder = [1, 1];
    Seq.Phase(4).CenterOfDephase = Seq.Read(1).CenterOfDephase;
    Seq.Phase(4).CenterOfRephase = Seq.Read(1).CenterOfRephase;
    Seq.Phase(4).GradDephaseLength = Seq.Read(1).GradDephaseLength;
    Seq.Phase(4).GradRephaseLength = Seq.Read(1).GradRephaseLength;
    Seq.Phase(4).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(1);
    Seq.Phase(4).GradTimeDelay = HW.Grad(Seq.AQSlice(1).iDevice).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(end);
    if Seq.CorrectPhase > 0 && Seq.CorrectPhaseSeparate
      Seq.Phase(4).UseAtRepetitionTime = ...
        setdiff(Seq.Slice(1).UseAtRepetitionTime(Seq.SteadyState_PreShots+2:end-Seq.SteadyState_PostShots), ...
        rfPulseFrequencytReps + nEPISegment + 1) - 1;
      Seq.Phase(4).UseAtRepetitionTime(end+1) = numel(Seq.tRep)-(Seq.SteadyState_PostShots+1)*(nEPISegment+1);
    else
      Seq.Phase(4).UseAtRepetitionTime = Seq.Slice(1).UseAtRepetitionTime(2:end)-1;
      Seq.Phase(4).UseAtRepetitionTime(end+1) = numel(Seq.tRep);
    end


    % spoiler encoding in phase(2) direction aligned with read rephase at last echo
    Seq.Phase(5).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(2);
    Seq.Phase(5).nPhase = 1;
    Seq.Phase(5).PhaseOS = 2;
    Seq.Phase(5).StepOrder = [1, 1];
    Seq.Phase(5).CenterOfDephase = Seq.Read(1).CenterOfDephase;
    Seq.Phase(5).CenterOfRephase = Seq.Read(1).CenterOfRephase;
    Seq.Phase(5).GradDephaseLength = Seq.Read(1).GradDephaseLength;
    Seq.Phase(5).GradRephaseLength = Seq.Read(1).GradRephaseLength;
    Seq.Phase(5).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(2);
    Seq.Phase(5).GradTimeDelay = HW.Grad(Seq.AQSlice(1).iDevice).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(end);
    Seq.Phase(5).UseAtRepetitionTime = Seq.Phase(4).UseAtRepetitionTime;


    % spoiler encoding in phase(3) or read direction aligned with read rephase at last echo
    Seq.Phase(6).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(3);
    Seq.Phase(6).nPhase = 1;
    Seq.Phase(6).PhaseOS = 2;
    Seq.Phase(6).StepOrder = [1, 1];
    Seq.Phase(6).CenterOfDephase = Seq.Read(1).CenterOfDephase;
    Seq.Phase(6).CenterOfRephase = Seq.Read(1).CenterOfRephase;
    Seq.Phase(6).GradDephaseLength = Seq.Read(1).GradDephaseLength;
    Seq.Phase(6).GradRephaseLength = Seq.Read(1).GradRephaseLength;
    Seq.Phase(6).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(3) * (Seq.AQSlice(1).nPhase(3)~=1) + ...
      Seq.AQSlice(1).ReadCoordinate * (Seq.AQSlice(1).nPhase(3)==1);
    Seq.Phase(6).GradTimeDelay = HW.Grad(Seq.AQSlice(1).iDevice).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(end);
    Seq.Phase(6).UseAtRepetitionTime = Seq.Phase(4).UseAtRepetitionTime;


    if (Seq.AQSlice(1).thickness<1000) && (Seq.CorrectPhase==1 && ~Seq.CorrectPhaseSeparate) && ...
        (Seq.AQSlice(1).SliceRephaseTimeOffset ~= 0)
      Seq.Slice(1).CenterOfRephase = Seq.Slice(1).CenterOfRephase + Seq.AQSlice(1).SliceRephaseTimeOffset;

      % generate slice timings
      Seq = get_SliceParameter(Seq, HW);
    end


    % generate phase parameters
    for t = 1:length(Seq.Phase)
      if ischar(Seq.Phase(t).UseAtRepetitionTime)
        Seq.Phase(t).UseAtRepetitionTime = eval(Seq.Phase(t).UseAtRepetitionTime);
      end
    end

    % apply all phase gradients to the same device
    [Seq.Phase(:).iDevice] = deal(Seq.AQSlice(1).iDevice);

    Seq = get_PhaseParameter(Seq, HW);

    if Seq.AQSlice(1).EPIGradientPhase
      Seq.Phase(1).CenterOfRephase = Seq.Phase(1).CenterOfRephase .* ones(size(Seq.Phase(1).UseAtRepetitionTimeRephase));
      Seq.Phase(2).CenterOfRephase = Seq.Phase(2).CenterOfRephase .* ones(size(Seq.Phase(2).UseAtRepetitionTimeRephase));
      Seq.Phase(3).CenterOfRephase = Seq.Phase(3).CenterOfRephase .* ones(size(Seq.Phase(3).UseAtRepetitionTimeRephase));

      % move rephase gradient (of image encoder) to dephase gradient in next tRep
      Seq.Phase(7:9) = Seq.Phase(1:3);
      Seq.Phase(7).UseAtRepetitionTimeRephase = Seq.Phase(1).UseAtRepetitionTimeRephase + 1;
      Seq.Phase(8).UseAtRepetitionTimeRephase = Seq.Phase(2).UseAtRepetitionTimeRephase + 1;
      Seq.Phase(9).UseAtRepetitionTimeRephase = Seq.Phase(3).UseAtRepetitionTimeRephase + 1;
      Seq.Phase(7).CenterOfRephase(:) = Seq.Phase(1).CenterOfDephase;
      Seq.Phase(8).CenterOfRephase(:) = Seq.Phase(2).CenterOfDephase;
      Seq.Phase(9).CenterOfRephase(:) = Seq.Phase(3).CenterOfDephase;
      Seq.Phase(7).GradRephaseLength(:) = Seq.Phase(1).GradDephaseLength;
      Seq.Phase(8).GradRephaseLength(:) = Seq.Phase(2).GradDephaseLength;
      Seq.Phase(9).GradRephaseLength(:) = Seq.Phase(3).GradDephaseLength;

      % but immediately re-phase at end of EPI segment
      [~, iLastInSegment] = ismember(Seq.kLinestReps(end,:), Seq.Phase(1).UseAtRepetitionTimeRephase);
      Seq.Phase(7).CenterOfRephase(iLastInSegment) = NaN;
      Seq.Phase(8).CenterOfRephase(iLastInSegment) = NaN;
      Seq.Phase(9).CenterOfRephase(iLastInSegment) = NaN;

      % remove re-phase gradients from moved positions
      % FIXME: Optionally remove from all tReps?
      iNotLastInSegment = true(size(Seq.kLinestReps));
      iNotLastInSegment(iLastInSegment) = false;
      Seq.Phase(1).CenterOfRephase(iNotLastInSegment) = NaN;
      Seq.Phase(2).CenterOfRephase(iNotLastInSegment) = NaN;
      Seq.Phase(3).CenterOfRephase(iNotLastInSegment) = NaN;

      % remove duplicate de-phase gradients
      Seq.Phase(7).CenterOfDephase(:) = NaN;
      Seq.Phase(8).CenterOfDephase(:) = NaN;
      Seq.Phase(9).CenterOfDephase(:) = NaN;

      Seq = get_PhaseParameter(Seq, HW);

      for iPhase = 7:9
        for iChannel = 1:numel(Seq.Phase(iPhase).GradDephase)
          Seq.Phase(iPhase).GradRephase(iChannel).Amp = Seq.Phase(iPhase).GradRephase(iChannel).Amp(:,1:numel(Seq.tRep));
          Seq.Phase(iPhase).GradRephase(iChannel).Time = Seq.Phase(iPhase).GradRephase(iChannel).Time(:,1:numel(Seq.tRep));
        end
      end
    end
  end

  Seq.AQSlice(1).EPIGradienReadIdx = [];
  Seq.AQSlice(1).EPIGradienReadPrePostIdx = [];
  if Seq.AQSlice(1).EPIGradientRead
    % compensate gradients at last tRep of each EPI segment
    % (must be done before the center timings are overridden in the next block)
    Seq.AQSlice(1).EPIGradienReadIdx = 3;
    Seq.Read(3) = Seq.Read(1);
    Seq.Read(3).Phase = 0;
    Seq.Read(3).PhaseIncrement = 0;
    Seq.Read(3).UseAtRepetitionTime = Seq.tRepEPISegment(end,:);
    Seq.Read(3).GradSign = ReadDirtRep(end,:);
    Seq.Read(3).UseAtRepetitionTimeDephase = Seq.Read(3).UseAtRepetitionTime+1;
    Seq.Read(3).UseAtRepetitionTimeRephase = Seq.Read(3).UseAtRepetitionTime;
    Seq.Read(3).CenterOfDephase = Seq.Read(3).CenterOfDephase + Seq.Read(1).GradTimeDelayOffset(end) - Seq.Read(1).GradTimeDelayOffset(1);
    Seq.Read(3).GradDephaseSign = -Seq.Read(1).GradDephaseSign;
    Seq.Read(3).GradTimeIntegralOffset = Seq.Read(3).GradTimeIntegralOffset(end);
    Seq.Read(3).GradTimeDelayOffset = Seq.Read(3).GradTimeDelayOffset(end);

    % move rephase gradient to dephase gradient in next tRep
    Seq.Read(1).UseAtRepetitionTimeRephase = Seq.Read(1).UseAtRepetitionTimeRephase + 1;
    Seq.Read(1).CenterOfRephase = Seq.Read(1).CenterOfDephase + Seq.Read(1).GradTimeDelayOffset(:).' - Seq.Read(1).GradTimeDelayOffset([2:end,1]);
    Seq.Read(1).GradRephaseLength = Seq.Read(1).GradDephaseLength;

    if Seq.SteadyState_PreShots > 0 || Seq.SteadyState_PostShots > 0
      % do the same in pre and post shots
      Seq.AQSlice(1).EPIGradienReadPrePostIdx = 4;
      Seq.Read(4) = Seq.Read(2);
      Seq.Read(4).Phase = 0;
      Seq.Read(4).PhaseIncrement = 0;
      Seq.Read(4).UseAtRepetitionTime = Seq.Read(4).UseAtRepetitionTime(end,:);
      if Seq.AQSlice(1).EPIReadBipolar
        Seq.Read(4).GradSign = readGradBipolar(end,:) * Seq.AQSlice(1).ReadGradSign(mod(Loop-1,numel(Seq.AQSlice(1).ReadGradSign))+1);
      else
        Seq.Read(4).GradSign = Seq.Read(2).GradSign(1);
      end
      Seq.Read(4).UseAtRepetitionTimeDephase = Seq.Read(4).UseAtRepetitionTime+1;
      Seq.Read(4).UseAtRepetitionTimeRephase = Seq.Read(4).UseAtRepetitionTime;
      Seq.Read(4).CenterOfDephase = Seq.Read(4).CenterOfDephase + Seq.Read(2).GradTimeDelayOffset(end) - Seq.Read(2).GradTimeDelayOffset(1);
      Seq.Read(4).GradDephaseSign = -Seq.Read(2).GradDephaseSign;
      Seq.Read(4).GradTimeIntegralOffset = Seq.Read(4).GradTimeIntegralOffset(end);
      Seq.Read(4).GradTimeDelayOffset = Seq.Read(4).GradTimeDelayOffset(end);

      Seq.Read(2).UseAtRepetitionTimeRephase = Seq.Read(2).UseAtRepetitionTimeRephase + 1;
      Seq.Read(2).CenterOfRephase = Seq.Read(2).CenterOfDephase + Seq.Read(2).GradTimeDelayOffset(:).' - Seq.Read(2).GradTimeDelayOffset([2:end,1]);
      Seq.Read(2).GradRephaseLength = Seq.Read(2).GradDephaseLength;
    end

    Seq = get_ReadParameter(Seq, HW);  % generate readout timings
  end

  if ~isZTE
    % Image shift in phase direction
    [~, b] = sort(Seq.kLines);
    Seq.Read(1).Phase(b) = reshape(Seq.Read(1).Phase(b), 1, []) + ...
      reshape(Seq.Phase(1).AQPhaseShift + ...
              Seq.Phase(2).AQPhaseShift + ...
              Seq.Phase(3).AQPhaseShift, 1, [])/pi*180;
    % update readout timings
    Seq = get_ReadParameter(Seq, HW);
  end


  Seq.CorrectPhaseRead = [];
  if Seq.CorrectPhase > 0
    Seq.CorrectPhaseRead = 5;
    Seq.Read(Seq.CorrectPhaseRead).useAQSlice = 2;
    Seq.Read(Seq.CorrectPhaseRead).HzPerPixelMin = Seq.AQSlice(2).HzPerPixMin;
    Seq.Read(Seq.CorrectPhaseRead).CenterOfReadout = Seq.CorrectPhaseAQtOffset ...
      + 0.5/Seq.Read(Seq.CorrectPhaseRead).HzPerPixelMin ...
      - Seq.CorrectPhaseSeparate * Seq.tEcho(1);
    if Seq.CorrectPhaseSeparate
      Seq.Read(Seq.CorrectPhaseRead).UseAtRepetitionTime = ...
        [Seq.rfPulsetReps(1:Seq.SteadyState_PreShots), rfPulseFrequencytReps, ...
         Seq.rfPulsetReps((end-Seq.SteadyState_PostShots+1):end)]+1;
      Seq.Read(Seq.CorrectPhaseRead).Phase = ...
        Seq.Slice(1).Pulse.Phase([1:Seq.SteadyState_PreShots, rfPulseFrequencyIdx, (end-Seq.SteadyState_PostShots+1):end]);
      Seq.Read(Seq.CorrectPhaseRead).sizeRead = Inf;
    else
      % Readout at FID
      Seq.Read(Seq.CorrectPhaseRead).PhaseIncrement = Seq.Read(1).PhaseIncrement;
      Seq.Read(Seq.CorrectPhaseRead).Phase = Seq.Slice(1).Pulse.Phase(:);
      Seq.Read(Seq.CorrectPhaseRead).UseAtRepetitionTime = Seq.AQSlice(2).UsetRep;
    end

    % generate readout timings
    Seq = get_ReadParameter(Seq, HW);

    % Check for collision between phase correction readout and dephase gradient
    % pulse.
    if ~Seq.CorrectPhaseSeparate && ...
        (Seq.Read(Seq.CorrectPhaseRead).CenterOfReadout + Seq.Read(Seq.CorrectPhaseRead).AcquisitionTime/2 > ...
         Seq.tRep(Seq.Read(Seq.CorrectPhaseRead).UseAtRepetitionTime(1)) + Seq.Read(1).CenterOfDephase - Seq.Read(1).GradDephaseLength/2)
      if ~all([abs(Seq.AQSlice(1).sizeRead)>1000, abs(Seq.AQSlice(1).sizePhase)>1000])
        error('tEcho too short for selected AQ length (1/HzPerPixelMin)');
      else
        if Seq.Read(Seq.CorrectPhaseRead).CenterOfReadout + Seq.Read(Seq.CorrectPhaseRead).AcquisitionTime/2 > ...
            Seq.tRep(Seq.Read(Seq.CorrectPhaseRead).UseAtRepetitionTime(1)) + Seq.Read(1).CenterOfReadout - Seq.Read(1).GradLength/2
          error('tEcho too short for selected AQ length (1/HzPerPixelMin)');
        end
      end
    end

    Seq.AQSlice(2).AcquisitionTime = Seq.Read(Seq.CorrectPhaseRead).AcquisitionTime;
    Seq.AQSlice(2).AcquisitionFrequency = Seq.Read(Seq.CorrectPhaseRead).AcquisitionFrequency;
    Seq.AQSlice(2).UsetRep = Seq.Read(Seq.CorrectPhaseRead).UseAtRepetitionTime;
  else
    if Seq.AQSlice(1).SliceRephaseImmediately && Seq.AQSlice(1).thickness < 1000
      % Check collision between end of slice rephase gradient (same slot like
      % dephase gradient) and start of read-out gradient.
      if Seq.Slice(1).CenterOfRephase + Seq.Slice(1).GradRephaseLength/2 > ...
          Seq.tRep(Seq.Slice(1).UseAtRepetitionTimeRephase(1)) + Seq.Read(1).CenterOfReadout - Seq.Read(1).GradLength/2 + Seq.Read(1).tRamp + 1e-9
        error('tEcho too short for selected AQ length (1/HzPerPixelMin)');
      end
    elseif isZTE  % never slice-selective
      % Check collision between end of rf pulse and start of acquisition window.
      % FIXME: Can this ever happen due to wrong user input? If it can't, we can
      % probably skip this kind of check for ZTE measurements.
      if Seq.Slice(1).CenterOfPulse + Seq.Slice(1).Pulse.MaxLength/2 > ...
          Seq.tRep(Seq.Slice(1).UseAtRepetitionTimeRephase(1)) + Seq.Read(1).CenterOfReadout - Seq.Read(1).AcquisitionTime/2 + HW.TX2RXdeadTime
        error('tEcho too short for selected AQ length (1/HzPerPixelMin)');
      end
    else
      % Check collision between end of rf pulse and start of dephase gradient
      % pulse.
      if Seq.Slice(1).CenterOfPulse + Seq.Slice(1).Pulse.MaxLength/2 > ...
          Seq.tRep(Seq.Slice(1).UseAtRepetitionTime(1)) + Seq.Read(1).CenterOfDephase - Seq.Read(1).GradDephaseLength/2 + 1e-9
        error('tEcho too short for selected AQ length (1/HzPerPixelMin)');
      end
    end

  end

  %% User function hook for sequence manipulation
  if ~isempty(Seq.Function_Prepare_Measurement)
    if iscell(Seq.Function_Prepare_Measurement)
      for iPrep = 1:numel(Seq.Function_Prepare_Measurement)
        [HW, Seq, AQ, TX, Grad] = Seq.Function_Prepare_Measurement{iPrep}(HW, Seq, AQ, TX, Grad);
      end
    else
      [HW, Seq, AQ, TX, Grad] = Seq.Function_Prepare_Measurement(HW, Seq, AQ, TX, Grad);
    end
  end


  %% merge all elements into the complete pulse program
  % add AQ
  AQ = add_AQ(AQ, Seq.Read(1).AQ);

  if Seq.CorrectPhase > 0  % readout for frequency tracking
    AQ = add_AQ(Seq.Read(Seq.CorrectPhaseRead).AQ, AQ);
  end
  % add TX pulses
  for t = 1:length(Seq.Slice)
    TX = add_TX(TX, Seq.Slice(t).TX);
  end

  if ~isZTE
    % add slice encoders
    if Seq.AQSlice(1).thickness < 1000
      if Seq.AQSlice(1).SliceDephase ~= 0
        Grad = add_Grad(Grad, Seq.Slice(1).GradDephase);
      end
      Grad = add_Grad(Grad, Seq.Slice(1).Grad);
      Grad = add_Grad(Grad, Seq.Slice(1).GradRephase);
    end
  end

  % add read encoders
  iLastRead = 1;
  Grad = add_Grad(Grad, Seq.Read(1).Grad);
  if ~isZTE
    Grad = add_Grad(Grad, Seq.Read(1).GradDephase);
    if Seq.AQSlice(1).RephaseReadEncoding
      Grad = add_Grad(Grad, Seq.Read(1).GradRephase);
    end
  end
  if ~isempty(Seq.AQSlice(1).ReadPrePostIdx) && ~(Seq.CorrectPhase > 0 && Seq.CorrectPhaseSeparate)
    iLastRead = max(iLastRead, Seq.AQSlice(1).ReadPrePostIdx);
    Grad = add_Grad(Grad, Seq.Read(Seq.AQSlice(1).ReadPrePostIdx).Grad);
    if ~isZTE
      Grad = add_Grad(Grad, Seq.Read(Seq.AQSlice(1).ReadPrePostIdx).GradDephase);
      if Seq.AQSlice(1).RephaseReadEncoding
        Grad = add_Grad(Grad, Seq.Read(Seq.AQSlice(1).ReadPrePostIdx).GradRephase);
      end
    end
  end
  if ~isZTE
    if ~isempty(Seq.AQSlice(1).EPIGradienReadIdx)
      iLastRead = max(iLastRead, Seq.AQSlice(1).EPIGradienReadIdx);
      Grad = add_Grad(Grad, Seq.Read(Seq.AQSlice(1).EPIGradienReadIdx).GradDephase);
      if Seq.AQSlice(1).RephaseReadEncoding
        Grad = add_Grad(Grad, Seq.Read(Seq.AQSlice(1).EPIGradienReadIdx).GradRephase);
      end
      if ~isempty(Seq.AQSlice(1).EPIGradienReadPrePostIdx)
        iLastRead = max(iLastRead, Seq.AQSlice(1).EPIGradienReadPrePostIdx);
        Grad = add_Grad(Grad, Seq.Read(Seq.AQSlice(1).EPIGradienReadPrePostIdx).GradDephase);
        if Seq.AQSlice(1).RephaseReadEncoding
          Grad = add_Grad(Grad, Seq.Read(Seq.AQSlice(1).EPIGradienReadPrePostIdx).GradRephase);
        end
      end
    end

    % add phase encoders
    iLastPhase = 6;
    Grad = add_Grad(Grad, Seq.Phase(1).GradDephase);
    Grad = add_Grad(Grad, Seq.Phase(2).GradDephase);
    Grad = add_Grad(Grad, Seq.Phase(3).GradDephase);
    if Seq.AQSlice(1).RephasePhaseEncoding
      Grad = add_Grad(Grad, Seq.Phase(1).GradRephase);
      Grad = add_Grad(Grad, Seq.Phase(2).GradRephase);
      Grad = add_Grad(Grad, Seq.Phase(3).GradRephase);
      if Seq.AQSlice(1).EPIGradientPhase
        iLastPhase = 9;
        Grad = add_Grad(Grad, Seq.Phase(7).GradRephase);
        Grad = add_Grad(Grad, Seq.Phase(8).GradRephase);
        Grad = add_Grad(Grad, Seq.Phase(9).GradRephase);
      end
    end
    Grad = add_Grad(Grad, Seq.Phase(4).GradRephase);
    Grad = add_Grad(Grad, Seq.Phase(5).GradRephase);
    Grad = add_Grad(Grad, Seq.Phase(6).GradRephase);

    % extra Gradients
    for t = iLastPhase+1:numel(Seq.Phase)
      Grad = add_Grad(Grad, Seq.Phase(t).GradDephase);
      Grad = add_Grad(Grad, Seq.Phase(t).GradRephase);
    end
    for t = iLastRead+(Seq.CorrectPhase > 0)+1:numel(Seq.Read)
      Grad = add_Grad(Grad, Seq.Read(t).Grad);
      Grad = add_Grad(Grad, Seq.Read(t).GradDephase);
      Grad = add_Grad(Grad, Seq.Read(t).GradRephase);
    end
    for t = (1 + (Seq.CorrectPhase>0 && Seq.CorrectPhaseSeparate) + 1):numel(Seq.Slice)
      Grad = add_Grad(Grad, Seq.Slice(t).GradDephase);
      Grad = add_Grad(Grad, Seq.Slice(t).Grad);
      Grad = add_Grad(Grad, Seq.Slice(t).GradRephase);
    end
  end


  %% extent last tRep for clamp coil signal (if it contains AQ window)
  if HW.RX(Seq.AQSlice(1).iDevice).ClampCoil.Enable
    lastAQ = find(~isnan(AQ.Start(:,end)), 1, 'last');
    if ~isempty(lastAQ)
      Seq.tRep(end) = max(Seq.tRep(end), ...
        AQ.Start(lastAQ,end) + AQ.nSamples(lastAQ,end)/AQ.fSample(lastAQ,end) + HW.RX(Seq.AQSlice(1).iDevice).ClampCoil.tPostset + 0.1e-3);
    end
  end


  if Seq.SingletRep > 0
    %% combine tReps of each EPI segment including excitation into singular logical tReps

    tRepsExcitation = Seq.Slice(1).UseAtRepetitionTime(:);

    temptRepStart = cumsum([0,Seq.tRep(1:end-1)]);
    temptRepEnd = cumsum(Seq.tRep);
    temptRepCombined = temptRepEnd([tRepsExcitation(2:end);end+1]-1) ...
      - temptRepStart(tRepsExcitation);

    % time to add to each un-combined tRep to merge it into a combined tRep
    temptRepAdd = zeros(size(Seq.tRep));
    temptRepAdd(tRepsExcitation) = temptRepAdd(tRepsExcitation) - [0,temptRepCombined(1:end-1)] ;
    temptRepAdd = temptRepStart + cumsum(temptRepAdd);

    maptRepToCombined = zeros(size(Seq.tRep));
    maptRepToCombined(tRepsExcitation) = 1;
    maptRepToCombined = cumsum(maptRepToCombined);

    Seq.tRep = temptRepCombined;


    % adjust tReps of k-lines for image reconstruction
    Seq.kLinestReps = maptRepToCombined(Seq.kLinestReps);

    if Seq.CorrectPhase
      % adjust tReps for phase correction
      Seq.AQSlice(2).UsetRep = ...
        reshape(maptRepToCombined(Seq.AQSlice(2).UsetRep), size(Seq.AQSlice(2).UsetRep));
      Seq.Slice(1).UseAtRepetitionTime = ...
        reshape(maptRepToCombined(Seq.Slice(1).UseAtRepetitionTime), ...
                size(Seq.Slice(1).UseAtRepetitionTime));
    end

    Seq.tRepEPISegment = ...
      reshape(maptRepToCombined(Seq.tRepEPISegment), size(Seq.tRepEPISegment));


    % merge pulse program elements into combined tReps
    for t = 1:numel(Grad)
      if ~isempty(Grad(t).Time)
        [Grad(t).Time, I] = ...
          sort(reshape(bsxfun(@plus, Grad(t).Time, temptRepAdd), [], nPulses), 1);
        I = bsxfun(@plus, I, 0:size(I,1):(size(I,1)*(size(I,2)-1)));
        Grad(t).Amp = reshape(Grad(t).Amp, [], nPulses);
        Grad(t).Amp = Grad(t).Amp(I);
        [row, ~] = find(sum(~isnan(Grad(t).Time),2), 1, 'last');
        Grad(t).Time = Grad(t).Time(1:max(row),:);
        Grad(t).Amp = Grad(t).Amp(1:max(row),:);
        % remove (trailing) lines of NaN
        AllIsNaN = all(isnan(Grad(t).Time), 2);
        Grad(t).Time(AllIsNaN,:) = [];
        Grad(t).Amp(AllIsNaN,:) = [];
      end
    end
    for t = 1:numel(TX)
      if ~isempty(TX(t).Start)
        [TX(t).Start, I] = ...
          sort(reshape(bsxfun(@plus, TX(t).Start, temptRepAdd), [], nPulses), 1);
        I = bsxfun(@plus, I, 0:size(I,1):(size(I,1)*(size(I,2)-1)));
        TX(t).Duration = reshape(TX(t).Duration, [], nPulses);
        TX(t).Duration = TX(t).Duration(I);
        TX(t).Amplitude = reshape(TX(t).Amplitude, [], nPulses);
        TX(t).Amplitude = TX(t).Amplitude(I);
        TX(t).Frequency = reshape(TX(t).Frequency, [], nPulses);
        TX(t).Frequency = TX(t).Frequency(I);
        TX(t).Phase = reshape(TX(t).Phase, [], nPulses);
        TX(t).Phase = TX(t).Phase(I);
        % remove (trailing) lines of NaN
        AllIsNaN = all(isnan(TX(t).Start), 2);
        TX(t).Start(AllIsNaN,:) = [];
        TX(t).Duration(AllIsNaN,:) = [];
        TX(t).Amplitude(AllIsNaN,:) = [];
        TX(t).Frequency(AllIsNaN,:) = [];
        TX(t).Phase(AllIsNaN,:) = [];
      end
    end
    for t = 1:numel(AQ)
      if ~isempty(AQ(t).Start)
        [AQ(t).Start, I] = ...
          sort(reshape(bsxfun(@plus, AQ(t).Start, temptRepAdd), [], nPulses), 1);
        I = bsxfun(@plus, I, 0:size(I,1):(size(I,1)*(size(I,2)-1)));
        AQ(t).fSample = reshape(AQ(t).fSample, [], nPulses);
        AQ(t).fSample = AQ(t).fSample(I);
        AQ(t).Frequency = reshape(AQ(t).Frequency, [], nPulses);
        AQ(t).Frequency = AQ(t).Frequency(I);
        AQ(t).Phase = reshape(AQ(t).Phase, [], nPulses);
        AQ(t).Phase = AQ(t).Phase(I);
        AQ(t).nSamples = reshape(AQ(t).nSamples, [], nPulses);
        AQ(t).nSamples = AQ(t).nSamples(I);
        AQ(t).SamplingFactor = reshape(AQ(t).SamplingFactor, [], nPulses);
        AQ(t).SamplingFactor = AQ(t).SamplingFactor(I);
        if ~isemptyfield(AQ(t), 'ResetPhases')
          AQ(t).ResetPhases = reshape(AQ(t).ResetPhases, [], nPulses);
          AQ(t).ResetPhases = any(AQ(t).ResetPhases(I), 1);
        end
        if ~isemptyfield(AQ(t), 'GetData')
          AQ(t).GetData = reshape(AQ(t).GetData, [], nPulses);
          AQ(t).GetData = any(AQ(t).GetData, 1);
        end
        % remove (trailing) lines of NaN
        AllIsNaN = all(isnan(AQ(t).Start), 2);
        AQ(t).Start(AllIsNaN,:) = [];
        AQ(t).fSample(AllIsNaN,:) = [];
        AQ(t).Frequency(AllIsNaN,:) = [];
        AQ(t).Phase(AllIsNaN,:) = [];
        AQ(t).nSamples(AllIsNaN,:) = [];
        AQ(t).SamplingFactor(AllIsNaN,:) = [];
      end
    end
    % "DigitalIO" is an optional field
    if ~isemptyfield(Seq, 'DigitalIO') && ...
        isfield(Seq.DigitalIO(1), 'SetTime') && ...
        isfield(Seq.DigitalIO(1), 'SetValue')
      for t = 1:numel(Seq.DigitalIO)
        if isempty(Seq.DigitalIO(t).SetTime) || isempty(Seq.DigitalIO(t).SetValue)
          continue;
        end
        [Seq.DigitalIO(t).SetTime, I] = ...
          sort(reshape(bsxfun(@plus, Seq.DigitalIO(t).SetTime, temptRepAdd), [], nPulses), 1);
        I = bsxfun(@plus, I, 0:size(I,1):(size(I,1)*(size(I,2)-1)));
        AllIsNaN = all(isnan(Seq.DigitalIO(t).SetTime), 2);
        Seq.DigitalIO(t).SetTime(AllIsNaN,:) = [];
        Seq.DigitalIO(t).SetValue = reshape(Seq.DigitalIO(t).SetValue, [], nPulses);
        Seq.DigitalIO(t).SetValue = Seq.DigitalIO(t).SetValue(I);
        Seq.DigitalIO(t).SetValue(AllIsNaN,:) = [];
        if ~isemptyfield(Seq.DigitalIO(t), 'SetLatenz')
          Seq.DigitalIO(t).SetLatenz = reshape(Seq.DigitalIO(t).SetLatenz, [], nPulses);
          Seq.DigitalIO(t).SetLatenz = Seq.DigitalIO(t).SetLatenz(I);
          Seq.DigitalIO(t).SetLatenz(AllIsNaN,:) = [];
        end
        if ~isemptyfield(Seq.DigitalIO(t), 'Repeat')
          Seq.DigitalIO(t).Repeat = reshape(Seq.DigitalIO(t).Repeat, 1, nPulses);
          Seq.DigitalIO(t).Repeat = all(Seq.DigitalIO(t).Repeat(I), 1);
        end
        % FIXME: Should we bother to remove settings that don't change the IO
        % state within a tRep?
      end
    end
  end


  % map AQ windows in tReps to k-lines of images
  Seq.AQSlice(1).UsetRep = zeros(size(Seq.kLinesImages));
  Seq.AQSlice(1).UseAQWindow = ones(size(Seq.kLinesImages)); % + Seq.CorrectPhase;
  for iImage = 1:Seq.AQSlice(1).nImages
    % identify tReps (and AQ windows) of k-lines of current image
    for iParts = 1:size(Seq.kLinesImages, 3)
      [~, ikLines] = ismember(Seq.kLinesImages(:,iImage,iParts), Seq.kLines);
      Seq.AQSlice(1).UsetRep(ikLines~=0,iImage,iParts) = Seq.kLinestReps(ikLines(ikLines~=0));
      if Seq.SingletRep
        % update UseAQWindow
        Seq.AQSlice(1).UseAQWindow(ikLines~=0,iImage,iParts) = ikLines(ikLines~=0);
      end
    end
  end

  if Seq.SingletRep
    % wrap AQ windows to number of AQ windows

    % all EPI segments have the same length
    Seq.AQSlice(1).UseAQWindow = mod(Seq.AQSlice(1).UseAQWindow, ...
      (Seq.AQSlice(1).oddEvenEchoes+1)*Seq.AQSlice(1).EPIFactor*Seq.AQSlice(1).nImages);
    Seq.AQSlice(1).UseAQWindow(Seq.AQSlice(1).UseAQWindow==0) = ...
      (Seq.AQSlice(1).oddEvenEchoes+1)*Seq.AQSlice(1).EPIFactor*Seq.AQSlice(1).nImages;
    Seq.AQSlice(1).UseAQWindow = Seq.AQSlice(1).UseAQWindow + (Seq.CorrectPhase && ~Seq.CorrectPhaseSeparate);
  end


  if ~any([ Seq.PreProcessSequence, ...
            Seq.StartSequence, ...
            Seq.PollPPGfast, ...
            Seq.GetRawData, ...
            Seq.PostProcessSequence])
    SeqLoop = Seq;
    if isa(HW, 'PD.HWClass')
      SeqLoop.HW = HW.ToStruct();
    else
      SeqLoop.HW = HW;
    end
    SeqLoop.AQ = AQ;
    SeqLoop.TX = TX;
    SeqLoop.Grad = Grad;
    break;
  end

  if Seq.PreProcessSequence
    tStartSequence = Seq.StartSequence;
    Seq.StartSequence = 0;
    Seq.PollPPGfast = 0;
    Seq.GetRawData = 0;
    Seq.PostProcessSequence = 0;

    if tStartSequence
      if Seq.Reinitialize
        Seq.Reinitialize = double((Loop==1)||~Seq.LoopsBreakExactly);
      end
      [~, SeqOut] = set_sequence(HW, Seq, AQ, TX, Grad);  % PreProcessSequence sequence
    else
      [~, SeqLoop] = set_sequence(HW, Seq, AQ, TX, Grad); % PreProcessSequence sequence only
      SeqLoop.TimeToNextSequence = SeqLoop.LoopsBreak + SeqLoop.tOffset(1);
      return;
    end
  end
  if ~isempty(SeqOut.LoopsRepetitionTime) && isempty(SeqOut.LoopsBreak)
    init.Seq.LoopsBreak = SeqOut.LoopsRepetitionTime - SeqOut.SequenceTime;
    SeqOut.LoopsBreak = init.Seq.LoopsBreak;
  end

  if SeqOut.LoopsBreakExactly
    if Loop == 1
      SeqOut.Reinitialize = 1;
      if Seq.Loops > 1
        SeqOut.TimeToNextSequence = SeqOut.LoopsBreak + SeqOut.tOffset(1);
        SeqOut.TimeFromLastSequence = [];
      else
        SeqOut.TimeToNextSequence = [];
        SeqOut.TimeFromLastSequence = [];
      end
    elseif (Loop>=2) && (Loop~=Seq.Loops)
      SeqOut.Reinitialize = 0;
      SeqOut.TimeFromLastSequence = SeqOut.LoopsBreak + SeqOut.tOffset(1);
      SeqOut.TimeToNextSequence = SeqOut.LoopsBreak + SeqOut.tOffset(1);
    elseif Loop == Seq.Loops
      SeqOut.Reinitialize = 0;
      SeqOut.TimeFromLastSequence = SeqOut.LoopsBreak + SeqOut.tOffset(1);
      SeqOut.TimeToNextSequence = [];
    end
  end

  if strcmp(Seq.LoopName{LoopCount}, 'normal')
    if Seq.Loops > 1
      if isempty(SeqOut.LoopsRepetitionTime)
        fprintf('Time to run remaining loops: % 10.1f sec.\n', ...
          SeqOut.SequenceTime*(SeqOut.Loops-Loop+1) + max([Seq.LoopsBreak, 1])*(SeqOut.Loops-Loop+1));
      else
        fprintf('Time to run remaining loops: % 10.1f sec.\n', ...
          SeqOut.LoopsRepetitionTime*(SeqOut.Loops-Loop+1) + SeqOut.SequenceTime*double(SeqOut.Loops==Loop));
      end
    end
  end

  % Run sequence on MMRT and collect data
  SeqOut.PreProcessSequence = 0;
  SeqOut.StartSequence = 1;
  SeqOut.PollPPGfast = 1;
  SeqOut.GetRawData = 1;
  SeqOut.PostProcessSequence = 1;
  [~, SeqOut, data, ~] = set_sequence(HW, SeqOut, AQ, TX, Grad);
  % talker.mySequency.exctractArrayToFile(talker.mySequency.getCommandArray, 'test.txt');

  if ~SeqOut.LoopsBreakExactly && ~isempty(SeqOut.LoopsBreak)
    init.Seq.StartSequenceTime = SeqOut.StartSequenceTime + SeqOut.SequenceTime + SeqOut.LoopsBreak;
  end

  if SeqOut.LoopsBreakExactly
    % Correct counters by timing error counter to get the corrected references
    % for the next loop (or timed measurement).
    init.Seq.EndTimeFPGA = SeqOut.EndTimeFPGA ...
      + SeqOut.TR_Error ./ [SeqOut.HW.MMRT(:).fSystem] * 4;
    init.Seq.LoopCountStart = SeqOut.LoopCountEnd + SeqOut.TR_Error;
  end

  if Seq.CorrectPhase > 0
    % correct frequency offset in encoded acquisition windows
    try
      % get data of (un-encoded) FIDs
      if (SeqOut.LoopPlot || SeqOut.Loops==Loop) && SeqOut.CorrectPlotFrequency
        [dataFID] = get_kSpaceAndImage(data, SeqOut.AQSlice(2));
        [dataFID] = plot_kSpaceAndImage(dataFID, SeqOut.AQSlice(2));
      end
      if Seq.LoopSaveAllData
        for iAQ = 1:numel(data)
          % FIXME: Can an array be transformed to a comma-separated-list without
          % duplicating the memory footprint with a temporary cell?
          data(iAQ).dataWithoutCorrectPhase = data(iAQ);
        end
      end
      [data, SeqOut] = get_CorrectedPhase(data, SeqOut);
    catch ME
      warning('PD:sequence_Flash:ImageReconstructionError', ...
        ['An error occurred during image reconstruction. Trying to continue.\n', ...
        'ATTENTION: The frequency correction (Seq.CorrectPhase) has not been applied yet.\n', ...
        'Error:\n%s'], ...
        getReport(ME));
    end
  end
  % dataCorr = data.data;
  % data.data = data.dataWithoutCorrectPhase;
  % data.data = dataCorr;


  % get data of encoded gradient echoes
  try
    data = get_kSpaceAndImage(data, SeqOut.AQSlice(1));
  catch ME
    warning('PD:sequence_Flash:ImageReconstructionError', ...
      'An error occurred during image reconstruction. Trying to continue.\nError:\n%s', ...
      getReport(ME));
    SeqLoop = SeqOut;
    if isa(HW, 'PD.HWClass')
      SeqLoop.HW = HW.ToStruct();
    else
      SeqLoop.HW = HW;
    end
    SeqLoop.AQ = AQ;
    SeqLoop.TX = TX;
    SeqLoop.Grad = Grad;
    SeqLoop.data = data;
    return;
  end

  if SeqOut.AQSlice(1).nImages > 1
    SeqOut.iLaplace1D.Problem = 'Decay';
  end

  iAQ = 1; % FIXME: Handle data from multiple devices.
  if isZTE
    % FIXME: Do we ever need a ZTE-EPI combination?
    data(iAQ).tImageZ = 0;
  else
    % find k-line without phase encoding (in image indexing)
    iNoPhase = floor(SeqOut.AQSlice(1).nPhase(3)*SeqOut.AQSlice(1).PhaseOS(3)/2)*prod(SeqOut.AQSlice(1).nPhase(1:2).*SeqOut.AQSlice(1).PhaseOS(1:2)) + ...
      floor(SeqOut.AQSlice(1).nPhase(2)*SeqOut.AQSlice(1).PhaseOS(2)/2)*SeqOut.AQSlice(1).nPhase(1)*SeqOut.AQSlice(1).PhaseOS(1) + ...
      floor(SeqOut.AQSlice(1).nPhase(1)*SeqOut.AQSlice(1).PhaseOS(1)/2)+1;

    kLinesNoPhase = Seq.kLinesImages(iNoPhase,:,:);

    % in tRep indexing
    tRepNoPhase = reshape(SeqOut.kLinestReps(kLinesNoPhase), size(kLinesNoPhase));
    if numel(Seq.AQSlice(1).UseAQWindow) == 1
      AQNoPhase = Seq.AQSlice(1).UseAQWindow;
    else
      AQNoPhase = reshape(Seq.AQSlice(1).UseAQWindow(kLinesNoPhase), size(kLinesNoPhase));
    end

    if Seq.SingletRep
      % find corresponding excitation in tRep indexing
      [~, iCol] = find(ismember(Seq.tRepEPISegment(2:end,:), tRepNoPhase));
      tRepExcitation = ...
        reshape(Seq.tRepEPISegment(1,iCol(AQNoPhase-(Seq.CorrectPhase && ~Seq.CorrectPhaseSeparate))), ...
        size(tRepNoPhase));

      tRepCum = cumsum([0, SeqOut.tRep]);
      % echo time for each image
      data(iAQ).tImageZ = ...
        reshape((tRepCum(tRepNoPhase) + SeqOut.Read(1).CenterOfReadout) - (tRepCum(tRepExcitation) + SeqOut.Slice(1).CenterOfPulse), ...
        size(tRepNoPhase, 2), size(tRepNoPhase, 3));

      temptRepAdd = reshape(temptRepAdd, [], nPulses);
      % add time of AQ window in tRep
      data(iAQ).tImageZ = data(iAQ).tImageZ + ...
        reshape(temptRepAdd(sub2ind(size(temptRepAdd), ...
                                    SeqOut.AQSlice(1).UseAQWindow(iNoPhase,:,:,:)+1-(Seq.CorrectPhase && ~Seq.CorrectPhaseSeparate), tRepNoPhase)), ...
        [], size(SeqOut.AQSlice(1).UsetRep, 3));
    else
      % find corresponding excitation in tRep indexing
      [~, iCol] = find(ismember(Seq.tRepEPISegment, tRepNoPhase));
      tRepExcitation = reshape(Seq.tRepEPISegment(1,iCol), size(tRepNoPhase));

      tRepCum = cumsum([0, SeqOut.tRep]);
      % echo time for each image
      data(iAQ).tImageZ = ...
        reshape((tRepCum(tRepNoPhase) + SeqOut.Read(1).CenterOfReadout) - (tRepCum(tRepExcitation) + SeqOut.Slice(1).CenterOfPulse), ...
        size(tRepNoPhase, 2), size(tRepNoPhase, 3));
    end
  end

  % data.Image(max(abs(data(iAQ).Image(:)))/4>abs(data(iAQ).Image(:))) = NaN;
  if Loop == 1
    SeqLoop = SeqOut;
    SeqLoop.data = data(iAQ);
    SeqLoop.data.data = zeros(size(data(iAQ).data));
    SeqLoop.AQSlice(1).raiseFigures = 1;
  end


  switch Seq.LoopName{SeqOut.LoopNameCount}
    case 'CSRLoopPlus'
      % correct slice rephase
      % in first loop just store the data
      SeqLoop.dataLoop(Loop) = deal(data(iAQ));  % FIXME: Why "deal"?

    case 'CSRLoopMinus'
      % correct slice rephase
      SeqLoop.dataLoop(Loop) = deal(data(iAQ));  % FIXME: Why "deal"?
      SeqLoop.AQSlice(1).plotB0ppm = 0;
      SeqLoop.AQSlice(1).plotB0Hz = 0;
      SeqLoop.AQSlice(1).plotB0HzPhase = 0;
      SeqLoop.AQSlice(1).plotB0PpmPhase = 0;
      if Seq.CorrectSliceRephasePlot
        SeqLoop.AQSlice(1).plotImage = 1;
        SeqLoop.AQSlice(1).plotkSpace = 1;
        SeqLoop.AQSlice(1).plotPhase = 1;
      else
        SeqLoop.AQSlice(1).plotImage = 0;
        SeqLoop.AQSlice(1).plotkSpace = 0;
      end
      [SeqLoop.data, SeqLoop.AQSlice(1)] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));

      if Seq.CorrectSliceRephasePlot
        figure(201);
        ax(1) = subplot(2,1,1);
        plot(ax(1), ...
          SeqLoop.data.Ticks(1).Read, abs(SeqLoop.dataLoop(1).Image), ...
          SeqLoop.data.Ticks(1).Read, abs(SeqLoop.dataLoop(2).Image));
        ylim(ax(1), [0, 1].*ylim(ax(1)));
        grid(ax(1), 'on');
        ax(2) = subplot(2,1,2);
        plot(ax(2), ...
          SeqLoop.data.Ticks(1).Read, unwrap3Dmiddle(angle(SeqLoop.dataLoop(1).Image)), ...
          SeqLoop.data.Ticks(1).Read, unwrap3Dmiddle(angle(SeqLoop.dataLoop(2).Image)));
        grid(ax(2), 'on');
        ylim(ax(2), [-2*pi, 2*pi]);
        linkaxes(ax, 'x');
      end

      % set ROI
      if isfield(SeqLoop, 'AQSlice') && isfield(SeqLoop.AQSlice, 'ReadLimMinMax') && ...
          ~isempty(SeqLoop.AQSlice(1).ReadLimMinMax)
        SeqLoop.AQSlice(1).SliceLim(1) = max(SeqLoop.AQSlice(1).SliceLimMinMax(1), ...
          SeqLoop.AQSlice(1).Center2OriginImage(1) - SeqLoop.AQSlice(1).thickness/2);
        SeqLoop.AQSlice(1).SliceLim(2) = min(SeqLoop.AQSlice(1).SliceLimMinMax(2), ...
          SeqLoop.AQSlice(1).Center2OriginImage(1) + SeqLoop.AQSlice(1).thickness/2);

        SeqLoop.AQSlice(1).SliceLim(1) = find(SeqLoop.data.Ticks(1).Read>SeqLoop.AQSlice(1).SliceLim(1), 1, 'first');
        SeqLoop.AQSlice(1).SliceLim(2) = find(SeqLoop.data.Ticks(1).Read<SeqLoop.AQSlice(1).SliceLim(2), 1, 'last');
        SeqLoop.AQSlice(1).SliceLimRoi = zeros(size(SeqLoop.dataLoop(1).Image));
        SeqLoop.AQSlice(1).SliceLimRoi(SeqLoop.AQSlice(1).SliceLim(1):SeqLoop.AQSlice(1).SliceLimI(2)) = 1;
      end
      if isfield(SeqLoop, 'AQSlice') && isemptyfield(SeqLoop.AQSlice, 'SliceLimRoi')
        SeqLoop.AQSlice(1).SliceLimRoi = ...
          abs(SeqLoop.dataLoop(1).Image) > (max(abs(SeqLoop.dataLoop(1).Image))/4) & ...
          abs(SeqLoop.dataLoop(2).Image) > (max(abs(SeqLoop.dataLoop(2).Image))/4);
      end

      SeqLoop.AQSlice(1).SliceLimRoi = conv(double(SeqLoop.AQSlice(1).SliceLimRoi), [1;1;1;1], 'same') == 4;
      SeqLoop.AQSlice(1).SliceLimDiffRoi = SeqLoop.AQSlice(1).SliceLimRoi(1:end-1);

      % calculate slope difference of image phase in read direction
      SeqLoop.data.DiffTicksSlice = diff(SeqLoop.data.Ticks(1).Read);
      SeqLoop.data.MeanDiffTicksSlice = mean(SeqLoop.data.DiffTicksSlice(SeqLoop.AQSlice(1).SliceLimDiffRoi));
      SeqLoop.dataLoop(1).DiffAngleSlice = diff(unwrap(angle(SeqLoop.dataLoop(1).Image)));
      SeqLoop.dataLoop(2).DiffAngleSlice = diff(unwrap(angle(SeqLoop.dataLoop(2).Image)));
      SeqLoop.dataLoop(1).SlopeAngleSlice = ...
        mean(SeqLoop.dataLoop(1).DiffAngleSlice(SeqLoop.AQSlice(1).SliceLimDiffRoi)) ./ SeqLoop.data.MeanDiffTicksSlice;
      SeqLoop.dataLoop(2).SlopeAngleSlice = ...
        mean(SeqLoop.dataLoop(2).DiffAngleSlice(SeqLoop.AQSlice(1).SliceLimDiffRoi)) ./ SeqLoop.data.MeanDiffTicksSlice;
      SeqLoop.data.MeanSlopeAngleSlice = mean([SeqLoop.dataLoop(2).SlopeAngleSlice, SeqLoop.dataLoop(1).SlopeAngleSlice]);
      SeqLoop.data.DiffSlopeAngleSlice = diff([SeqLoop.dataLoop(2).SlopeAngleSlice, SeqLoop.dataLoop(1).SlopeAngleSlice]);
      % calculate offset to slice rephase encoder to compensate slope
      % SeqLoop.data.SliceReadGradTimeIntegralOffset = SeqLoop.data.MeanSlopeAngleSlice / SeqLoop.AQSlice(1).Gamma;
      % SeqLoop.data.SliceReadGradTimeIntegralOffset = 0;
      SeqLoop.data.SliceGradTimeIntegralRephaseOffset = -SeqLoop.data.DiffSlopeAngleSlice / SeqLoop.AQSlice(1).Gamma / 2;
      % SeqLoop.data.SliceGradTimeIntegralRephaseOffset
      if isnan(SeqLoop.data.SliceGradTimeIntegralRephaseOffset)
        warning('isnan(SeqLoop.data.SliceGradTimeIntegralRephaseOffset) set to 0');
        SeqLoop.data.SliceGradTimeIntegralRephaseOffset = 0;
      end

      if Seq.CorrectSliceRephasePlot
        hf = figure(202);
        ax(1) = subplot(2,1,1, 'Parent', hf);
        plot(ax(1), ...
          SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).SliceLimRoi), ...
          abs(SeqLoop.dataLoop(1).Image(SeqLoop.AQSlice(1).SliceLimRoi)), ...
          SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).SliceLimRoi), ...
          abs(SeqLoop.dataLoop(2).Image(SeqLoop.AQSlice(1).SliceLimRoi)));
        ylim(ax(1), [0, 1].*ylim(ax(1)));
        grid(ax(1), 'on');
        ax(2) = subplot(2,1,2, 'Parent', hf);
        plot(ax(2), ...
          SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).SliceLimRoi), ...
          unwrap3Dmiddle(angle(SeqLoop.dataLoop(1).Image(SeqLoop.AQSlice(1).SliceLimRoi))), ...
          SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).SliceLimRoi), ...
          unwrap3Dmiddle(angle(SeqLoop.dataLoop(2).Image(SeqLoop.AQSlice(1).SliceLimRoi))));
        grid(ax(2), 'on');
        linkaxes(ax, 'x');
        title(ax(2), ['SliceGradTimeIntegralRephaseOffset = ', num2str(SeqLoop.data.SliceGradTimeIntegralRephaseOffset), ' T s / m']);
      end

    case 'CRRLoop'
      % correct read rephase

      % don't plot unless requested for this correction
      SeqLoop.AQSlice(1).plotB0ppm = 0;
      SeqLoop.AQSlice(1).plotB0Hz = 0;
      SeqLoop.AQSlice(1).plotB0HzPhase = 0;
      SeqLoop.AQSlice(1).plotB0PpmPhase = 0;
      if Seq.CorrectReadRephasePlot
        SeqLoop.AQSlice(1).plotImage = 1;
        SeqLoop.AQSlice(1).plotkSpace = 1;
        SeqLoop.AQSlice(1).plotPhase = 1;
      else
        SeqLoop.AQSlice(1).plotImage = 0;
        SeqLoop.AQSlice(1).plotkSpace = 0;
      end
      [SeqLoop.data, SeqLoop.AQSlice(1)] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));

      if Seq.CorrectReadRephasePlot
        hf = clf(figure(201));
        ax(1) = subplot(2,1,1, 'Parent', hf);
        plot(ax(1), SeqLoop.data.Ticks(1).Read, squeeze(abs(mean(SeqLoop.data.Image, 6))));
        ylim(ax(1), [0, 1].*ylim(ax(1)));
        grid(ax(1), 'on');
        ax(2) = subplot(2,1,2, 'Parent', hf);
        plot(ax(2), SeqLoop.data.Ticks(1).Read, unwrap3Dmiddle(squeeze(angle(mean(SeqLoop.data.Image, 6)))));
        grid(ax(2), 'on');
        ylim(ax(2), [-2*pi, 2*pi]);
        linkaxes(ax, 'x');
      end

      % set ROI
      tempIm = reshape(SeqLoop.data.Image,SeqLoop.AQSlice(1).nRead,SeqLoop.AQSlice(1).nImages,[]);
      if ~isemptyfield(SeqLoop.AQSlice(1), 'ReadLimMinMax')
        SeqLoop.AQSlice(1).ReadLim(1) = max(SeqLoop.AQSlice(1).ReadLimMinMax(1), ...
          SeqLoop.AQSlice(1).Center2OriginImage(3) - SeqLoop.AQSlice(1).sizeRead/2);
        SeqLoop.AQSlice(1).ReadLim(2) = min(SeqLoop.AQSlice(1).ReadLimMinMax(2), ...
          SeqLoop.AQSlice(1).Center2OriginImage(3) + SeqLoop.AQSlice(1).sizeRead/2);

        SeqLoop.AQSlice(1).ReadLimI(1) = find(SeqLoop.data.Ticks(1).Read>SeqLoop.AQSlice(1).ReadLim(1), 1, 'first');
        SeqLoop.AQSlice(1).ReadLimI(2) = find(SeqLoop.data.Ticks(1).Read<SeqLoop.AQSlice(1).ReadLim(2), 1, 'last');
        SeqLoop.AQSlice(1).ReadLimRoi = zeros(size(tempIm));
        SeqLoop.AQSlice(1).ReadLimRoi(SeqLoop.AQSlice(1).ReadLimI(1):SeqLoop.AQSlice(1).ReadLimI(2),:,:) = 1;
      end
      if isemptyfield(SeqLoop.AQSlice(1), 'ReadLimRoi')
        SeqLoop.AQSlice(1).ReadLimRoi = bsxfun(@gt, abs(tempIm), max(abs(tempIm),[],1)/2);
      end

      SeqLoop.AQSlice(1).ReadLimRoi = convn(double(SeqLoop.AQSlice(1).ReadLimRoi), [1;1;1;1], 'same') == 4;
      SeqLoop.AQSlice(1).ReadLimDiffRoi = SeqLoop.AQSlice(1).ReadLimRoi(1:end-1,:,:);

      % calculate slope of image phase in read direction
      SeqLoop.data.DiffTicksRead = diff(SeqLoop.data.Ticks(1).Read);
      tempIm(~SeqLoop.AQSlice(1).ReadLimRoi) = NaN;
      SeqLoop.data.DiffAngleRead = get_MeanPhaseDiffWeighted(tempIm, 1, 'omitnan');

      for iGradEcho = 1:size(SeqLoop.AQSlice(1).ReadLimDiffRoi, 2)
        SeqLoop.data.MeanDiffTicksRead(iGradEcho) = mean(SeqLoop.data.DiffTicksRead(SeqLoop.AQSlice(1).ReadLimDiffRoi(:,iGradEcho)));
      end
      SeqLoop.data.SlopeAngleRead = bsxfun(@rdivide, SeqLoop.data.DiffAngleRead, SeqLoop.data.MeanDiffTicksRead);
      % calculate offset to read encoder to compensate slope
      SeqLoop.data.ReadGradTimeIntegralOffset = SeqLoop.data.SlopeAngleRead / SeqLoop.HW.GammaDef .* SeqLoop.AQSlice(1).ReadDir(1,:,:);
      SeqLoop.data.ReadGradTimeIntegralOffset = permute(SeqLoop.data.ReadGradTimeIntegralOffset, [3,2,1]); % odd-even are in third dimension
      SeqLoop.data.ReadGradTimeDelayOffset = ...
        SeqLoop.data.ReadGradTimeIntegralOffset(:).' ./ max(abs3D(SeqLoop.Read(1).GradAmp));
      if Seq.CorrectReadRephasePlot
        hf = figure(202);
        clf(hf);
        ax(1) = subplot(2,1,1, 'Parent', hf);
        plot(ax(1), SeqLoop.data.Ticks(1).Read, abs(tempIm(:,:)));
        ylabel(ax(1), 'amplitude');
        ylim(ax(1), [0, 1].*ylim(ax(1)));
        grid(ax(1), 'on');
        ax(2) = subplot(2,1,2, 'Parent', hf);
        plot(ax(2), SeqLoop.data.Ticks(1).Read, unwrap3Dmiddle(angle(tempIm(:,:))));
        ylabel(ax(2), 'phase in rad');
        grid(ax(2), 'on');
        linkaxes(ax, 'x');
        title(ax(2), ...
          {['ReadGradTimeIntegralOffset = ', num2str(SeqLoop.data.ReadGradTimeIntegralOffset(:).'), ' T s / m'], ...
           ['ReadGradTimeDelayOffset = ', num2str(SeqLoop.data.ReadGradTimeDelayOffset*1e6), ' ', char(181) 's ']})
      end

      if any(isnan(SeqLoop.data.ReadGradTimeIntegralOffset(:)))
        warning('isnan(SeqLoop.data.ReadGradTimeIntegralOffset). Setting to 0');
        SeqLoop.data.ReadGradTimeIntegralOffset(isnan(SeqLoop.data.ReadGradTimeIntegralOffset)) = 0;
        SeqLoop.data.ReadGradTimeDelayOffset(isnan(SeqLoop.data.ReadGradTimeIntegralOffset)) = 0;
      end

      if any(abs(SeqLoop.data.ReadGradTimeDelayOffset) > Seq.AQSlice(1).DephasePreEddyCurrentTime)
        warning('ReadGradTimeIntegralOffset is very large. Setting to 0');
        SeqLoop.data.ReadGradTimeIntegralOffset(abs(SeqLoop.data.ReadGradTimeDelayOffset) > Seq.AQSlice(1).DephasePreEddyCurrentTime) = 0;
        SeqLoop.data.ReadGradTimeDelayOffset(abs(SeqLoop.data.ReadGradTimeDelayOffset) > Seq.AQSlice(1).DephasePreEddyCurrentTime) = 0;
      end

    case 'B0map_tEcho1'
      % measure (and correct for) B0 deviation map
      % in first loop just store the data
      SeqLoop.dataLoop(Loop) = deal(data(iAQ));  % FIXME: Why "deal"?
      if SeqLoop.CorrectB0Read.Plot
        if SeqLoop.CorrectB0Read.Plot == 1
          SeqLoop.CorrectB0Read.Plot = 200;
        end
        AQSliceB0step = SeqLoop.AQSlice(1);
        AQSliceB0step = SetCheckFigureHandle(AQSliceB0step, 'plotImage', 50, SeqLoop.CorrectB0Read.Plot);
        AQSliceB0step = SetCheckFigureHandle(AQSliceB0step, 'plotImagePhase', 51, SeqLoop.CorrectB0Read.Plot);
        AQSliceB0step = SetCheckFigureHandle(AQSliceB0step, 'plotkSpace', 52, SeqLoop.CorrectB0Read.Plot);
        AQSliceB0step = SetCheckFigureHandle(AQSliceB0step, 'plotkSpacePhase', 53, SeqLoop.CorrectB0Read.Plot);
        AQSliceB0step = SetCheckFigureHandle(AQSliceB0step, 'plotImageOs', 54, SeqLoop.CorrectB0Read.Plot);
        AQSliceB0step = SetCheckFigureHandle(AQSliceB0step, 'plotImageOsPhase', 55, SeqLoop.CorrectB0Read.Plot);
        % FIXME: Should we zero the plotB0* fields?
        plot_kSpaceAndImage(data(iAQ), AQSliceB0step);
      end

    case 'B0map_tEcho2'
      % measure (and correct for) B0 deviation map
      SeqLoop.dataLoop(Loop) = deal(data(iAQ));  % FIXME: Why "deal"?
      if SeqLoop.CorrectB0Read.Plot
        SeqLoop.CorrectB0Read.Plot = SeqLoop.CorrectB0Read.Plot+100;
        AQSliceB0step = SeqLoop.AQSlice(1);
        AQSliceB0step = SetCheckFigureHandle(AQSliceB0step, 'plotImage', 50, SeqLoop.CorrectB0Read.Plot);
        AQSliceB0step = SetCheckFigureHandle(AQSliceB0step, 'plotImagePhase', 51, SeqLoop.CorrectB0Read.Plot);
        AQSliceB0step = SetCheckFigureHandle(AQSliceB0step, 'plotkSpace', 52, SeqLoop.CorrectB0Read.Plot);
        AQSliceB0step = SetCheckFigureHandle(AQSliceB0step, 'plotkSpacePhase', 53, SeqLoop.CorrectB0Read.Plot);
        AQSliceB0step = SetCheckFigureHandle(AQSliceB0step, 'plotImageOs', 54, SeqLoop.CorrectB0Read.Plot);
        AQSliceB0step = SetCheckFigureHandle(AQSliceB0step, 'plotImageOsPhase', 55, SeqLoop.CorrectB0Read.Plot);
        % FIXME: Should we zero the plotB0* fields?
        plot_kSpaceAndImage(data(iAQ), AQSliceB0step);
      end

      dataB0 = get_B0_map_from_image(HW, SeqLoop, SeqLoop.dataLoop(Loop-1), SeqLoop.dataLoop(Loop));

      if SeqOut.CorrectB0Read.WaitForSample
        disp('Waiting for sample exchange...');
        waitfor(msgbox('Please, change the sample and press "OK" to continue.'));
      end

    case 'normal'
      if Seq.LoopSaveAllData
        % save complete data structure for each averaging step
        if SeqOut.CorrectB0Read.Use
          if ~isfield(SeqOut.AQSlice(1), 'CorrectAmplitude')
            SeqOut.AQSlice(1).CorrectAmplitude = [];
          end
          [data(iAQ), SeqOut.AQSlice(1)] = correct_read_B0(data(iAQ), ...
            SeqOut.AQSlice(1), dataB0, max(sum(SeqOut.Read(1).GradAmp.^2,1).^0.5));
        end
        data(iAQ).StartSequenceTime = SeqOut.StartSequenceTime;
        data(iAQ).fCenter = SeqOut.AQSlice(1).AcquisitionFrequency;
        SeqLoop.dataLoop(Loop) = data(iAQ);
      else
        SeqLoop.dataLoop(Loop).StartSequenceTime = SeqOut.StartSequenceTime;
        SeqLoop.dataLoop(Loop).fCenter = SeqOut.HW.fLarmor;
        SeqLoop.dataLoop(Loop).Image = single(data(iAQ).Image);
      end

      if isfield(SeqOut, 'Correct_foffset')
        SeqLoop.dataLoop(Loop).Correct_foffset = SeqOut.Correct_foffset;
      else
        SeqLoop.dataLoop(Loop).Correct_foffset = 0;
      end
      SeqLoop.dataLoop(Loop).Correct_fCenter = SeqOut.HW.fLarmor;
      if (Seq.CorrectPhase>0) && Seq.LoopSaveAllData
        SeqLoop.dataLoop(Loop).dataWithoutCorrectPhase = data(iAQ).dataWithoutCorrectPhase;
      end

      % average over loops
      SeqLoop.data.data = (SeqLoop.data.data*(Loop-1) + data(iAQ).data)/Loop;

      if Seq.LoopPlot || (Seq.Loops==Loop)  % always plot at last averaging step
        if Seq.LoopPlotAverages || (Seq.Loops==Loop)
          % plot the (current) averaged image
          if Seq.Loops > 1
            try
              SeqLoop.data = get_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
            catch ME
              warning('PD:sequence_Flash:ImageReconstructionError', ...
                'An error occurred during image reconstruction. Trying to continue.\nError:\n%s', ...
                getReport(ME));
              continue;
            end
          end
          if SeqOut.CorrectB0Read.Use
            SeqLoop.data = correct_read_B0(SeqLoop.data, ...
              SeqOut.AQSlice(1), dataB0, max(sum(SeqOut.Read(1).GradAmp.^2,1).^0.5));
          end
          SeqLoop.data.RoI = [];  % re-calculate RoI with averaged data
          try
            [SeqLoop.data, SeqLoop.AQSlice(1)] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
          catch ME
            warning('PD:sequence_Flash:ImageReconstructionError', ...
              'An error occurred while trying to display the results. Trying to continue.\nError:\n%s', ...
              getReport(ME));
            continue;
          end
          SeqLoop.AQSlice(1).raiseFigures = Seq.AQSlice(1).raiseFigures;
        elseif Seq.LoopSaveAllData
          % [SeqLoop.data]=get_kSpaceAndImage(SeqLoop.dataLoop(Loop),SeqLoop.AQSlice(1));
          try
            [~, SeqLoop.AQSlice(1)] = plot_kSpaceAndImage(SeqLoop.dataLoop(Loop), SeqLoop.AQSlice(1));
          catch ME
            warning('PD:sequence_Flash:ImageReconstructionError', ...
              'An error occurred while trying to display the results. Trying to continue.\nError:\n%s', ...
              getReport(ME));
            continue;
          end
          SeqLoop.AQSlice(1).raiseFigures = Seq.AQSlice(1).raiseFigures;
        end
      end

    otherwise
      % FIXME: Issue warning?

  end
% disp('end')


% if Seq.LoopSaveAllData; SeqLoop.dataLoop(Loop).data=data.data; end
% SeqLoop.dataLoop(Loop).Image=single(data.Image);
% SeqLoop.dataLoop(Loop).StartSequenceTime=SeqOut.StartSequenceTime;
% if isfield(SeqOut,'Correct_foffset'); SeqLoop.dataLoop(Loop).Correct_foffset=SeqOut.Correct_foffset;else SeqLoop.dataLoop(Loop).Correct_foffset=0;end
% SeqLoop.dataLoop(Loop).Correct_fCenter=SeqOut.HW.fLarmor;
% if and(Seq.CorrectPhase>0,Seq.LoopSaveAllData); SeqLoop.dataLoop(Loop).dataWithoutCorrectPhase=data.dataWithoutCorrectPhase; end;
% SeqLoop.data.data=(SeqLoop.data.data*(Loop-1)+data.data)/Loop;
% if or(Seq.LoopPlot,Seq.Loops==Loop)
%   [SeqLoop.data]=get_kSpaceAndImage(SeqLoop.data,SeqLoop.AQSlice(1));
%   [SeqLoop.data]=plot_kSpaceAndImage(SeqLoop.data,SeqOut.AQSlice(1));
% end

% %%
% if 0;%isfield(data,'ImageCsi')
% %%
% index=[4,6,4];
% figure(2)
% plot(data.ImageCsiFrequency(:,index(1),index(2),index(3))/HW.fLarmor*1e6-1e6,abs(data.ImageCsi(:,index(1),index(2),index(3))))%%
% xlim([-20 20])
% %%
% index=[5,4,4];
% figure(3)
% plot(data.ImageCsiFrequencyZero(:,index(1),index(2),index(3))/HW.fLarmor*1e6-1e6,abs(data.ImageCsiRawZero(:,index(1),index(2),index(3))))
% xlim([-20 20])
% %%
% end

end

if SeqOut.CorrectB0Read.Get || SeqOut.CorrectB0Read.Use
  SeqLoop.dataB0 = dataB0;
end

end


function AQSlice = SetCheckFigureHandle(AQSlice, field, default, increment)
if ~isemptyfield(AQSlice, field)
  if AQSlice.(field) == 1
    AQSlice.(field) = default;
  end
  if AQSlice.(field) || ishghandle(AQSlice.plotImage, 'figure')
    AQSlice.(field) = double(AQSlice.(field))+increment;
  end
end

end
