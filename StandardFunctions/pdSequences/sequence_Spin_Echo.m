function [SeqLoop, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave)
%% Highly customizable (Turbo) Spin Echo sequence for imaging in 1, 2, or 3 dimensions
%
%     [SeqLoop, mySave] = sequence_Spin_Echo(HW, Seq, AQ, TX, Grad, mySave)
%
% This function executes a (Turbo) Spin Echo imaging experiment and reconstructs
% the corresponding images. Additionally to a highly configurable imaging block,
% the pulse program can be modified or extended (e.g. with encoding blocks)
% using an optional "prepare function" (see Seq.Function_Prepare_Measurement)
% which is executed just before the actual measurement is started.
%
% Optionally, correction measurements can be done immediately before the actual
% image acquisition. Also the image acquisition can be repeated to increase SNR.
% These measurements run in loops.
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
%               Echo time in seconds. (Default: 5e-3)
%
%       T1Estimated
%               Estimated T1 time of the sample in seconds. It is used to set
%               sensible recovery times (breaks) between excitations or echo
%               trains (e.g. RARE). If the particular parameters (see
%               Seq.RepetitionTime, Seq.AQSlice(1).TurboBreak, Seq.LoopsBreak,
%               ...) are also set, they take precedence.
%
%       T2Estimated
%               Estimated T2 time of the sample in seconds. It is used for
%               sensible partioning of echo trains into echo blocks (e.g. RARE
%               or Turbo Spin Echo). If the particular parameters (see
%               Seq.AQSlice(1).TurboFactor, Seq.AQSlice(1).TurboBlocks, ...) are
%               also set, they take precedence.
%
%       RepetitionTime
%               Repetition time (i.e. the time between excitation pulses) in
%               seconds.
%
%       Loops
%               Scalar number of repeated measurements (default: 1). This can be
%               used to increase SNR.
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
%       SteadyState_PreShots90
%               Scalar number of excitations (including refocusing pulses or
%               echo trains) before starting to acquire data.
%
%       SteadyState_PostShots90
%               Scalar number of excitations (including refocusing pulses or
%               echo trains) after data has been acquired.
%
%       SteadyState_PreShots180
%               Scalar number of refocusing (inversion) pulses before starting
%               to acquire data.
%
%       SteadyState_PostShots180
%               Scalar number of refocusing (inversion) pulses after data has
%               been acquired.
%
%       SteadyState_PreShots180SpoilerGrad
%               Boolean value whether to use spoiler gradients around inversion
%               pulses that aren't followed by acquisition windows (see
%               Seq.SteadyState_PreShots180).
%               (Default: 1)
%
%       SteadyState_PreShots180SpoilerGradFirst
%               Boolean value whether to use spoiler gradients around the entire
%               block of inversion pulses that aren't followed by acquisition
%               windows (instead of around each of these inversion pulses). It
%               might be desired to switched this on if
%               Seq.SteadyState_PreShots180SpoilerGrad is set to 0.
%               (Default: Seq.SteadyState_PreShots180SpoilerGrad)
%
%       SingletRep
%               Boolean value to indicate whether the pulse program for each
%               echo train should be realized as single tReps. (default: 0)
%
%       plotSeq
%               control plotting the pulse program. See function "plotSeq" for
%               further details and additional settings.
%
%       plotSeqTR
%               The same as before, but wrap the pulse program at the end of
%               each tRep. See also function "plotSeqTR".
%
%       plotSeqAQ
%               The same as before, but the function tries to wrap the pulse
%               program at sensible points and to display just the crucial part.
%               See also function "plotSeq".
%
%       LoopPlot
%               Plot the image for each loop step (default: 1).
%
%       LoopPlotAverages
%               Calculate and plot the average image at each step (default: 1).
%               Seq.LoopPlot must be true for this to take any effect.
%
%       LoopPlotAll
%               Plot the image aquired in each loop (default: 0). Seq.LoopPlot
%               must be true and Seq.LoopPlotAverages must be false for this to
%               take any effect.
%
%       LoopPlotLastAverage
%               Calculate and plot the average image at the last loop (default:
%               1). If true, this overrides Seq.LoopPlotAll at the last loop
%               step.
%
%       LoopSaveAllData
%               Save the complete data structure (see function "get_data") at
%               each loop step iStep in SeqLoop.dataLoop(iStep). (default:
%               ~Seq.LoopPlotAverages )
%
%       LoopSaveAllSeq
%               Save the actually used Seq structure for each loop step iLoop in
%               SeqLoop.SeqLoop(iLoop). (default: 0).
%
%       CorrectAQWindowPhase
%               Boolean value to indicate whether the phase of the acquisition
%               windows should be corrected (default: 1). This is helpful if the
%               mixer frequency of the acquisition windows doesn't match the
%               Larmor frequency. See also function "get_data".
%
%       CorrectPhase
%               Boolean value. If true, a short acquisition window is placed
%               immediately after the excitation pulse. This signal is used to
%               correct the phase of the acquisition window at the encoded
%               spin echo to compensate for the drift of the Larmor frequency.
%               (default: 0)
%
%       CorrectPhaseDuration
%               Duration of the acquisition window used for Seq.CorrectPhase in
%               seconds. (default: min(0.5e-3, Seq.tEcho/4))
%
%       CorrectPhaseAQtOffset
%               Minimum offset from excitation pulse center to start of the
%               acquisition window for Seq.CorrectPhase in seconds. If
%               necessary, this offset is extented to acommodate for the
%               excitation pulse (+ dead time) or the slice rephase
%               gradient (+ eddy current time). (Default: 0)
%
%       CorrectPhaseNRead
%               Number of samples taken at the acquisition window used for
%               Seq.CorrectPhase. (default: 16)
%
%       CorrectPhase_Set_fLarmorOfNextLoop
%               Boolean value. If true, the Larmor frequency that is used in the
%               next imaging loop is adjusted with the frequency offset that was
%               determined from the last frequency tracking window
%               (Seq.CorrectPhase) of the current imaging loop. (Default: 0)
%
%       CorrectPhaseFrequencyGamma
%               Gyromagnetic ratio in rad/s/T used to calculate the frequency
%               that is used for the correction of a potential frequency drift.
%               See Seq.CorrectPhase. If the set gyromagnetic ratio differs from
%               HW.GammaDef, a dual-frequency rf pulse is used for the
%               excitation pulse. (Default: HW.GammaDef)
%
%       CorrectPhaseFlipAngle
%               Flip angle in degrees of the rf pulse that is used for frequency
%               tracking. This only applies if the frequency used for frequency
%               tracking differs from the frequency used for image acquisition
%               (i.e., if Seq.CorrectPhaseFrequencyGamma does not equal
%               HW.GammaDef). (Default: 90)
%
%       CorrectSliceRephase
%               In an additional pair of loops prior to the actual measurement,
%               use a read out gradient parallel to the slice gradient to detect
%               a delay of the slice gradient pulses (e.g. due to eddy currrents
%               in the pole shoes). See also GradTimeIntegralRephaseOffset in
%               function "get_SliceParameter" (default: 0).
%
%       CorrectSliceRephasePlot
%               Show a plot with the results of the slice rephase correction
%               loops (default: 0).
%
%       CorrectPhaseRephase
%               In an additional pair of loops prior to the actual measurement,
%               use a read out gradient parallel to the phase gradient to detect
%               a delay of the phase gradient pulses (e.g. due to eddy currrents
%               in the pole shoes). See also GradTimeIntegralRephaseOffset in
%               function "get_PhaseParameter" (default: 0).
%
%       CorrectPhaseRephasePlot
%               Show a plot with the results of the phase rephase correction
%               loops (default: 0).
%
%       CorrectReadRephase
%               In an additional loop prior to the actual measurement, detect a
%               delay of the read out gradient by measuring the slope of the
%               phase of the signal during the acquisition window (e.g. due to
%               eddy currrents in the pole shoes). See also GradTimeDelayOffset
%               in function "get_ReadParameter" (default: 0).
%
%       CorrectReadRephasePlot
%               Show a plot with the results of the read rephase correction loop
%               (default: 0).
%
%       DephaseBefore180
%               Move the dephasing gradient pulse to before the refocussing
%               (inversion) pulse. In "non-Turbo" Spin Echo experiments, this
%               allows for a longer read out (acquisition) window. (default: 0)
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
%       CorrectB0Read
%               structure with settings for B0 measurement and read correction.
%               This structure might contain the following fields (default
%               values apply if the fields are omitted or empty):
%           Use
%                   Boolean value (default: false). If it is true the B0 map is
%                   used to compensate the offset in read direction due to B0
%                   deviation. This can also be used with a B0 map that is
%                   acquired in a separate measurement. In this case, the B0
%                   data has to be supplied in Seq.dataB0 (see description of
%                   output values below).
%           Get
%                   Boolean value to set if the B0 map should be acquired
%                   (default: Seq.CorrectB0Read.Use && isempty(Seq.dataB0)). If
%                   it is true, two additional measurements are performed where
%                   the acquisition window is moved forth and back,
%                   respectively. The phase difference between the resulting
%                   images is used to determine the B0 deviation. The results
%                   are returned in SeqLoop.dataB0.
%           tReadoutDiff
%                   difference between the centers of the read out windows of
%                   the two images in seconds (default: 1e-3). The readout
%                   windows and encoding gradients are moved to the front in the
%                   first measurement and to the back in the second measurement
%                   (i.e. the readout is centered at half of this value before
%                   tEcho in the first image, and at half of this value after
%                   Seq.tEcho in the second image).
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
%           thicknessInversion
%                   Thickness of the slice in meter (default: Inf). During the
%                   refocusing (inversion) pulses a gradient is applied such
%                   that the bandwidth of the pulse corresponds to the set
%                   thickness. For a definition of the bandwidth see the used
%                   pulse shape function (see: Seq.AQSlice(1).inversionPulse).
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
%           inversionPulse
%                   Function handle to the pulse shape function used for the
%                   refocusing (inversion) pulses (default: @Pulse_RaisedCos).
%                   Type "Pulse_" followed by hitting the tabulator key to see a
%                   list of available pulse shape functions (all in folder
%                   "pdFunctions"). If you want to add your own pulse shape
%                   function, see function Pulse_Rect_Composite180. Note that
%                   not all pulse shape functions are sensible inversion pulses.
%
%           sizeRead
%                   scalar with the size of the image in read direction in meter
%                   (default: HW.Grad.ImageVol(6)-HW.Grad.ImageVol(5)).
%
%           sizePhase
%                   1x3 vector with the size of the image in phase1, phase2, and
%                   phase3 directions in meter (default according to
%                   HW.Grad.ImageVol). By default, phase3 is parallel to the
%                   read direction and the three phase directions form a right
%                   handed coordinate system.
%
%           nRead
%                   scalar with the integer number of pixels in read direction.
%
%           nPhase
%                   1x3 vector with integer numbers of pixels in phase1, phase2,
%                   and phase3 direction in meter. For a 2d image, set nRead and
%                   nPhase(2) to values > 1. For a 3d image, set nRead,
%                   nPhase(1), and nPhase(2) to values > 1. For a 3d image
%                   (CSI), set nPhase(1), nPhase(2), and nPhase(3) to values >
%                   1.
%
%           Resolution
%                   1x3 vector with the resolution of a 3d image in voxels per
%                   meter (phase1 x phase2 < phase3/read).
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
%                   scalar integer with number of images within (each) Echo
%                   train (Multi-Spin-Echo).
%
%           TurboBlocks
%                   The total number of echoes is evenly distributed over this
%                   scalar integer number of echo trains.
%
%           TurboFactor
%                   Number of image k-lines per excitation.
%
%           TurboBreak
%                   scalar with the duration of the break between turbo blocks
%                   in seconds.
%
%           oddEvenEchoes
%                   Boolean value to acquire separate images consisting of only
%                   odd or even Echoes, respectively. This is helpful to reduce
%                   ghost images when acquiring Multi-Spin-Echo images. (Correct
%                   pulse length errors?). This effectively doubles the number
%                   of echoes per acquired k-line.
%
%           phaseCycling
%                   Boolean value to acquire separate images with the inversion
%                   pulses rotated. This is helpful to reduce the signal of
%                   spins excited by the refocusing (inversion) pulse. This
%                   effectively multiplies the number of echo trains (i.e. also
%                   the number of echoes acquired per k-line).
%
%           phaseCycleSteps
%                   scalar integer with the number of steps in the phase cycle
%                   (default: 2). The steps are evenly distributed for a
%                   complete cycle of 360 degrees between excitation and
%                   refocusing (inversion). If this value is set, it takes
%                   precedence over Seq.AQSlice(1).phaseCycling.
%
%           phaseCycleAngles
%                   list of phase cycle angles in degrees. If this vector is
%                   set, it takes precedence over
%                   Seq.AQSlice(1).phaseCycleSteps.
%                   (Default: (1:Seq.AQSlice(1).phaseCycleSteps).' *
%                             360/Seq.AQSlice(1).phaseCycleSteps )
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
%           excitationFlipAngle
%                   Flip angle of the excitation pulse in degrees (default: 90).
%
%           inversionFlipAngle
%                   Flip angle of the refocusing (inversion) pulse in degrees
%                   (default: 180).
%
%           readoutPhaseInversion
%                   Boolean value to indicate that the phase of the read out
%                   window matches the phase of the preceding inversion pulse.
%                   Otherwise, the read-out phase corresponds to the phase of
%                   the excitation pulse flipped by the inversion pulse.
%                   (Default: false)
%
%           DephaseLengthFactor
%                   scalar value with a factor for the duration of the dephase
%                   pulses (default: Seq.AQSlice.SpoilLengthFactor). A factor of
%                   1 means that the dephase pulse amplitude is equal to the
%                   read-out gradient amplitude (i.e. the duration is
%                   approximately half the acquisition time).
%
%           DephaseTimeOffset
%                   scalar value for the time offset of the dephase pulses in
%                   seconds (default: Seq.AQSlice.SpoilTimeOffset). A value of 0
%                   means that the dephase pulses are immediately before the
%                   read out gradient pulse (both ramps overlap). Positive
%                   values move the dephase pulse away from the read out pulse.
%
%           RephaseLengthFactor
%           RephaseTimeOffset
%                   Same as before but for the rephase gradient pulses. If
%                   RephaseLengthFactor is set to 0, the corresponding gradient
%                   pulses are omitted.
%
%           SpoilFactor
%                   1x3 vector with spoil factors in slice/phase1, phase2, and
%                   read/phase3 direction. A factor of 1 means that the spoilers
%                   dephase the spins over the thickness/pixels/voxels by 360
%                   degrees. By default, a spoiler in phase1/slice direction
%                   is used.
%                   If the spoil factor in read direction is Inf, a spoil factor
%                   is used that leads to the same amplitude for the read
%                   dephase/spoiler pulse and the read out gradient. The total
%                   spoiler size depends on the respective length factors in
%                   this case.
%
%           sizePhaseSpoil
%                   1x3 vector with spoil sizes in slice/phase1, phase2, and
%                   read/phase3 direction in m. Spins dephase by 180 degrees
%                   over this distance. If Seq.AQSlice.SpoilFactor is set at the
%                   same time, this value takes precedence.
%                   Setting sizePhaseSpoil in read direction to 0 has the same
%                   effect like setting SpoilFactor in read direction to Inf.
%
%           SpoilLengthFactor
%                   scalar value with a factor for the duration of the spoiler
%                   pulses (default: 1). A factor of 1 means that the spoiler
%                   has the same duration as the phase gradient pulses (i.e.
%                   approximately half the duration of the readout gradients).
%                   This only applies if Seq.AQSlice.sizePhaseSpoil(3) is 0.
%
%           SpoilFactorEnd
%                   1x3 vector with spoil factors in slice/phase1, phase2, and
%                   read/phase3 direction after the last read-out window in an
%                   echo train. The definition is the same as for
%                   Seq.AQSlice.SpoilFactor.
%                   (Default: Seq.AQSlice.SpoilFactor)
%
%           sizePhaseSpoilEnd
%                   1x3 vector with spoil sizes in slice/phase1, phase2, and
%                   read/phase3 direction in m after the last read-out window in
%                   an echo train. The definition is the same as for
%                   Seq.AQSlice.sizePhaseSpoil.
%                   (Default: Seq.AQSlice.sizePhaseSpoil)
%
%           SpoilTimeOffset
%                   scalar value for the time offset of the spoiler pulses in
%                   seconds (default: 0). A value of 0 means that the spoiler
%                   pulses start immediately after the read out pulse. Positive
%                   values move the spoilers away from the read out pulse.
%
%           SpoilDephaseLengthFactor
%           SpoilDephaseTimeOffset
%                   The same as immediately before but only applies to the
%                   dephase spoilers (i.e. the ones before the refocusing
%                   pulse). By default, the dephase spoilers move with the read
%                   rephase gradients.
%
%           SpoilRephaseLengthFactor
%           SpoilRephaseTimeOffset
%                   The same as immediately before but only applies to the
%                   rephase spoilers (i.e. the ones after the refocusing pulse).
%                   By default, the rephase spoilers move with the read dephase
%                   gradients.
%
%           SliceRephaseLengthFactor
%                   scalar value with a factor to adjust the duration of the
%                   excitation slice rephase pulse. By default, this pulse has
%                   the same duration as the read dephase pulse.
%                   (Default: Seq.AQSlice.DephaseLengthFactor)
%
%           SliceRephaseTimeOffset
%                   time offset of the excitation slice rephase gradient pulse
%                   in seconds (default: Seq.AQSlice.SpoilTimeOffset).
%
%           dualNuclearImage
%                   Boolean value to indicate whether the image should be
%                   acquired simultaneously at two different Gamma values. All
%                   size and resolution settings apply to the (primary)
%                   Seq.AQSlice(1).Gamma. The size and resolution of the image
%                   at the secondary frequency depend on the the ratio between
%                   Seq.AQSlice(1).Gamma and Seq.AQSlice(1).GammaX.
%                   The data for the image at the primary frequency is returned
%                   in SeqLoop.data. The data for the image at the secondary
%                   frequency is returned in SeqLoop.dataX.
%                   (Default: false)
%
%           GammaX
%                   Gyromagnetic ratio of the secondary nucleus in rad/T/s. Only
%                   used if Seq.AQSlice(1).dualNuclearImage is set to true.
%                   (Default: HW.GammaDef)
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
%               Time (in s) measured between excitation pulse and k-line without
%               phase encoding for each image in the echo train.
%
%       StartSequenceTime
%               Approximate time (in Matlab units, see "now") when the
%               experiment started.
%
%     dataLoop
%             Vector of structures containing the data of each (averaging loop)
%             step. Dependend on the settings, this structure contains the
%             following fields:
%             - If Seq.LoopSaveAllData is true:
%                 All fields also contained in SeqLoop.data for each loop step
%             - If Seq.LoopSaveAllData is false:
%       Image, StartSequenceTime, fCenter:
%                 Same as the corresponding fields in SeqLoop.data for each loop
%                 step.
%
%   mySave
%           mySave structure containing (most notably) the current Larmor
%           frequency and timing information for Larmor frequency sweeps.
%
%
% PROGRAMMING NOTES:
%
%   The pulse program is constructed using logical units that correspond to
%   slice, read, and phase encoders.
%   The following functions translate the logical units to the actual pulse
%   program:
%     For slice units:      get_SliceParameter
%     For read units:       get_ReadParameter
%     For phase units:      get_PhaseParameter
%   Please, see the documentation of these functions for further details.
%
%   When manipulating the pulse program with Seq.Function_Prepare_Measurement,
%   the settings for these units can be changed. Currently, the following units
%   are used:
%     Seq.Read(1)
%       Read out at tEcho
%     Seq.Read(2)
%       Dummy read outs at tEcho of 180 degrees pre and post shots
%     Seq.Read(3)
%       Dummy read outs at tEcho of 90 degrees pre and post shots
%     Seq.Read(4)
%       dummy read to move read dephase after excitation pulse to slice rephase
%       for DephaseBefore180
%     Seq.Read(5)
%       not-encoded read out for frequency tracking after excitation pulses for
%       phase correction of encoded readout windows (Seq.CorrectPhase > 0)
%     Seq.Slice(1)
%       90 degrees excitation (with optional slice selection)
%     Seq.Slice(2)
%       180 degrees refocusing (with optional slice selection)
%     Seq.Slice(3)
%       dummy slice used to move the slice dephase pulse of the refocussing
%       pulse to the slice rephase pulse of the excitation pulse for the first
%       inverion.
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
%       Spoiler/crusher at refocusing pulses in phase(1)/slice direction aligned
%       with read dephase and rephase.
%     Seq.Phase(5)
%       Spoiler/crusher at refocusing pulses in phase(2) direction aligned with
%       read dephase and rephase.
%     Seq.Phase(6)
%       Spoiler/crusher at refocusing pulses in phase(3)/read direction aligned
%       with read dephase and rephase.
%     Seq.Phase(7)
%       Spoiler/crusher at end of echo train in phase(1)/slice direction aligned
%       with read rephase (but might be of differing length).
%     Seq.Phase(8)
%       Spoiler/crusher at end of echo train in phase(2) direction aligned with
%       read rephase (but might be of differing length).
%     Seq.Phase(9)
%       Spoiler/crusher at end of echo train in phase(3)/read direction aligned
%       with read rephase (but might be of differing length).
%     Seq.Phase(10)
%       Dummy spoiler that moves the spoiler in phase(1)/slice direction before
%       the first inversion pulse with spoilers to before the first inversion
%       pulse (aligned with excitation slice rephase pulse).
%     Seq.Phase(11)
%       Dummy spoiler that moves the spoiler in phase(2) direction before the
%       first inversion pulse with spoilers to before the first inversion pulse
%       (aligned with excitation slice rephase pulse).
%     Seq.Phase(12)
%       Dummy spoiler that moves the spoiler in phase(3)/read direction before
%       the first inversion pulse with spoilers to before the first inversion
%       pulse (aligned with excitation slice rephase pulse).
%
%
%   See the documentation of "set_sequence" for information about exact timing
%   between subsequent measurements.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2013-2024 Toni Driessle, Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------
%
% See also: get_ReadParameter, get_PhaseParameter, get_SliceParameter,
%           set_sequence, get_data, get_kSpaceAndImage, get_kSpaceAndImageTicks

% FIXME: The documentation of the input still misses some default values.

% LoadSystem                                              % Load system parameters (Reset to default: HW Seq AQ TX)
% % Seq.Loops=1;                                            % Number of loop averages    1...
% % Seq.LoopsBreak=[];                                      % Pause between two loop averages in seconds ([]= fast as possible)
% % Seq.LoopSaveAllData=0;                                  % Save all data of the Loops
% % Seq.LoopPlot = 1;                                       % Plot image after each loop
% % Seq.Find_Frequency_interval = 100;                      % Time between two frequency search
% %
% % Seq.average=1;                                          % Number of averages 1...
% % Seq.averageBreak=1;                                     % Pause between two averages in seconds
% %
% % Seq.plotSeqTR=1:3;                                      % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
% Seq.plotSeq=1:3;                                        % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
% % Seq.plotSeqStart=8;                                     % Plot sequence start with tRep x
% % Seq.plotSeqEnd=12;                                      % Plot sequence stop with tRep x
% %
% Seq.tEcho=5e-3;                                         % Echo time in seconds eg. 5e-3
% Seq.tRep=100e-3;                                        % Repetition time in seconds only used if Seq.TurboBreak=[], auto Seq.tRep=[] and Seq.AQSlice(1).TurboBreak= x sec;
% %
% %
% % % Pixels and size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.AQSlice(1).nRead=32;                                              % Number of Pixels in read, if nRead>1 nPhase(1)=1
% Seq.AQSlice(1).nPhase(1)=1;                                           % Number of Pixels in phase(1), use for 3D, Slice direction
% Seq.AQSlice(1).nPhase(2)=32;                                          % Number of Pixels in phase(2), use for 2D
% Seq.AQSlice(1).nPhase(3)=1;                                           % Number of Pixels in CSI, if nPhase(3)>1 nRead=1 sizeRead=1e12
% Seq.AQSlice(1).HzPerPixMin=1000;                                      % Bandwith per pixel in Hz (1/HzPerPixMin= duration of AQ)
% Seq.AQSlice(1).sizeRead=0.01;                                         % Image size in read in meter (for CSI set to 1e12)
% Seq.AQSlice(1).sizePhase(1)=0.01;                                     % Image size in phase in meter
% Seq.AQSlice(1).sizePhase(2)=0.01;                                     % Image size in phase(2) in meter
% Seq.AQSlice(1).sizePhase(3)=0.01;                                     % Image size in phase(3) in meter
% Seq.AQSlice(1).thickness=0.005;                                       % Image thickness in slice direction  used for 2D and 3D! ([] for no Slice) in meter
% % Seq.AQSlice(1).thicknessInversion=Seq.AQSlice(1).thickness+1e12;    % Invert slice thickness in slice direction  used for 2D and 3D! in meter ([] for no Slice)
% Seq.AQSlice(1).excitationPulse=@Pulse_Sinc_2;                         % excitation pulse function (type "Pulse_" than press tab for selection of pulses)
% Seq.AQSlice(1).inversionPulse=@Pulse_Rect;                            % inversion pulse function
%
% % % Oversampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Seq.AQSlice(1).ReadOS=16;                                           % Oversampling read ([] for automatic, recommended >=2) 1...
% % Seq.AQSlice(1).PhaseOS(1)=1;                                        % Oversampling phase(1)  1... (set to 1 if nPhase(1)=1;
% Seq.AQSlice(1).PhaseOS(2)=2;                                          % Oversampling phase(2)  1... (set to 1 if nPhase(2)=1;
% % Seq.AQSlice(1).PhaseOS(3)=2;                                        % Oversampling phase(3)  1... (set to 1 if nPhase(3)=1;
%
% % % Turbo factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.AQSlice(1).TurboFactor=prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS)/1; % number of echoes per excitation
% Seq.AQSlice(1).TurboBreak=0.15;                                       % break between last echo and next excitation
%
% % % Plot        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Seq.LoopPlot=0;                                            % Plot every Loop
% Seq.AQSlice(1).plotkSpace=0;                               % Plot k-Space
% Seq.AQSlice(1).plotImage=1;                                % Plot Image
% Seq.AQSlice(1).plotPhase=0;                                % Plot phase of k-Space or Image
% Seq.AQSlice(1).ZeroFillWindowSize=1.4;                     % Zero Fill Window Size (k-Space)
% Seq.AQSlice(1).ZeroFillFactor=2;                           % Zero fill resolution factor
% % Seq.AQSlice(1).plotImageHandle=111+Loop;                              % figure handel of the 2D plot
% Seq.AQSlice(1).plotB0HzPhase                                % 3D phase encoded field map
% Seq.AQSlice(1).plotB0PpmPhase                                % 3D phase encoded field map
%
% % % Orientation in Space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Seq.AQSlice(1).alfa=0.0*pi;                     % first rotation around x axis in RAD
% Seq.AQSlice(1).phi=0.0*pi;                      % second rotation around y axis in RAD
% Seq.AQSlice(1).theta=0.0*pi;                    % third rotation around z axis in RAD
% % Seq.AQSlice(1).Center2OriginImage = [0,0,0];  % coordinate of the origin (of the gradient system) in the image coordinate system (slice or phase(1) x phase(2) x read or phase(3)) in m
% % Seq.AQSlice(1).SliceCoordinate=1;             % direction of Slice:   x = 1,  y = 2, z = 3
% % Seq.AQSlice(1).PhaseCoordinate(1)=1;          % direction of Phase(1):   x = 1,  y = 2, z = 3
% % Seq.AQSlice(1).PhaseCoordinate(2)=2;          % direction of Phase(2):   x = 1,  y = 2, z = 3
% % Seq.AQSlice(1).PhaseCoordinate(3)=3;          % direction of Phase(3):   x = 1,  y = 2, z = 3
% % Seq.AQSlice(1).ReadCoordinate=3;              % direction of Read:  x = 1,  y = 2, z = 3
%
% % % Some corrections    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Seq.Function_Prepare_Measurement=[];                                  % Function to change something of the sequence
% % Seq.AQSlice(1).sizePhaseSpoil=0.001;                                  % coding length of spoiler
% % HW.Grad.SliceTimeDelay=[0,0,0]*1e-6;                                  % delay only for the Slice gradient
% % HW.Grad.ReadTimeDelay=[0,0,0]*1e-6;                                   % delay only for the read gradient
% % HW.Grad.PhaseTimeDelay=[0,0,0]*1e-6;                                  % delay only for the phase gradient
% % Seq.LoopsBreakExactly=0;                                              % Lock Time between two loops
% % Seq.AQSlice(1).ReadOSUsedForImage=17;                                 % Number of samples used in CSI
% % Seq.AQSlice(1).TurboFactor=16;                                        % number of echoes per excitation
% % Seq.AQSlice(1).AQWindowDirection=[1,1];                               % [1,-1] to reverse even readouts
% Seq.CorrectAQWindowPhase=1                                              % Correct phase of AQ Window if the AQ.Frequency is not HW.fLarmor
% % Seq.MaxGradAmpSlice=0.1;                                              % limit slice gradient strength
% % Seq.tReadoutOffset=0;                                                 % time shift of readout and gradients (read and phase)
% % Seq.AQSlice(1).excitationFlipAngle=90;                                % excitation Flip Angle in deg
% % Seq.AQSlice(1).inversionFlipAngle=180;                                % inversion Flip Angle in deg
% % Seq.AQSlice(1).SpoilLengthFactor=2;                                   % Spoil Length Factor
% % Seq.CorrectSliceRephase=0;                                            % Correct SliceGradTimeIntegralRephaseOffset
% % Seq.CorrectPhaseRephase=1;                                            % Correct PhaseGradTimeIntegralRephaseOffset
% % Seq.CorrectReadRephase=0;                                             % Correct ReadGradTimeIntegralOffset
% %  Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset=SeqLoop.AQSlice(1).SliceGradTimeIntegralRephaseOffset;
% %  Seq.AQSlice(1).PhaseGradTimeIntegralRephaseOffset=SeqLoop.AQSlice(1).PhaseGradTimeIntegralRephaseOffset;
% %  Seq.AQSlice(1).ReadGradTimeIntegralOffset=SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset;
% Seq.CorrectSliceRephasePlot
% Seq.CorrectPhaseRephasePlot
% Seq.CorrectReadRephasePlot
%
% Spin Echo V.3       29.07.2016

if isemptyfield(Seq, 'PreProcessSequence'), Seq.PreProcessSequence  = 1;  end
if isemptyfield(Seq, 'StartSequence'),      Seq.StartSequence       = 1;  end
if isemptyfield(Seq, 'PollPPGfast'),        Seq.PollPPGfast         = 1;  end
if isemptyfield(Seq, 'GetRawData'),         Seq.GetRawData          = 1;  end
if isemptyfield(Seq, 'PostProcessSequence'),Seq.PostProcessSequence = 1;  end

if Seq.PreProcessSequence
  if ~isfield(Seq, 'AQSlice'),  Seq.AQSlice = struct(); end
  if isemptyfield(Seq, 'Loops'),              Seq.Loops               = 1;  end
  if isemptyfield(Seq, 'SteadyState_PreShotsLoops'), Seq.SteadyState_PreShotsLoops = 0;  end
  if isemptyfield(Seq, 'SteadyState_PreShotsLoopsName'), Seq.SteadyState_PreShotsLoopsName = 'normal_PreShots';  end
  if ~isfield(Seq, 'T1Estimated'),            Seq.T1Estimated         = []; end % estimated T1 in seconds, e.g. oil T1 ~150e-3 sec and water ~2000e-3 sec
  if ~isfield(Seq, 'T2Estimated'),            Seq.T2Estimated         = []; end % estimated T2 in seconds, e.g. oil T2 ~100e-3 sec and water ~2000e-3 sec
  if ~isfield(Seq, 'LoopsBreak'),             Seq.LoopsBreak          = []; end % break between loops in seconds
  if ~isfield(Seq, 'LoopsRepetitionTime'),    Seq.LoopsRepetitionTime = []; end % repetition time of loops in seconds
  if isemptyfield(Seq, 'LoopsBreakExactly'),  Seq.LoopsBreakExactly   = 0;  end
  if ~isfield(Seq, 'StartSequenceTime'),      Seq.StartSequenceTime   = []; end
  if ~isfield(Seq, 'RepetitionTime'),         Seq.RepetitionTime      = []; end % repetition time of excitations in seconds
  if ~isfield(Seq.AQSlice(1), 'TurboBreak'),  Seq.AQSlice(1).TurboBreak = []; end % repetition time in seconds reduced by number of Echoes in train + 1
  if isemptyfield(Seq.AQSlice(1), 'Gamma'),   Seq.AQSlice(1).Gamma = HW.GammaDef; end % Gamma in RAD
  if isemptyfield(Seq, 'CorrectSliceRephasePlot'), Seq.CorrectSliceRephasePlot = 0; end
  if isemptyfield(Seq, 'CorrectPhaseRephasePlot'), Seq.CorrectPhaseRephasePlot = 0; end
  if isemptyfield(Seq, 'CorrectReadRephasePlot'),  Seq.CorrectReadRephasePlot  = 0; end
  if isemptyfield(Seq, 'LoopPlot'),           Seq.LoopPlot            = 1;  end
  if isemptyfield(Seq, 'LoopPlotAverages'),   Seq.LoopPlotAverages    = 1;  end
  if isemptyfield(Seq, 'LoopSeqPlot'),        Seq.LoopSeqPlot         = 0;  end
  if isemptyfield(Seq, 'LoopSaveAllData'),    Seq.LoopSaveAllData     = ~Seq.LoopPlotAverages; end % Save data of each loop
  if isemptyfield(Seq, 'LoopSaveAllSeq'),     Seq.LoopSaveAllSeq      = 0;  end % Save Seq structure used at each loop in field "SeqLoop"
  if isemptyfield(Seq, 'LoopPlotAll'),        Seq.LoopPlotAll         = Seq.LoopSaveAllData; end % plot data of each loop (must set Seq.LoopPlot = 1;)
  if isemptyfield(Seq, 'LoopPlotLastAverage'),Seq.LoopPlotLastAverage = 1;  end % plot average in last loop instead of single measurement
  if isemptyfield(Seq, 'CorrectAQWindowPhase'), Seq.CorrectAQWindowPhase = 1; end  % Correct phase of AQ Window if the AQ.Frequency is not HW.fLarmor
  if isemptyfield(Seq, 'tEcho'),              Seq.tEcho               = 5e-3; end
  Seq.AQSlice(1).tEcho = Seq.tEcho;

  if isemptyfield(Seq.AQSlice(1), 'iDevice'), Seq.AQSlice(1).iDevice  = 1;    end

  if isemptyfield(Seq, 'average'),            Seq.average             = 1;    end
  if isemptyfield(Seq, 'T1'),                 Seq.T1                  = 100e-3; end  % FIXME: Is this used anywhere?

  if isemptyfield(Seq, 'CorrectPhase')
    % use (un-encoded) acquisitions to track magnet frequency
    Seq.CorrectPhase = 0;
  end
  if isemptyfield(Seq, 'CorrectPhase_Set_fLarmorOfNextLoop')
    % use (un-encoded) acquisitions to track magnet frequency
    Seq.CorrectPhase_Set_fLarmorOfNextLoop = 0;
  end
  if isemptyfield(Seq, 'CorrectPhaseDuration')
    Seq.CorrectPhaseDuration = min(0.5e-3, Seq.tEcho/4);
  end
  if isemptyfield(Seq, 'CorrectPhaseAQtOffset'), Seq.CorrectPhaseAQtOffset = 0; end
  if isemptyfield(Seq, 'CorrectPhaseNRead'),  Seq.CorrectPhaseNRead = 16; end
  if isemptyfield(Seq, 'CorrectPhaseFrequencyGamma')
    Seq.CorrectPhaseFrequencyGamma = HW.GammaDef;
  end
  if isemptyfield(Seq, 'CorrectPhaseFlipAngle')
    Seq.CorrectPhaseFlipAngle = 90;
  end
  if isemptyfield(Seq, 'CorrectPlotFrequency'), Seq.CorrectPlotFrequency = 0; end

  if isemptyfield(Seq, 'CorrectRemanence'),   Seq.CorrectRemanence    = 0;    end
  if isemptyfield(Seq, 'CorrectSliceRephase'),Seq.CorrectSliceRephase = 0;    end
  if isemptyfield(Seq, 'CorrectPhaseRephase'),Seq.CorrectPhaseRephase = 0;    end
  if isemptyfield(Seq, 'CorrectReadRephase'), Seq.CorrectReadRephase  = 0;    end
  % if isemptyfield(Seq, 'CorrectPlotFrequency'), Seq.CorrectPlotFrequency = (Seq.CorrectPhase > 0); end

  % structure with data for B0 correction
  if ~isfield(Seq, 'dataB0'), Seq.dataB0 = []; end
  % settings for B0 measurement and correction
  if ~isfield(Seq, 'CorrectB0Read'), Seq.CorrectB0Read = struct(); end
  % correct read offset with B0 map
  if isemptyfield(Seq.CorrectB0Read, 'Use'), Seq.CorrectB0Read.Use = false; end
  % get data for B0 correction (i.e. measure B0 map)
  if isemptyfield(Seq.CorrectB0Read, 'Get'), Seq.CorrectB0Read.Get = Seq.CorrectB0Read.Use && isempty(Seq.dataB0); end
  % readout difference at echo between measurements in s (centered around tEcho)
  if isemptyfield(Seq.CorrectB0Read, 'tReadoutDiff')
    if ~isemptyfield(Seq, 'tReadoutDiff')
      % legacy
      Seq.CorrectB0Read.tReadoutDiff = Seq.tReadoutDiff;
    else
      Seq.CorrectB0Read.tReadoutDiff = 1e-3;
    end
  end
  % maximum frequency offset of B0 map in Hz
  if isemptyfield(Seq.CorrectB0Read, 'MaxFreqOffset'), Seq.CorrectB0Read.MaxFreqOffset = 2000; end
  % minimum amplitude relative to maximum amplitude
  if isemptyfield(Seq.CorrectB0Read, 'MinRelAmp'), Seq.CorrectB0Read.MinRelAmp = 0.2; end
  % maximum relative amplitude deviation from 1.0 between measurements 1 and 2
  if isemptyfield(Seq.CorrectB0Read, 'MaxRelAmpDiff'), Seq.CorrectB0Read.MaxRelAmpDiff = 0.25; end
  % after measuring B0 map, wait for sample exchange
  if isemptyfield(Seq.CorrectB0Read, 'WaitForSample'), Seq.CorrectB0Read.WaitForSample = false; end
  % ZeroFillWindowSize used for B0map
  % FIXME: This value is very low. Increase default value?
  if isemptyfield(Seq.CorrectB0Read, 'ZeroFillWindowSize'), Seq.CorrectB0Read.ZeroFillWindowSize = 0.4; end
  if Seq.CorrectB0Read.Use && ~Seq.CorrectB0Read.Get && isempty(Seq.dataB0)
    warning('PD:sequence_Flash:CorrectB0ReadImpossible', ...
      ['To be able to correct the read error caused by B0 deviations either the ', ...
      'B0 map must be measured (Seq.CorrectB0Read.Get) or a previously measured ', ...
      'B0 map must be provided (Seq.dataB0).\nB0 read correction will be skipped.']);
    Seq.CorrectB0Read.Use = false;
  end

  if isemptyfield(Seq, 'MaxGradAmpSlice'),    Seq.MaxGradAmpSlice     = HW.Grad(Seq.AQSlice(1).iDevice).MaxAmpSlice; end
  if isemptyfield(Seq, 'MaxGradAmpInversion'),    Seq.MaxGradAmpInversion = Seq.MaxGradAmpSlice; end
  if isemptyfield(Seq, 'tReadoutOffset'),     Seq.tReadoutOffset      = 0;    end
  if isemptyfield(Seq, 'SingletRep'),         Seq.SingletRep          = 0;    end
  if ~isfield(Seq,'Function_Prepare_Measurement'), Seq.Function_Prepare_Measurement = []; end
  if ~isfield(Seq.AQSlice(1), 'SliceLimMinMax'),  Seq.AQSlice(1).SliceLimMinMax = [];   end
  if ~isfield(Seq.AQSlice(1), 'PhaseLimMinMax'),  Seq.AQSlice(1).PhaseLimMinMax = [];   end
  if ~isfield(Seq.AQSlice(1), 'ReadLimMinMax'),   Seq.AQSlice(1).ReadLimMinMax  = [];   end
  if isemptyfield(Seq.AQSlice(1), 'AmplitudeUnit'),  Seq.AQSlice(1).AmplitudeUnit  = HW.RX(Seq.AQSlice(1).iDevice).AmplitudeUnit; end
  if isemptyfield(Seq.AQSlice(1), 'AmplitudeUnitScale'), Seq.AQSlice(1).AmplitudeUnitScale = HW.RX(Seq.AQSlice(1).iDevice).AmplitudeUnitScale; end
  if isemptyfield(Seq.AQSlice(1), 'LengthUnit'),  Seq.AQSlice(1).LengthUnit = HW.Grad(Seq.AQSlice(1).iDevice).LengthUnit; end
  if isemptyfield(Seq.AQSlice(1), 'LengthUnitScale'),  Seq.AQSlice(1).LengthUnitScale = HW.Grad(Seq.AQSlice(1).iDevice).LengthUnitScale; end
  if isemptyfield(Seq.AQSlice(1), 'ReadGradSign'),  Seq.AQSlice(1).ReadGradSign = 1; end
  if isemptyfield(Seq.AQSlice(1), 'ReadCoordinate'),  Seq.AQSlice(1).ReadCoordinate = 3; end
  if isemptyfield(Seq.AQSlice(1), 'PhaseCoordinate'),  Seq.AQSlice(1).PhaseCoordinate = [1 2 3]; end
  if isemptyfield(Seq.AQSlice(1), 'PhaseGradSign'),  Seq.AQSlice(1).PhaseGradSign = 1; end
  if isemptyfield(Seq.AQSlice(1), 'SliceGradSign'),  Seq.AQSlice(1).SliceGradSign = 1; end
  if isemptyfield(Seq.AQSlice(1), 'SliceCoordinate'),  Seq.AQSlice(1).SliceCoordinate = 1; end
  if isemptyfield(Seq.AQSlice(1), 'SliceCoordinateInvert'),  Seq.AQSlice(1).SliceCoordinateInvert = Seq.AQSlice(1).SliceCoordinate; end
  if isemptyfield(Seq.AQSlice(1), 'GradSign'),  Seq.AQSlice(1).GradSign = 1;    end

  % fields for plot_kSpaceAndImage
  if isemptyfield(Seq.AQSlice(1), 'plotImageHandle'),  Seq.AQSlice(1).plotImageHandle = 110; end
  if ~isfield(Seq.AQSlice(1), 'iSlice'),             Seq.AQSlice(1).iSlice         = [];   end
  if ~isfield(Seq.AQSlice(1), 'plotImage'),          Seq.AQSlice(1).plotImage      = [];   end
  if ~isfield(Seq.AQSlice(1), 'plotImagePhase'),     Seq.AQSlice(1).plotImagePhase = [];   end
  if ~isfield(Seq.AQSlice(1), 'plotkSpace'),         Seq.AQSlice(1).plotkSpace     = [];   end
  if ~isfield(Seq.AQSlice(1), 'plotkSpacePhase'),    Seq.AQSlice(1).plotkSpacePhase = [];  end
  if ~isfield(Seq.AQSlice(1), 'plotImageOs'),        Seq.AQSlice(1).plotImageOs    = [];   end
  if ~isfield(Seq.AQSlice(1), 'plotImageOsPhase'),   Seq.AQSlice(1).plotImageOsPhase= [];  end
  if isemptyfield(Seq.AQSlice(1), 'plotImagehAxes'), Seq.AQSlice(1).plotImagehAxes = cell(1, 4); end
  if isemptyfield(Seq.AQSlice(1), 'plotImageOshAxes'), Seq.AQSlice(1).plotImageOshAxes = cell(1, 4); end
  if isemptyfield(Seq.AQSlice(1), 'plotDatahAxes'), Seq.AQSlice(1).plotDatahAxes = cell(1, 4); end
  if ~isfield(Seq.AQSlice(1), 'plotB0ppm'),          Seq.AQSlice(1).plotB0ppm      = [];   end
  if ~isfield(Seq.AQSlice(1), 'plotB0PpmPhase'),     Seq.AQSlice(1).plotB0PpmPhase = [];   end
  if ~isfield(Seq.AQSlice(1), 'plotB0ppmGradient'),  Seq.AQSlice(1).plotB0ppmGradient = []; end
  if ~isfield(Seq.AQSlice(1), 'plotB0Gradient'),  Seq.AQSlice(1).plotB0Gradient = []; end
  if isemptyfield(Seq.AQSlice(1), 'plotB0ppmhAxes'), Seq.AQSlice(1).plotB0ppmhAxes = cell(1, 3); end
  if isemptyfield(Seq.AQSlice(1), 'plotB0GradientppmhAxes'),  Seq.AQSlice(1).plotB0GradienthAxes = cell(1, 3); end
  if ~isfield(Seq.AQSlice(1), 'plotB0Hz'),           Seq.AQSlice(1).plotB0Hz       = [];   end
  if ~isfield(Seq.AQSlice(1), 'plotB0HzPhase'),      Seq.AQSlice(1).plotB0HzPhase  = [];   end
  if ~isfield(Seq.AQSlice(1), 'plotB0HzGradient'),   Seq.AQSlice(1).plotB0HzGradient = []; end
  if ~isfield(Seq.AQSlice(1), 'plotB0GradientPhase'),Seq.AQSlice(1).plotB0GradientPhase = []; end
  if isemptyfield(Seq.AQSlice(1), 'plotFft1_data'),  Seq.AQSlice(1).plotFft1_data    = []; end
  if isemptyfield(Seq.AQSlice(1), 'plotData'),       Seq.AQSlice(1).plotData    = []; end
  if isemptyfield(Seq.AQSlice(1), 'plotB0HzhAxes'),  Seq.AQSlice(1).plotB0HzhAxes  = cell(1, 3); end
  if ~isfield(Seq.AQSlice(1), 'plotB1percent'),      Seq.AQSlice(1).plotB1percent  = [];   end
  if isemptyfield(Seq.AQSlice(1), 'plotB1hAxes'),    Seq.AQSlice(1).plotB1hAxes    = cell(1); end
  if ~isfield(Seq.AQSlice(1), 'RoiRelativeValue'),   Seq.AQSlice(1).RoiRelativeValue = []; end
  if ~isfield(Seq.AQSlice(1), 'RoiCutOffPercentile'),Seq.AQSlice(1).RoiCutOffPercentile = []; end
  if ~isfield(Seq.AQSlice(1), 'ZeroFillFactor'),     Seq.AQSlice(1).ZeroFillFactor = [];   end
  if ~isfield(Seq.AQSlice(1), 'ZeroFillWindowSize'), Seq.AQSlice(1).ZeroFillWindowSize = []; end
  if ~isfield(Seq.AQSlice(1), 'sliceomaticProps'),   Seq.AQSlice(1).sliceomaticProps = []; end
  if isemptyfield(Seq.AQSlice(1), 'SliceCartesianAxis'), Seq.AQSlice(1).SliceCartesianAxis = cell(1);    end
  if isemptyfield(Seq.AQSlice(1), 'ReadCartesianAxis'),  Seq.AQSlice(1).ReadCartesianAxis  = cell(1);    end
  if isemptyfield(Seq.AQSlice(1), 'PhaseCartesianAxis'), Seq.AQSlice(1).PhaseCartesianAxis = cell(3,1);  end

  if isemptyfield(Seq.AQSlice(1), 'UseAQWindow'),    Seq.AQSlice(1).UseAQWindow    = 1;    end
  if isemptyfield(Seq.AQSlice(1), 'PhaseOS'),        Seq.AQSlice(1).PhaseOS        = [1 1 1]; end
  if isemptyfield(Seq.AQSlice(1), 'SliceGradTimeIntegralRephaseOffset'), Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset = 0; end
  if isemptyfield(Seq.AQSlice(1), 'PhaseGradTimeIntegralRephaseOffset'), Seq.AQSlice(1).PhaseGradTimeIntegralRephaseOffset = 0; end
  if isemptyfield(Seq.AQSlice(1), 'ReadGradTimeIntegralOffset'), Seq.AQSlice(1).ReadGradTimeIntegralOffset = 0; end
  if isemptyfield(Seq.AQSlice(1), 'ReadTimeDelay'),  Seq.AQSlice(1).ReadTimeDelay  = HW.Grad(Seq.AQSlice(1).iDevice).ReadTimeDelay; end
  if isemptyfield(Seq.AQSlice(1), 'SliceTimeDelay'), Seq.AQSlice(1).SliceTimeDelay = HW.Grad(Seq.AQSlice(1).iDevice).SliceTimeDelay; end
  if isemptyfield(Seq.AQSlice(1), 'PhaseTimeDelay'), Seq.AQSlice(1).PhaseTimeDelay = HW.Grad(Seq.AQSlice(1).iDevice).PhaseTimeDelay; end
  if ~isfield(Seq.AQSlice(1), 'readOutPhase') ...
      || isempty(Seq.AQSlice(1).readOutPhase) ...
      || isnan(Seq.AQSlice(1).readOutPhase)
    Seq.AQSlice(1).readOutPhase = 0;
  end
  if ~isfield(Seq.AQSlice(1), 'readOutPhaseIncrement') ...
      || isempty(Seq.AQSlice(1).readOutPhaseIncrement) ...
      || isnan(Seq.AQSlice(1).readOutPhaseIncrement)
    Seq.AQSlice(1).readOutPhaseIncrement = 0;
  end
  if ~isfield(Seq.AQSlice(1), 'thickness') ...
      || isempty(Seq.AQSlice(1).thickness) ...
      || isnan(Seq.AQSlice(1).thickness)
    Seq.AQSlice(1).thickness = [];
  end
  if ~isfield(Seq.AQSlice(1), 'excitationFlipAngle') ...
      || isempty(Seq.AQSlice(1).excitationFlipAngle) ...
      || any(isnan(Seq.AQSlice(1).excitationFlipAngle))
    Seq.AQSlice(1).excitationFlipAngle = 90;
  end
  if ~isfield(Seq.AQSlice(1), 'excitationFlipAngleIncrement') ...
      || isempty(Seq.AQSlice(1).excitationFlipAngleIncrement) ...
      || isnan(Seq.AQSlice(1).excitationFlipAngleIncrement)
    Seq.AQSlice(1).excitationFlipAngleIncrement = 0;
  end
  if isemptyfield(Seq.AQSlice(1), 'readoutPhaseInversion')
    Seq.AQSlice(1).readoutPhaseInversion = false;
  end
  if ~isfield(Seq.AQSlice(1), 'excitationPhase') ...
      || isempty(Seq.AQSlice(1).excitationPhase) ...
      || isnan(Seq.AQSlice(1).excitationPhase)
    if Seq.AQSlice(1).readoutPhaseInversion
      Seq.AQSlice(1).excitationPhase = 90;
    else
      Seq.AQSlice(1).excitationPhase = 0;
    end
  end
  if ~isfield(Seq.AQSlice(1), 'excitationPhaseIncrement') ...
      || isempty(Seq.AQSlice(1).excitationPhaseIncrement) ...
      || isnan(Seq.AQSlice(1).excitationPhaseIncrement)
    Seq.AQSlice(1).excitationPhaseIncrement = 0;
  end
  if ~isfield(Seq.AQSlice(1), 'inversionFlipAngle') ...
      || isempty(Seq.AQSlice(1).inversionFlipAngle) ...
      || any(isnan(Seq.AQSlice(1).inversionFlipAngle))
    Seq.AQSlice(1).inversionFlipAngle = 180;
  end
  if ~isfield(Seq.AQSlice(1), 'inversionFlipAngleIncrement') ...
      || isempty(Seq.AQSlice(1).inversionFlipAngleIncrement) ...
      || isnan(Seq.AQSlice(1).inversionFlipAngleIncrement)
    Seq.AQSlice(1).inversionFlipAngleIncrement = 0;
  end
  if ~isfield(Seq.AQSlice(1), 'inversionPhase') ...
      || isempty(Seq.AQSlice(1).inversionPhase) ...
      || any(isnan(Seq.AQSlice(1).inversionPhase(:)))
    Seq.AQSlice(1).inversionPhase = 0;
  end
  if ~isfield(Seq.AQSlice(1), 'inversionPhaseIncrement') ...
      || isempty(Seq.AQSlice(1).inversionPhaseIncrement(:)) ...
      || any(isnan(Seq.AQSlice(1).inversionPhaseIncrement(:)))
    Seq.AQSlice(1).inversionPhaseIncrement = 0;
  end
  if ~isfield(Seq.AQSlice(1), 'thicknessInversion') ...
      || isempty(Seq.AQSlice(1).thicknessInversion) ...
      || isnan(Seq.AQSlice(1).thicknessInversion)
    Seq.AQSlice(1).thicknessInversion = Inf;
  end
  if isemptyfield(Seq.AQSlice(1), 'dualNuclearImage')
    Seq.AQSlice(1).dualNuclearImage = false;
  end
  if Seq.AQSlice(1).dualNuclearImage
    if isemptyfield(Seq.AQSlice(1), 'GammaX')
      Seq.AQSlice(1).GammaX = HW.GammaX;  % gyromagnetic ratio of secondary nucleus in rad/T/s
    end
  end

  if isemptyfield(Seq.AQSlice(1), 'plotPhase'),      Seq.AQSlice(1).plotPhase      = 1;    end
  if isemptyfield(Seq.AQSlice(1), 'plotkSpace'),     Seq.AQSlice(1).plotkSpace     = 1;    end
  if isemptyfield(Seq.AQSlice(1), 'plotImage'),      Seq.AQSlice(1).plotImage      = 1;    end
  if isemptyfield(Seq.AQSlice(1), 'raiseFigures'),   Seq.AQSlice(1).raiseFigures   = 0;    end
  if isemptyfield(Seq.AQSlice(1), 'CenterRot'),      Seq.AQSlice(1).CenterRot      = [0,0,0]; end
  if isemptyfield(Seq.AQSlice(1), 'alfa'),           Seq.AQSlice(1).alfa           = 0;    end
  if isemptyfield(Seq.AQSlice(1), 'phi'),            Seq.AQSlice(1).phi            = 0;    end
  if isemptyfield(Seq.AQSlice(1), 'theta'),          Seq.AQSlice(1).theta          = 0;    end
  if isemptyfield(Seq.AQSlice(1), 'angle2Turns'),    Seq.AQSlice(1).angle2Turns = 1/(2*pi); end
  if isemptyfield(Seq.AQSlice(1), 'AcquisitionTime'),Seq.AQSlice(1).AcquisitionTime = Seq.tEcho/4; end
  if isemptyfield(Seq.AQSlice(1), 'HzPerPixMin'),    Seq.AQSlice(1).HzPerPixMin    = 1/Seq.AQSlice(1).AcquisitionTime; end

  % image size
  if isemptyfield(Seq.AQSlice(1), 'sizeRead'),  Seq.AQSlice(1).sizeRead = 0; end
  if isemptyfield(Seq.AQSlice(1), 'sizePhase'), Seq.AQSlice(1).sizePhase = [0,0,0]; end
  if Seq.AQSlice(1).sizePhase(1) == 0
    Seq.AQSlice(1).sizePhase(1) = HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(2) - HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(1);
  end
  if numel(Seq.AQSlice(1).sizePhase) == 1 || Seq.AQSlice(1).sizePhase(2) == 0
    Seq.AQSlice(1).sizePhase(2) = HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(4) - HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(3);
  end
  if numel(Seq.AQSlice(1).sizePhase) == 2 || Seq.AQSlice(1).sizePhase(3) == 0
    Seq.AQSlice(1).sizePhase(3) = HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(6) - HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(5);
  end
  if Seq.AQSlice(1).sizeRead == 0
    Seq.AQSlice(1).sizeRead = HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(6) - HW.Grad(Seq.AQSlice(1).iDevice).ImageVol(5);
  end

  % Resolution
  if isemptyfield(Seq.AQSlice(1), 'Resolution')
    Seq.AQSlice(1).Resolution = [];
  else
    switch numel(Seq.AQSlice(1).Resolution)
      case 1
        if ~isemptyfield(Seq.AQSlice(1), 'thickness') && Seq.AQSlice(1).thickness < 1000
          Seq.AQSlice(1).Resolution=[inf, repmat(Seq.AQSlice(1).Resolution(1), 1, 2)];      % 2D iso
        else
          Seq.AQSlice(1).Resolution=repmat(Seq.AQSlice(1).Resolution(1), 1, 3);             % 1D 2D 3D CSI iso
        end
      case 2
          Seq.AQSlice(1).Resolution=[inf, Seq.AQSlice(1).Resolution(1),Seq.AQSlice(1).Resolution(2)]; % 2D
    end
    if Seq.AQSlice(1).sizePhase(1)<1000;
      Seq.AQSlice(1).nPhase(1)=round(Seq.AQSlice(1).sizePhase(1)/Seq.AQSlice(1).Resolution(1));
      Seq.AQSlice(1).sizePhase(1)=Seq.AQSlice(1).nPhase(1)*Seq.AQSlice(1).Resolution(1);
    else
      if isemptyfield(Seq.AQSlice(1), 'thickness');Seq.AQSlice(1).thickness=Seq.AQSlice(1).Resolution(1); end
      Seq.AQSlice(1).nPhase(1)=1;
    end
    Seq.AQSlice(1).nPhase(2)=round(Seq.AQSlice(1).sizePhase(2)/Seq.AQSlice(1).Resolution(2));
    Seq.AQSlice(1).sizePhase(2)=Seq.AQSlice(1).nPhase(2)*Seq.AQSlice(1).Resolution(2);
    if Seq.AQSlice(1).sizePhase(3)<Seq.AQSlice(1).sizeRead;
      Seq.AQSlice(1).nPhase(3)=round(Seq.AQSlice(1).sizePhase(3)/Seq.AQSlice(1).Resolution(3));
      Seq.AQSlice(1).sizePhase(3)=Seq.AQSlice(1).nPhase(3)*Seq.AQSlice(1).Resolution(3);
      Seq.AQSlice(1).nRead=1;
    else
      Seq.AQSlice(1).nRead=round(Seq.AQSlice(1).sizeRead/Seq.AQSlice(1).Resolution(3));
      Seq.AQSlice(1).sizeRead=Seq.AQSlice(1).nRead*Seq.AQSlice(1).Resolution(3);
      Seq.AQSlice(1).nPhase(3)=1;
    end
    Seq.AQSlice(1).nPhase(isinf(Seq.AQSlice(1).nPhase)|isnan(Seq.AQSlice(1).nPhase)|(Seq.AQSlice(1).nPhase)==0)=1;
    Seq.AQSlice(1).nRead(isinf(Seq.AQSlice(1).nRead)|isnan(Seq.AQSlice(1).nRead)|(Seq.AQSlice(1).nRead)==0)=1;
    Seq.AQSlice(1).sizePhase(isinf(Seq.AQSlice(1).sizePhase)|isnan(Seq.AQSlice(1).sizePhase)|(Seq.AQSlice(1).sizePhase)==0)=inf;
    Seq.AQSlice(1).sizeRead(isinf(Seq.AQSlice(1).sizeRead)|isnan(Seq.AQSlice(1).sizeRead)|(Seq.AQSlice(1).sizeRead)==0)=inf;
  end

  if isemptyfield(Seq.AQSlice(1), 'nRead'), Seq.AQSlice(1).nRead = 32; end
  if isemptyfield(Seq.AQSlice(1), 'nPhase'), Seq.AQSlice(1).nPhase = [1 32 1]; end
  % dimensions of AQSlice.nPhase
  if Seq.AQSlice(1).nPhase(1)==0,                                    Seq.AQSlice(1).nPhase(1) = 1; end
  if numel(Seq.AQSlice(1).nPhase)==1 || Seq.AQSlice(1).nPhase(2)==0, Seq.AQSlice(1).nPhase(2) = 1; end
  if numel(Seq.AQSlice(1).nPhase)==2 || Seq.AQSlice(1).nPhase(3)==0, Seq.AQSlice(1).nPhase(3) = 1; end
  % dimensions of AQSlice.PhaseOS
  if isemptyfield(Seq.AQSlice(1), 'PhaseOS'), Seq.AQSlice(1).PhaseOS = [1 1 1]; end
  if Seq.AQSlice(1).PhaseOS(1)==0, Seq.AQSlice(1).PhaseOS(1) = 1; end
  if numel(Seq.AQSlice(1).PhaseOS)==1 || Seq.AQSlice(1).PhaseOS(2)==0, Seq.AQSlice(1).PhaseOS(2) = 1; end
  if numel(Seq.AQSlice(1).PhaseOS)==2 || Seq.AQSlice(1).PhaseOS(3)==0, Seq.AQSlice(1).PhaseOS(3) = 1; end

  Seq.AQSlice(1).sizePhase(Seq.AQSlice(1).nPhase==1) = Inf;
  Seq.AQSlice(1).sizeRead(Seq.AQSlice(1).nRead==1) = Inf;
  if isemptyfield(Seq.AQSlice(1), 'thickness'), Seq.AQSlice(1).thickness = Inf; end

  % (approximate) extension of the RoI for B0 map in meters
  if isemptyfield(Seq.CorrectB0Read, 'RoIExtension')
    Seq.CorrectB0Read.RoIExtension = 0.15 * [Seq.AQSlice(1).sizeRead, Seq.AQSlice(1).sizePhase(1:2)];
  end
  if isscalar(Seq.CorrectB0Read.RoIExtension)
    Seq.CorrectB0Read.RoIExtension(1:3) = Seq.CorrectB0Read.RoIExtension;
  end

  Seq.AQSlice(1).PhaseOS=round(Seq.AQSlice(1).nPhase.*Seq.AQSlice(1).PhaseOS)./Seq.AQSlice(1).nPhase; % integer number of phase steps
  if ~isfield(Seq, 'plotSeqEnd'), Seq.plotSeqEnd = []; end
  % if isemptyfield(Seq.AQSlice(1), 'UsetRep')
  %   Seq.AQSlice(1).UsetRep = Seq.SteadyState_PreShots+(1:Seq.nEchos.*prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS));
  % end
  autoReadOS = false;
  if isemptyfield(Seq.AQSlice(1), 'ReadOS')
    Seq.AQSlice(1).ReadOS = max(2, ...
      ceil(HW.RX(Seq.AQSlice(1).iDevice).fSample / ...
        min(8000, HW.RX(Seq.AQSlice(1).iDevice).CIC_Decimation_Max) / ...
        (Seq.AQSlice(1).HzPerPixMin*Seq.AQSlice(1).nRead)));
    autoReadOS = true;
  end
  if isemptyfield(Seq.AQSlice(1), 'SamplingFactor')
    % integer sampling factor that is used to (down-)sample the received signal
    % (before image reconstruction)
    minReadOS = ...
      ceil(HW.RX(Seq.AQSlice(1).iDevice).fSample / ...
        min(8000, HW.RX(Seq.AQSlice(1).iDevice).CIC_Decimation_Max) / ...
        (Seq.AQSlice(1).HzPerPixMin*Seq.AQSlice(1).nRead*Seq.AQSlice(1).ReadOS))*Seq.AQSlice(1).ReadOS;
    if minReadOS > Seq.AQSlice(1).ReadOS
      Seq.AQSlice(1).SamplingFactor = ceil(minReadOS/Seq.AQSlice(1).ReadOS);
    else
      Seq.AQSlice(1).SamplingFactor = 1;
    end
  end
  if autoReadOS
    Seq.AQSlice(1).ReadOS = ceil(Seq.AQSlice(1).ReadOS / Seq.AQSlice(1).SamplingFactor);
  end
  if isemptyfield(Seq.AQSlice(1), 'ReadOSUsedForImage'), Seq.AQSlice(1).ReadOSUsedForImage = min(Seq.AQSlice(1).ReadOS, 17); end

  numKLines = prod(Seq.AQSlice(1).nPhase) * prod(Seq.AQSlice(1).PhaseOS);  % k-lines per image
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
        error('PD:sequence_Spin_Echo:incompatiblekLineOrderType', ...
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
      error('PD:sequence_Spin_Echo:incompatiblekLineOrder', ...
        '"AQSlice.kLineOrder" must be a permutation or selection vector with the numbers of kLines to be acquired.');
    end
  end

  if ~isfield(Seq.AQSlice(1), 'excitationPulse'),              Seq.AQSlice(1).excitationPulse = [];   end
  if ~isfield(Seq.AQSlice(1), 'excitationFlipAngleComposite'), Seq.AQSlice(1).excitationFlipAngleComposite = [];   end
  if ~isfield(Seq.AQSlice(1), 'excitationFlipPhaseComposite'), Seq.AQSlice(1).excitationFlipPhaseComposite = [];   end
  if isemptyfield(Seq.AQSlice(1), 'inversionPulse'), Seq.AQSlice(1).inversionPulse = @Pulse_RaisedCos; end
  if ~isfield(Seq.AQSlice(1), 'inversionFlipAngleComposite'),  Seq.AQSlice(1).inversionFlipAngleComposite = [];   end
  if ~isfield(Seq.AQSlice(1), 'inversionFlipPhaseComposite'),  Seq.AQSlice(1).inversionFlipPhaseComposite = [];   end

  % Spoil factor in slice, phase and read direction. A spoil factor of 1 means
  % that the spins at the border of the slice/pixel/voxel de-phase 2*pi with
  % respect to each other.
  if isemptyfield(Seq.AQSlice(1), 'SpoilFactor')
    if Seq.AQSlice(1).nPhase(1) > 1
      Seq.AQSlice(1).SpoilFactor = [1, 0, 0];
    else
      % Spoil only in slice direction
      Seq.AQSlice(1).SpoilFactor = [min(20,min(Seq.AQSlice(1).thickness,Seq.AQSlice(1).thicknessInversion)/0.0002), 0, 0];
    end
  end
  if isemptyfield(Seq.AQSlice(1), 'sizePhaseSpoil')
    % Seq.AQSlice(1).sizePhaseSpoil = 0.5.*min(abs([Seq.AQSlice(1).sizePhase(:)./Seq.AQSlice(1).nPhase(:)./(double(Seq.AQSlice(1).nPhase(:)~=1)), ...
    %   [Seq.AQSlice(1).thickness./double(Seq.AQSlice(1).thickness<1000);Inf;Seq.AQSlice(1).sizeRead./Seq.AQSlice(1).nRead./(double(Seq.AQSlice(1).nRead~=1))]]),[],2);
    Seq.AQSlice(1).sizePhaseSpoil = ...
      1 ./ ...
      ((Seq.AQSlice(1).nPhase>1).*(Seq.AQSlice(1).nPhase)./Seq.AQSlice(1).sizePhase*2 + ...
      (Seq.AQSlice(1).nPhase==1).*[max(2 / min(Seq.AQSlice(1).thickness, Seq.AQSlice(1).thicknessInversion), ...
                                       2 / mean(abs(HW.Grad(Seq.AQSlice(1).iDevice).ImageVol))), ...
                                   0, Seq.AQSlice(1).nRead/Seq.AQSlice(1).sizeRead*2]) ...
      ./ Seq.AQSlice(1).SpoilFactor;
  elseif numel(Seq.AQSlice(1).sizePhaseSpoil) == 1
    Seq.AQSlice(1).sizePhaseSpoil(1:3) = Seq.AQSlice(1).sizePhaseSpoil;
    Seq.AQSlice(1).SpoilFactor = NaN(1, 3);
  else
    Seq.AQSlice(1).SpoilFactor = NaN(1, 3);
  end
  if isemptyfield(Seq.AQSlice(1), 'SpoilLengthFactor'), Seq.AQSlice(1).SpoilLengthFactor = 1; end

  % Optionally differing values for the spoiler after the last read-out window
  % in each echo train.
  SpoilFactorEndDefined = true;
  if isemptyfield(Seq.AQSlice(1), 'SpoilFactorEnd')
    Seq.AQSlice(1).SpoilFactorEnd = Seq.AQSlice(1).SpoilFactor;
    SpoilFactorEndDefined = false;
  end
  if isemptyfield(Seq.AQSlice(1), 'sizePhaseSpoilEnd')
    if ~SpoilFactorEndDefined
      % If user didn't specifically define SpoilFactorEnd, use identical
      % spoilers during the echo train and at the end.
      Seq.AQSlice(1).sizePhaseSpoilEnd = Seq.AQSlice(1).sizePhaseSpoil;
    else
      Seq.AQSlice(1).sizePhaseSpoilEnd = ...
        1 ./ ...
        ((Seq.AQSlice(1).nPhase>1).*(Seq.AQSlice(1).nPhase)./Seq.AQSlice(1).sizePhase*2 + ...
        (Seq.AQSlice(1).nPhase==1).*[max(2 / min(Seq.AQSlice(1).thickness, Seq.AQSlice(1).thicknessInversion), ...
                                         2 / mean(abs(HW.Grad(Seq.AQSlice(1).iDevice).ImageVol))), ...
                                     0, Seq.AQSlice(1).nRead/Seq.AQSlice(1).sizeRead*2]) ...
        ./ Seq.AQSlice(1).SpoilFactorEnd;
    end
  elseif numel(Seq.AQSlice(1).sizePhaseSpoilEnd) == 1
    Seq.AQSlice(1).sizePhaseSpoilEnd(1:3) = Seq.AQSlice(1).sizePhaseSpoilEnd;
    Seq.AQSlice(1).SpoilFactorEnd = NaN(1, 3);
  else
    Seq.AQSlice(1).SpoilFactorEnd = NaN(1, 3);
  end

  if isemptyfield(Seq.AQSlice(1), 'SpoilTimeOffset'), Seq.AQSlice(1).SpoilTimeOffset = 0; end
  % read, phase and slice rephase and dephase pulses adapt to SpoilLengthFactor
  % and SpoilTimeOffset by default
  if isemptyfield(Seq.AQSlice(1), 'DephaseLengthFactor')
    Seq.AQSlice(1).DephaseLengthFactor = Seq.AQSlice(1).SpoilLengthFactor;
  end
  if isemptyfield(Seq.AQSlice(1), 'DephaseTimeOffset')
    Seq.AQSlice(1).DephaseTimeOffset = Seq.AQSlice(1).SpoilTimeOffset;
  end
  if isemptyfield(Seq.AQSlice(1), 'RephaseLengthFactor')
    Seq.AQSlice(1).RephaseLengthFactor = Seq.AQSlice(1).SpoilLengthFactor;
  end
  if isemptyfield(Seq.AQSlice(1), 'RephaseTimeOffset')
    Seq.AQSlice(1).RephaseTimeOffset = Seq.AQSlice(1).SpoilTimeOffset;
  end
  if isemptyfield(Seq.AQSlice(1), 'SliceRephaseTimeOffset')
    Seq.AQSlice(1).SliceRephaseTimeOffset = Seq.AQSlice(1).SpoilTimeOffset;
  end
  if isemptyfield(Seq.AQSlice(1), 'SliceRephaseLengthFactor')
    Seq.AQSlice(1).SliceRephaseLengthFactor = 1;
  end
  % spoiler pulses move with rephase and dephase pulses by default
  if isemptyfield(Seq.AQSlice(1), 'SpoilDephaseLengthFactor')
    Seq.AQSlice(1).SpoilDephaseLengthFactor = Seq.AQSlice(1).SpoilLengthFactor/Seq.AQSlice(1).RephaseLengthFactor;
  end
  if isemptyfield(Seq.AQSlice(1), 'SpoilDephaseTimeOffset')
    Seq.AQSlice(1).SpoilDephaseTimeOffset = Seq.AQSlice(1).RephaseTimeOffset;
  end
  if isemptyfield(Seq.AQSlice(1), 'SpoilRephaseLengthFactor')
    Seq.AQSlice(1).SpoilRephaseLengthFactor = Seq.AQSlice(1).SpoilLengthFactor/Seq.AQSlice(1).DephaseLengthFactor;
  end
  if isemptyfield(Seq.AQSlice(1), 'SpoilRephaseTimeOffset')
    Seq.AQSlice(1).SpoilRephaseTimeOffset = Seq.AQSlice(1).DephaseTimeOffset;
  end
  if isemptyfield(Seq.AQSlice(1), 'SpoilDephaseLengthFactorEnd')
    Seq.AQSlice(1).SpoilDephaseLengthFactorEnd = Seq.AQSlice(1).SpoilLengthFactor/Seq.AQSlice(1).RephaseLengthFactor;
  end


  if isemptyfield(Seq, 'DephaseBefore180'), Seq.DephaseBefore180 = 0; end
  Seq.DephaseBefore180 = (Seq.DephaseBefore180 > 0);  % convert to logical
  if isemptyfield(Seq.AQSlice(1), 'SpoilFactorInversionBlock'), Seq.AQSlice(1).SpoilFactorInversionBlock = 0; end
  if isemptyfield(Seq.AQSlice(1), 'sizeSpoilInversionBlock')
    Seq.AQSlice(1).sizeSpoilInversionBlock = Seq.AQSlice(1).thicknessInversion/Seq.AQSlice(1).SpoilFactorInversionBlock;
  else
    Seq.AQSlice(1).SpoilFactorInversionBlock = Seq.AQSlice(1).thicknessInversion/Seq.AQSlice(1).sizeSpoilInversionBlock;
  end
  % acquire separate images consisting of only odd or even Echoes respectively (correct pulse length errors?)
  if isemptyfield(Seq.AQSlice(1), 'oddEvenEchoes'),  Seq.AQSlice(1).oddEvenEchoes = 0; end
  Seq.AQSlice(1).oddEvenEchoes = double(Seq.AQSlice(1).oddEvenEchoes>0);
  % acquire separate images with the inversion pulses rotated
  if isemptyfield(Seq.AQSlice(1), 'phaseCycling'),   Seq.AQSlice(1).phaseCycling = 0; end
  % number of acquire separate images with the inversion pulses rotated by 360/phaseCycleSteps
  if isemptyfield(Seq.AQSlice(1), 'phaseCycleSteps')
    if ~isemptyfield(Seq.AQSlice(1), 'phaseCycleAngles')
      Seq.AQSlice(1).phaseCycleSteps = numel(Seq.AQSlice(1).phaseCycleAngles);
    else
      if Seq.AQSlice(1).phaseCycling
        Seq.AQSlice(1).phaseCycleSteps = 2;
      else
        Seq.AQSlice(1).phaseCycleSteps = 1;
      end
    end
  end
  % explicit list of phase cycle angles
  if isemptyfield(Seq.AQSlice(1), 'phaseCycleAngles')
    Seq.AQSlice(1).phaseCycleAngles = (1:Seq.AQSlice(1).phaseCycleSteps).' * 360/Seq.AQSlice(1).phaseCycleSteps;
  end
  Seq.AQSlice(1).phaseCycleAngles = Seq.AQSlice(1).phaseCycleAngles(:);
  Seq.AQSlice(1).phaseCycleSteps = numel(Seq.AQSlice(1).phaseCycleAngles);

  if ~isemptyfield(Seq, 'nEchos') % number of images within Echo train (deprecated)
    if ~isemptyfield(Seq.AQSlice(1), 'TurboFactor') && Seq.AQSlice(1).TurboFactor > 1
      warning('PD:sequence_Spin_Echo:nEchos_TurboFactor', ...
        'The definition of the "TurboFactor" in combination with "nEchos" has changed from prior versions.\nConsider using "AQSlice.nImages" instead of "nEchos".');
    end
    if isemptyfield(Seq.AQSlice(1), 'nImages')
      Seq.AQSlice(1).nImages = Seq.nEchos;
    else
      warning('PD:sequence_Spin_Echo:nEchos_TurboFactor', ...
        'The usage of "nEchos" is deprecated. If "AQSlice.nImages" is also set, the value of "nEchos" is ignored.');
    end
  end
  if isemptyfield(Seq.AQSlice(1), 'nImages')  % number of images within Echo train (Multi-Spin-Echo)
    Seq.AQSlice(1).nImages = 1;
  end

  if Seq.AQSlice(1).RephaseLengthFactor == 0 && ...
      (Seq.AQSlice(1).TurboFactor > 1 || Seq.AQSlice(1).oddEvenEchoes || Seq.AQSlice(1).nImages > 1)
    error('PD:sequence_Spin_Echo:rephaseNeeded', ...
      'Rephasing the spin system cannot be skipped when there are multiple echoes in the echo train.');
  end

  if ~isempty(Seq.T2Estimated) && ~isemptyfield(Seq.AQSlice(1), 'TurboBlocks')
    % FIXME: Do we need this warning?
    warning('PD:sequence_Spin_Echo:AmbiguousTurboBlocks', ...
      'Seq.AQSlice(1).TurboBlocks takes precedence over Seq.T2Estimated.');
  end
  nPhase = round(prod(Seq.AQSlice(1).nPhase.*Seq.AQSlice(1).PhaseOS));
  % Turbo factor
  if isemptyfield(Seq.AQSlice(1), 'TurboBlocks')
    if ~isempty(Seq.T2Estimated)
      Seq.AQSlice(1).TurboBlocks = ceil(nPhase*Seq.AQSlice(1).nImages*(Seq.AQSlice(1).oddEvenEchoes+1)*Seq.tEcho/Seq.T2Estimated); % break Echo train before Seq.T2Estimated
    else
      Seq.AQSlice(1).TurboBlocks = nPhase;  % one k-line of each image per Echo train (Spin Echo)
    end
  end
  Seq.AQSlice(1).TurboBlocks = min([nPhase*Seq.AQSlice(1).nImages, Seq.AQSlice(1).TurboBlocks]);  % limit to nAQperTrain
  Seq.AQSlice(1).TurboBlocks = Seq.AQSlice(1).TurboBlocks-1+find(~rem(nPhase./(Seq.AQSlice(1).TurboBlocks:nPhase),1),1,'first'); % find next integer divisor
  if isemptyfield(Seq.AQSlice(1), 'TurboFactor')  % number of image k-lines per excitation
    Seq.AQSlice(1).TurboFactor = numel(Seq.AQSlice(1).kLineOrder) / Seq.AQSlice(1).TurboBlocks;
    Seq.AQSlice(1).TurboFactorAll = nPhase / Seq.AQSlice(1).TurboBlocks;
  else
    Seq.AQSlice(1).TurboFactorAll = Seq.AQSlice.TurboFactor;
  end
  Seq.AQSlice(1).TurboBlocks = numel(Seq.AQSlice(1).kLineOrder) / Seq.AQSlice(1).TurboFactor;
  Seq.AQSlice(1).TurboBlocksAll = nPhase / Seq.AQSlice(1).TurboFactorAll;

  if (Seq.AQSlice(1).TurboFactor ~= 1 || Seq.AQSlice(1).nImages ~= 1 || Seq.AQSlice(1).oddEvenEchoes) && ...
      Seq.DephaseBefore180
    error('PD:sequence_SpinEcho:TurboPreventsDephaseBefore180', ...
      '"DephaseBefore180" cannot be used if "AQSlice.TurboFactor" or "AQSlice.nImages" differs from 1 or with oddEvenEchoes set to true.');
  end

  if isemptyfield(Seq, 'SteadyState_PreShots90')
    Seq.SteadyState_PreShots90 = 2*((prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS)~=Seq.AQSlice(1).TurboFactor) || Seq.AQSlice(1).phaseCycling || Seq.CorrectPhase);
  end
  if isemptyfield(Seq, 'SteadyState_PostShots90')
    Seq.SteadyState_PostShots90 = 0*((prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS)~=Seq.AQSlice(1).TurboFactor) || Seq.AQSlice(1).phaseCycling) ...
      + (Seq.CorrectPhase > 0);
  end
  if (Seq.CorrectPhase > 0) && Seq.SteadyState_PostShots90 < 1
    warning('PD:sequence_Spin_Echo:CorrectPhaseNoPost', ...
      ['There should be at least one post-shot (Seq.SteadyState_PostShots90 = %d) ', ...
      'when using Seq.CorrectPhase'], ...
      Seq.SteadyState_PostShots90);
  end
  if isemptyfield(Seq, 'SteadyState_PreShots180')
    Seq.SteadyState_PreShots180 = 2*((Seq.AQSlice(1).TurboFactor~=1) || Seq.AQSlice(1).oddEvenEchoes);
  end
  if isemptyfield(Seq, 'SteadyState_PreShots180ReadGrad')
    Seq.SteadyState_PreShots180ReadGrad = 1;
  end
  if isemptyfield(Seq, 'SteadyState_PreShots180SpoilerGrad')
    Seq.SteadyState_PreShots180SpoilerGrad = 1;
  end
  if isemptyfield(Seq, 'SteadyState_PreShots180SpoilerGradFirst')
    Seq.SteadyState_PreShots180SpoilerGradFirst = Seq.SteadyState_PreShots180SpoilerGrad;
  end
  if isemptyfield(Seq, 'SteadyState_PostShots180')
    Seq.SteadyState_PostShots180 = 0*((Seq.AQSlice(1).TurboFactor~=1) || Seq.AQSlice(1).oddEvenEchoes);
  end

  % settings for sequence plot
  if isemptyfield(Seq, 'plotSeqAQ'), Seq.plotSeqAQ = []; end
  if isemptyfield(Seq, 'plotSeqEchoTrain'), Seq.plotSeqEchoTrain = Seq.plotSeqAQ; end
  if ~isempty(Seq.plotSeqEchoTrain)
    Seq.plotSeq = Seq.plotSeqEchoTrain;
    if isemptyfield(Seq,'plotSeqStart')
      Seq.plotSeqStart = 1 + Seq.SteadyState_PreShots90 * ...
        (1 + Seq.SteadyState_PreShots180 + Seq.SteadyState_PostShots180 + ...
        (prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS)*(Seq.AQSlice(1).oddEvenEchoes+1))*Seq.AQSlice(1).nImages/Seq.AQSlice(1).TurboBlocks);
    end
    if isemptyfield(Seq,'plotSeqEnd')
      Seq.plotSeqEnd = Seq.plotSeqStart + ...
        (1 + Seq.SteadyState_PreShots180 + Seq.SteadyState_PostShots180 + ...
        (prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS)*(Seq.AQSlice(1).oddEvenEchoes+1))*Seq.AQSlice(1).nImages/Seq.AQSlice(1).TurboBlocks) * ...
        Seq.AQSlice(1).TurboBlocks * Seq.AQSlice(1).phaseCycleSteps - 1;
    end
    if ~isfield(Seq, 'plotSequence'), Seq.plotSequence = struct(); end
    if isemptyfield(Seq.plotSequence, 'xLim')
      Seq.plotSequence.xLim = [-Inf, ...
        Seq.tEcho * (1 + Seq.SteadyState_PreShots180 + Seq.SteadyState_PostShots180 + ...
        (prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS)*(Seq.AQSlice(1).oddEvenEchoes+1))*Seq.AQSlice(1).nImages/Seq.AQSlice(1).TurboBlocks )];
    end
    if isemptyfield(Seq.plotSequence, 'wraps'), Seq.plotSequence.wraps = Seq.AQSlice(1).TurboBlocks*Seq.AQSlice(1).phaseCycleSteps;end
  end

  % repetition time (of excitations)
  tTurboBlock = ...
    Seq.tEcho * (Seq.SteadyState_PreShots180 + Seq.SteadyState_PostShots180 + ...
    (prod(Seq.AQSlice(1).nPhase)*prod(Seq.AQSlice(1).PhaseOS)*(Seq.AQSlice(1).oddEvenEchoes+1)*Seq.AQSlice(1).nImages) / Seq.AQSlice(1).TurboBlocks) + ...
    max(Seq.tEcho, ...
        Seq.tEcho/2 + ...  % ~ center of acquisition window
        HW.RX(Seq.AQSlice(1).iDevice).ClampCoil.Enable * ...
          (1/2/Seq.AQSlice(1).HzPerPixMin + ...  % ~ end of AQ window (in last tRep of echo train)
           HW.RX(Seq.AQSlice(1).iDevice).ClampCoil.tPostset + 10e-6 + ...  % ~ end of clamp coil signal
           0*HW.RX(Seq.AQSlice(1).iDevice).ClampCoil.tOffset));  %
  if ~isempty(Seq.RepetitionTime) && ~isempty(Seq.AQSlice(1).TurboBreak)
    if Seq.RepetitionTime ~= Seq.AQSlice(1).TurboBreak + tTurboBlock
      error('PD:sequence_Spin_Echo:RepetitionTimeOrTurboBreak', ...
        'Please set only one of the variables Seq.RepetitionTime or Seq.AQSlice(1).TurboBreak.');
    end
  elseif isempty(Seq.RepetitionTime) && isempty(Seq.AQSlice(1).TurboBreak)
    if isempty(Seq.T1Estimated)
      error('PD:sequence_Spin_Echo:T1Estimated', ...
        'Please set one of the variables Seq.T1Estimated, Seq.RepetitionTime or Seq.AQSlice(1).TurboBreak.');
    else
      Seq.AQSlice(1).TurboBreak = Seq.T1Estimated;
    end
  end
  if ~isempty(Seq.AQSlice(1).TurboBreak)
    Seq.RepetitionTime = Seq.AQSlice(1).TurboBreak + tTurboBlock;
  elseif ~isempty(Seq.RepetitionTime)
    Seq.AQSlice(1).TurboBreak = Seq.RepetitionTime - tTurboBlock;
  end
  if Seq.AQSlice(1).TurboBreak < 1e-3
    error('PD:sequence_Spin_Echo:TurboBreakTooShort', ...
      'Seq.AQSlice(1).TurboBreak must be at least 1 ms. (It is %.3f ms.) Increase TurboBreak (or Seq.RepetitionTime).', ...
      Seq.AQSlice(1).TurboBreak*1e3);
  end

  % repetition time (of loops)
  if ~isempty(Seq.LoopsRepetitionTime) && ~isempty(Seq.LoopsBreak)
    % This is an error. We don't know Seq.tOffset(1) yet, so we can't decide
    % whether the combination of values will match.
    error('PD:sequence_Spin_Echo:LoopsRepetitionTimeOrLoopsBreak', ...
      'Please set only one of the variables Seq.LoopsRepetitionTime or Seq.LoopsBreak.');
  elseif isempty(Seq.LoopsRepetitionTime) && isempty(Seq.LoopsBreak)
    % Loops maintain the RepetitionTime between sequences
    Seq.LoopsBreak = Seq.AQSlice(1).TurboBreak;
  elseif ~isempty(Seq.LoopsRepetitionTime)
    % calculate LoopsBreak from Seq.LoopsRepetitionTime. (But not the reverse.)
    Seq.LoopsBreak = Seq.LoopsRepetitionTime - ...
      Seq.RepetitionTime * (Seq.SteadyState_PreShots90 + Seq.SteadyState_PostShots90 + Seq.AQSlice(1).TurboBlocks) - Seq.AQSlice(1).TurboBreak + Seq.tEcho;
  end

  if isemptyfield(Seq, 'Find_Frequency_interval'), Seq.Find_Frequency_interval = HW.FindFrequencySweep.maxTime; end
  if isemptyfield(Seq.AQSlice(1), 'Center2OriginImage'), Seq.AQSlice(1).Center2OriginImage = [0,0,0]; end
  Seq.AQSlice(1).nPhase3D=0;

  if isempty(Seq.AQSlice(1).excitationPulse)
    if Seq.AQSlice(1).thickness>=1
      Seq.AQSlice(1).excitationPulse=@Pulse_Rect;
    else
      Seq.AQSlice(1).excitationPulse=@Pulse_RaisedCos;
    end
  end
end

if Seq.AQSlice(1).sizeRead>=1;Seq.CorrectRemanence=0;end
CRLoops = ones(1, Seq.CorrectRemanence);
Seq.LoopName = repmat({'CRLoop'}, 1, Seq.CorrectRemanence);

if min(Seq.AQSlice(1).thickness,Seq.AQSlice(1).thicknessInversion)>=1;Seq.CorrectSliceRephase=0;end
CSRLoops = repmat([1,2], 1, Seq.CorrectSliceRephase);
Seq.LoopName = [Seq.LoopName, repmat({'CSRLoopPlus', 'CSRLoopMinus'}, 1, Seq.CorrectSliceRephase)];

CPRLoops = repmat([1,2], 1, Seq.CorrectPhaseRephase);
Seq.LoopName = [Seq.LoopName, repmat({'CPRLoopPlus', 'CPRLoopMinus'}, 1, Seq.CorrectPhaseRephase)];

if Seq.AQSlice(1).sizeRead>=1;Seq.CorrectReadRephase=0;end
CRRLoops = ones(1, Seq.CorrectReadRephase);
Seq.LoopName = [Seq.LoopName, repmat({'CRRLoop'}, 1, Seq.CorrectReadRephase)];

CB0Loops = repmat([1,2], 1, Seq.CorrectB0Read.Get);
Seq.LoopName = [Seq.LoopName, repmat({'B0map_tEcho1', 'B0map_tEcho2'}, 1, Seq.CorrectB0Read.Get)];

PreShotsLoops = ones(1, Seq.SteadyState_PreShotsLoops);
Seq.LoopName = [Seq.LoopName, repmat({Seq.SteadyState_PreShotsLoopsName}, 1, Seq.SteadyState_PreShotsLoops)];

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
  % Caution: Using high values over a long time will damage the gradient coils
  % and amplifiers!
  [Grad(1:HW.Grad(Seq.AQSlice(1).iDevice).n).Shim] = deal([]);
  % if there is a one in the array the gradients of the prior TR are used.
  % (reduces data traffic between the PC and the MRI device)
  [Grad(1:HW.Grad(Seq.AQSlice(1).iDevice).n).Repeat] = deal([]);
end

init.AQ = AQ;
init.TX = TX;
init.Seq = Seq;
init.Grad = Grad;
LoopCount = 0;
dataB0 = Seq.dataB0;
for Loop = [CRLoops, CSRLoops, CPRLoops, CRRLoops, CB0Loops, PreShotsLoops, 1:Seq.Loops]
  LoopCount = LoopCount+1;
  AQ = init.AQ;
  TX = init.TX;
  Grad = init.Grad;
  Seq = init.Seq;
  Seq.Loop = Loop;
  Seq.LoopNameCount = LoopCount;

  %% collect offset from previous correction measurements
  % Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset=0;
  % Seq.AQSlice(1).PhaseGradTimeIntegralRephaseOffset=0;
  if exist('SeqLoop', 'var') && isfield(SeqLoop, 'AQSlice') && isfield(SeqLoop, 'data')
    if isfield(SeqLoop.AQSlice, 'SliceGradTimeIntegralRephaseOffset')
      if isfield(SeqLoop.data, 'SliceGradTimeIntegralRephaseOffset')
        Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset=SeqLoop.AQSlice(1).SliceGradTimeIntegralRephaseOffset+SeqLoop.data.SliceGradTimeIntegralRephaseOffset;
      else
        Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset=SeqLoop.AQSlice(1).SliceGradTimeIntegralRephaseOffset;
      end
    end
    if isfield(SeqLoop.AQSlice, 'PhaseGradTimeIntegralRephaseOffset')
      if isfield(SeqLoop.data, 'PhaseGradTimeIntegralRephaseOffset')
        Seq.AQSlice(1).PhaseGradTimeIntegralRephaseOffset=SeqLoop.AQSlice(1).PhaseGradTimeIntegralRephaseOffset+SeqLoop.data.PhaseGradTimeIntegralRephaseOffset;
      else
        Seq.AQSlice(1).PhaseGradTimeIntegralRephaseOffset=SeqLoop.AQSlice(1).PhaseGradTimeIntegralRephaseOffset;
      end
    end
  end

  %% prepare loop type
  switch Seq.LoopName{Seq.LoopNameCount}

    case 'CRLoop'
      % correct remanence in pole shoes
      Seq.AQSlice(1).nPhase(1) = 1;     % number of pixels in phase(1)
      Seq.AQSlice(1).nPhase(2) = 1;     % number of pixels in phase(2)
      Seq.AQSlice(1).nPhase(3) = 1;     % number of pixels in phase(3) for CSI, if nPhase(3)>1 nRead=1 sizeRead=1e12
      Seq.AQSlice(1).PhaseOS(1) = 1;    % oversampling phase(1)  1...
      Seq.AQSlice(1).PhaseOS(2) = 1;    % oversampling phase(2)  1...
      Seq.AQSlice(1).PhaseOS(3) = 1;    % oversampling phase(3)  1... (set to 1 if nPhase(3)=1;
      % no turbo
      Seq.AQSlice(1).TurboFactor = 1;   % number of echoes per excitation
      Seq.AQSlice(1).TurboFactorAll = 1;
      Seq.AQSlice(1).TurboBlocks = 1;
      Seq.AQSlice(1).TurboBlocksAll = 1;
      Seq.AQSlice(1).inversionFlipAngleIncrement = 0;
      Seq.AQSlice(1).inversionFlipAngle = 180;
      Seq.AQSlice(1).excitationFlipAngleIncrement = 0;
      Seq.AQSlice(1).excitationFlipAngle = 90;

      Seq.AQSlice(1).nImages = 1;
      Seq.AQSlice(1).oddEvenEchoes = 0;
      Seq.AQSlice(1).phaseCycleAngles = 0;
      Seq.tReadoutOffset=0;
      Seq.SteadyState_PreShots90 = 0;
      if ~Seq.LoopSeqPlot
        Seq.plotSeqTR = [];     % plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeq = [];       % plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeqStart = [];
        Seq.plotSeqEnd = [];
        Seq.plotSequence = [];
      end

    case 'CSRLoopPlus'
      % correct slice rephase
      Seq.AQSlice(1).SliceGradSign = [1, -1];

      % remove actual read encoder
      % Seq.AQSlice(1).ReadGradTimeIntegralOffset=0;
      if exist('SeqLoop', 'var') && isfield(SeqLoop, 'AQSlice') && isfield(SeqLoop, 'data')
        if ~isfield(SeqLoop.data, 'SliceReadGradTimeIntegralOffset')
          Seq.AQSlice(1).ReadGradTimeIntegralOffset=0;
        elseif isfield(SeqLoop.AQSlice, 'ReadGradTimeIntegralOffset')
          % maintain previously determined correction value
          if isfield(SeqLoop.data, 'SliceReadGradTimeIntegralOffset')
            Seq.AQSlice(1).ReadGradTimeIntegralOffset=SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset+SeqLoop.data.SliceReadGradTimeIntegralOffset;
          else
            Seq.AQSlice(1).ReadGradTimeIntegralOffset=SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset;
          end
        end
      end

      % singular phase step
      Seq.AQSlice(1).nPhase(1) = 1;     % number of pixels in phase(1)
      Seq.AQSlice(1).nPhase(2) = 1;     % number of pixels in phase(2)
      Seq.AQSlice(1).nPhase(3) = 1;     % number of pixels in phase(3) for CSI, if nPhase(3)>1 nRead=1 sizeRead=1e12
      Seq.AQSlice(1).PhaseOS(1) = 1;    % oversampling phase(1)  1...
      Seq.AQSlice(1).PhaseOS(2) = 1;    % oversampling phase(2)  1...
      Seq.AQSlice(1).PhaseOS(3) = 1;    % oversampling phase(3)  1... (set to 1 if nPhase(3)=1;
      if Seq.AQSlice(1).ReadCoordinate ~= Seq.AQSlice(1).SliceCoordinate
        % sizePhaseSpoil = 0 has a special meaning for the read direction.
        % We can't leave the spoiler size like this if we change the read
        % direction. The safest option is to remove that spoiler.
        if Seq.AQSlice(1).sizePhaseSpoil(Seq.AQSlice(1).ReadCoordinate) == 0
          Seq.AQSlice(1).sizePhaseSpoil(Seq.AQSlice(1).ReadCoordinate) = Inf;
        end
        if Seq.AQSlice(1).sizePhaseSpoilEnd(Seq.AQSlice(1).ReadCoordinate) == 0
          Seq.AQSlice(1).sizePhaseSpoilEnd(Seq.AQSlice(1).ReadCoordinate) = Inf;
        end
      end
      % read encoding parallel to slice encoding
      Seq.AQSlice(1).ReadCoordinate = Seq.AQSlice(1).SliceCoordinate;  % direction of read:  x = 1,  y = 2, z = 3
      % no turbo
      Seq.AQSlice(1).TurboFactor = 1;   % number of echoes per excitation
      Seq.AQSlice(1).TurboFactorAll = 1;
      Seq.AQSlice(1).TurboBlocks = 1;
      Seq.AQSlice(1).TurboBlocksAll = 1;
      Seq.AQSlice(1).inversionFlipAngleIncrement = 0;
      Seq.AQSlice(1).excitationFlipAngleIncrement = 0;
      Seq.AQSlice(1).nImages = 1;
      Seq.AQSlice(1).oddEvenEchoes = 0;
      Seq.AQSlice(1).phaseCycleAngles = 0;
      Seq.AQSlice(1).nRead = 16;        % number of pixels in read, if nRead>1 nPhase(1)=1
      Seq.AQSlice(1).ReadOS = max(2, ...
        ceil(HW.RX(Seq.AQSlice(1).iDevice).fSample / ...
          min(8000, HW.RX(Seq.AQSlice(1).iDevice).CIC_Decimation_Max) / ...
          (Seq.AQSlice(1).HzPerPixMin*Seq.AQSlice(1).nRead)));
      Seq.AQSlice(1).sizeRead = min(Seq.AQSlice(1).thickness, Seq.AQSlice(1).thicknessInversion)*2;
      Seq.SteadyState_PreShots90 = 0;
      Seq.tReadoutOffset = 0;
      if ~Seq.LoopSeqPlot
        Seq.plotSeqTR = [];     % plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeq = [];       % plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeqStart = [];
        Seq.plotSeqEnd = [];
        Seq.plotSequence = [];
      end

    case 'CSRLoopMinus'
      % correct slice rephase
      Seq.AQSlice(1).SliceGradSign = [1, -1];
      % copy encoding size from CSRLoopPlus
      Seq.AQSlice(1).ReadGradTimeIntegralOffset = SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset;
      Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset = SeqLoop.AQSlice(1).SliceGradTimeIntegralRephaseOffset;
      % singular phase step
      Seq.AQSlice(1).nPhase(1) = 1;     % number of pixels in phase(1)
      Seq.AQSlice(1).nPhase(2) = 1;     % number of pixels in phase(2)
      Seq.AQSlice(1).nPhase(3) = 1;     % number of pixels in phase(3) for CSI, if nPhase(3)>1 nRead=1 sizeRead=1e12
      Seq.AQSlice(1).PhaseOS(1) = 1;    % oversampling phase(1)  1...
      Seq.AQSlice(1).PhaseOS(2) = 1;    % oversampling phase(2)  1...
      Seq.AQSlice(1).PhaseOS(3) = 1;    % oversampling phase(3)  1... (set to 1 if nPhase(3)=1;
      if Seq.AQSlice(1).ReadCoordinate ~= Seq.AQSlice(1).SliceCoordinate
        % sizePhaseSpoil = 0 has a special meaning for the read direction.
        % We can't leave the spoiler size like this if we change the read
        % direction. The safest option is to remove that spoiler.
        if Seq.AQSlice(1).sizePhaseSpoil(Seq.AQSlice(1).ReadCoordinate) == 0
          Seq.AQSlice(1).sizePhaseSpoil(Seq.AQSlice(1).ReadCoordinate) = Inf;
        end
        if Seq.AQSlice(1).sizePhaseSpoilEnd(Seq.AQSlice(1).ReadCoordinate) == 0
          Seq.AQSlice(1).sizePhaseSpoilEnd(Seq.AQSlice(1).ReadCoordinate) = Inf;
        end
      end
      % read encoding parallel to slice encoding
      Seq.AQSlice(1).ReadCoordinate = Seq.AQSlice(1).SliceCoordinate;  % direction of read:  x = 1,  y = 2, z = 3
      % no turbo
      Seq.AQSlice(1).TurboFactor = 1;   % number of echoes per excitation
      Seq.AQSlice(1).TurboFactorAll = 1;
      Seq.AQSlice(1).TurboBlocks = 1;
      Seq.AQSlice(1).TurboBlocksAll = 1;
      Seq.AQSlice(1).inversionFlipAngleIncrement = 0;
      Seq.AQSlice(1).excitationFlipAngleIncrement = 0;
      Seq.AQSlice(1).nImages = 1;
      Seq.AQSlice(1).oddEvenEchoes = 0;
      Seq.AQSlice(1).phaseCycleAngles = 0;
      Seq.AQSlice(1).nRead = 16;        % number of pixels in read, if nRead>1 nPhase(1)=1
      Seq.AQSlice(1).ReadOS =  max(2, ...
        ceil(HW.RX(Seq.AQSlice(1).iDevice).fSample / ...
          min(8000, HW.RX(Seq.AQSlice(1).iDevice).CIC_Decimation_Max) / ...
          (Seq.AQSlice(1).HzPerPixMin*Seq.AQSlice(1).nRead)));
      Seq.AQSlice(1).sizeRead=min(Seq.AQSlice(1).thickness,Seq.AQSlice(1).thicknessInversion)*2;
      Seq.SteadyState_PreShots90 = 0;
      Seq.tReadoutOffset = 0;
      if ~Seq.LoopSeqPlot
        Seq.plotSeqTR = [];     % plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeq = [];       % plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeqStart = [];
        Seq.plotSeqEnd = [];
        Seq.plotSequence = [];
      end

    case 'CPRLoopPlus'
      % correct phase rephase
      Seq.AQSlice(1).PhaseGradSign = [1, -1];

      % Use the same encoding size for read as will be used for the phase
      % encoding.
      % Seq.AQSlice(1).ReadGradTimeIntegralOffset=0;
      if exist('SeqLoop', 'var') && isfield(SeqLoop, 'AQSlice') && isfield(SeqLoop, 'data')
        if ~isfield(SeqLoop.data, 'PhaseReadGradTimeIntegralOffset')
          Seq.AQSlice(1).ReadGradTimeIntegralOffset = 0;
        elseif isfield(SeqLoop.AQSlice, 'ReadGradTimeIntegralOffset')
          if isfield(SeqLoop.data, 'PhaseReadGradTimeIntegralOffset')
            Seq.AQSlice(1).ReadGradTimeIntegralOffset = ...
              SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset + SeqLoop.data.PhaseReadGradTimeIntegralOffset;
          else
            Seq.AQSlice(1).ReadGradTimeIntegralOffset = ...
              SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset;
          end
        end
      end

      % singular phase step
      Seq.AQSlice(1).nPhase(1) = 1;     % number of pixels in phase(1)
      Seq.AQSlice(1).nPhase(2) = 1;     % number of pixels in phase(2)
      Seq.AQSlice(1).nPhase(3) = 1;     % number of pixels in phase(3) for CSI, if nPhase(3)>1 nRead=1 sizeRead=1e12
      Seq.AQSlice(1).PhaseOS(1) = 1;    % oversampling phase(1)  1...
      Seq.AQSlice(1).PhaseOS(2) = 1;    % oversampling phase(2)  1...
      Seq.AQSlice(1).PhaseOS(3) = 1;    % oversampling phase(3)  1... (set to 1 if nPhase(3)=1;
      if Seq.AQSlice(1).ReadCoordinate ~= Seq.AQSlice(1).PhaseCoordinate(2)
        % sizePhaseSpoil = 0 has a special meaning for the read direction.
        % We can't leave the spoiler size like this if we change the read
        % direction. The safest option is to remove that spoiler.
        if Seq.AQSlice(1).sizePhaseSpoil(Seq.AQSlice(1).ReadCoordinate) == 0
          Seq.AQSlice(1).sizePhaseSpoil(Seq.AQSlice(1).ReadCoordinate) = Inf;
        end
        if Seq.AQSlice(1).sizePhaseSpoilEnd(Seq.AQSlice(1).ReadCoordinate) == 0
          Seq.AQSlice(1).sizePhaseSpoilEnd(Seq.AQSlice(1).ReadCoordinate) = Inf;
        end
      end
      % read encoding parallel to phase(2) encoding
      Seq.AQSlice(1).ReadCoordinate = Seq.AQSlice(1).PhaseCoordinate(2);  % direction of read:  x = 1,  y = 2, z = 3
      % no turbo
      Seq.AQSlice(1).TurboFactor = 1;   % number of echoes per excitation
      Seq.AQSlice(1).TurboFactorAll = 1;
      Seq.AQSlice(1).TurboBlocks = 1;
      Seq.AQSlice(1).TurboBlocksAll = 1;
      Seq.AQSlice(1).inversionFlipAngleIncrement = 0;
      Seq.AQSlice(1).excitationFlipAngleIncrement = 0;
      Seq.AQSlice(1).nImages = 1;
      Seq.AQSlice(1).oddEvenEchoes = 0;
      Seq.AQSlice(1).phaseCycleAngles = 0;
      Seq.SteadyState_PreShots90 = 0;
      Seq.tReadoutOffset = 0;
      if ~Seq.LoopSeqPlot
        Seq.plotSeqTR = [];     % plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeq = [];       % plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeqStart = [];
        Seq.plotSeqEnd = [];
        Seq.plotSequence = [];
      end

    case 'CPRLoopMinus'
      % correct phase rephase
      Seq.AQSlice(1).PhaseGradSign = [1, -1];
      % copy encoding size from CPRLoopPlus
      Seq.AQSlice(1).ReadGradTimeIntegralOffset = SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset;
      Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset = SeqLoop.AQSlice(1).SliceGradTimeIntegralRephaseOffset;
      % singular phase step
      Seq.AQSlice(1).nPhase(1) = 1;     % number of pixels in phase(1)
      Seq.AQSlice(1).nPhase(2) = 1;     % number of pixels in phase(2)
      Seq.AQSlice(1).nPhase(3) = 1;     % number of pixels in phase(3) for CSI, if nPhase(3)>1 nRead=1 sizeRead=1e12
      Seq.AQSlice(1).PhaseOS(1) = 1;    % oversampling phase(1)  1...
      Seq.AQSlice(1).PhaseOS(2) = 1;    % oversampling phase(2)  1...
      Seq.AQSlice(1).PhaseOS(3) = 1;    % oversampling phase(3)  1... (set to 1 if nPhase(3)=1;
      if Seq.AQSlice(1).ReadCoordinate ~= Seq.AQSlice(1).PhaseCoordinate(2)
        % sizePhaseSpoil = 0 has a special meaning for the read direction.
        % We can't leave the spoiler size like this if we change the read
        % direction. The safest option is to remove that spoiler.
        if Seq.AQSlice(1).sizePhaseSpoil(Seq.AQSlice(1).ReadCoordinate) == 0
          Seq.AQSlice(1).sizePhaseSpoil(Seq.AQSlice(1).ReadCoordinate) = Inf;
        end
        if Seq.AQSlice(1).sizePhaseSpoilEnd(Seq.AQSlice(1).ReadCoordinate) == 0
          Seq.AQSlice(1).sizePhaseSpoilEnd(Seq.AQSlice(1).ReadCoordinate) = Inf;
        end
      end
      % read encoding parallel to phase(2) encoding
      Seq.AQSlice(1).ReadCoordinate = Seq.AQSlice(1).PhaseCoordinate(2);  % direction of read:  x = 1,  y = 2, z = 3
      % no turbo
      Seq.AQSlice(1).TurboFactor = 1;   % number of echoes per excitation
      Seq.AQSlice(1).TurboFactorAll = 1;
      Seq.AQSlice(1).TurboBlocks = 1;
      Seq.AQSlice(1).TurboBlocksAll = 1;
      Seq.AQSlice(1).inversionFlipAngleIncrement = 0;
      Seq.AQSlice(1).excitationFlipAngleIncrement = 0;
      Seq.AQSlice(1).nImages = 1;
      Seq.AQSlice(1).oddEvenEchoes = 0;
      Seq.AQSlice(1).phaseCycleAngles = 0;
      Seq.SteadyState_PreShots90 = 0;
      Seq.tReadoutOffset = 0;
      if ~Seq.LoopSeqPlot
        Seq.plotSeqTR = [];     % plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeq = [];       % plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeqStart = [];
        Seq.plotSeqEnd = [];
        Seq.plotSequence = [];
      end

    case 'CRRLoop'
      % correct read rephase
      % Seq.AQSlice(1).ReadGradTimeIntegralOffset=0;
      if exist('SeqLoop', 'var') && isfield(SeqLoop, 'AQSlice') && isfield(SeqLoop, 'data')
        if isfield(SeqLoop.data, 'SliceReadGradTimeIntegralOffset') || ...
            isfield(SeqLoop.data, 'PhaseReadGradTimeIntegralOffset')
          Seq.AQSlice(1).ReadGradTimeIntegralOffset=0;
        elseif isfield(SeqLoop.AQSlice, 'ReadGradTimeIntegralOffset')
          if isfield(SeqLoop.data, 'ReadGradTimeIntegralOffset')
            Seq.AQSlice(1).ReadGradTimeIntegralOffset = ...
              SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset + SeqLoop.data.ReadGradTimeIntegralOffset;
          else
            Seq.AQSlice(1).ReadGradTimeIntegralOffset = ...
              SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset;
          end
        end
      end

      Seq.AQSlice(1).nPhase(1) = 1;     % number of pixels in phase(1)
      Seq.AQSlice(1).nPhase(2) = 1;     % number of pixels in phase(2)
      Seq.AQSlice(1).nPhase(3) = 1;     % number of pixels in phase(3) for CSI, if nPhase(3)>1 nRead=1 sizeRead=1e12
      Seq.AQSlice(1).PhaseOS(1) = 1;    % oversampling phase(1)  1...
      Seq.AQSlice(1).PhaseOS(2) = 1;    % oversampling phase(2)  1...
      Seq.AQSlice(1).PhaseOS(3) = 1;    % oversampling phase(3)  1... (set to 1 if nPhase(3)=1;
      % no turbo
      Seq.AQSlice(1).TurboFactor = 1;   % number of echoes per excitation
      Seq.AQSlice(1).TurboFactorAll = 1;
      Seq.AQSlice(1).TurboBlocks = 1;
      Seq.AQSlice(1).TurboBlocksAll = 1;
      Seq.AQSlice(1).inversionFlipAngleIncrement = 0;
      Seq.AQSlice(1).excitationFlipAngleIncrement = 0;
      Seq.AQSlice(1).nImages = 1;
      Seq.AQSlice(1).oddEvenEchoes = 0;
      Seq.AQSlice(1).phaseCycleAngles = 0;
      Seq.SteadyState_PreShots90 = 0;
      Seq.tReadoutOffset = 0;
      if ~Seq.LoopSeqPlot
        Seq.plotSeqTR = [];     % plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeq = [];       % plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
        Seq.plotSeqStart = [];
        Seq.plotSeqEnd = [];
        Seq.plotSequence = [];
      end

    case 'B0map_tEcho1'
      if Seq.AQSlice(1).nRead < 2
        continue;
      end
      if Seq.CorrectB0Read.tReadoutDiff == 0
        warning('PD:sequence_Spin_Echo:B0ZeroReadoutDiff', ...
          'Seq.CorrectB0Read.tReadoutDiff must not equal 0 for B0map_tEcho. Skipping B0 measurement...')
        continue;
      end
      Seq.tReadoutOffset = Seq.tReadoutOffset - Seq.CorrectB0Read.tReadoutDiff/2;
      Seq.Find_Frequency_interval = 0;
      Seq.AQSlice(1).ZeroFillWindowSize = Seq.CorrectB0Read.ZeroFillWindowSize;
      Seq.CorrectB0Read.Use = false;

    case 'B0map_tEcho2'
      if Seq.AQSlice(1).nRead < 2 || Seq.CorrectB0Read.tReadoutDiff == 0
        continue;
      end
      Seq.tReadoutOffset = Seq.tReadoutOffset + Seq.CorrectB0Read.tReadoutDiff/2;
      Seq.Find_Frequency_interval = Inf;  % FIXME: Find frequency again for second measurement?
      Seq.AQSlice(1).ZeroFillWindowSize = Seq.CorrectB0Read.ZeroFillWindowSize;
      Seq.CorrectB0Read.Use = false;

    case 'normal'
      % Seq.AQSlice(1).ReadGradTimeIntegralOffset=0;
      if exist('SeqLoop', 'var') && isfield(SeqLoop, 'AQSlice') && isfield(SeqLoop, 'data')
        if isfield(SeqLoop.data, 'SliceReadGradTimeIntegralOffset') || ...
            isfield(SeqLoop.data, 'PhaseReadGradTimeIntegralOffset')

        else
          if isfield(SeqLoop.AQSlice, 'ReadGradTimeIntegralOffset')
            % copy (cumulative) correction values from previous loops
            if isfield(SeqLoop.data, 'ReadGradTimeIntegralOffset')
              Seq.AQSlice(1).ReadGradTimeIntegralOffset = ...
                SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset + SeqLoop.data.ReadGradTimeIntegralOffset;
            else
              Seq.AQSlice(1).ReadGradTimeIntegralOffset = ...
                SeqLoop.AQSlice(1).ReadGradTimeIntegralOffset;
            end
          end
        end
      end

    otherwise
      % do nothing

  end

  %% Larmor frequency
  if Seq.StartSequence && (Seq.LoopNameCount==1 || Seq.LoopsBreakExactly==0)
    if isemptyfield(mySave, 'lastTime'), mySave.lastTime = 0; end
    if isinf(Seq.Find_Frequency_interval) || ...
        (now*24*3600-mySave.lastTime < Seq.Find_Frequency_interval-0.01)
      % Use last results of frequency sweep (i.e. mySave.HW.B0)
      [HW, mySave] = Find_Frequency_Sweep(HW, mySave, Inf);
    else
      if ~isempty(Seq.StartSequenceTime)
        % delay execution for timed loops
        sleep(max(0, Seq.StartSequenceTime - now*24*3600));
      end

      oldFindFrequencyPause = HW.FindFrequencyPause;
      hwGuard = onCleanup(@() setfield(HW, 'FindFrequencyPause', oldFindFrequencyPause)); %#ok<SFLD>
      HW.FindFrequencyPause = 0;
      [HW, mySave] = Find_Frequency_Sweep(HW, mySave, 0);  % Find magnet frequency
      delete(hwGuard);

      % Update start time for sequence to account for the magnetization lost due
      % to the frequency sweep
      Seq.StartSequenceTime = now*24*3600 + max([Seq.LoopsBreak, HW.FindFrequencyPause]);
    end
  end

  if Seq.CorrectPhase && Seq.CorrectPhase_Set_fLarmorOfNextLoop && Seq.StartSequence ...
      && Seq.LoopNameCount>1 && strcmp(Seq.LoopName{Seq.LoopNameCount-1}, 'normal')  % second normal Loop
    % set fLarmor and B0 of HW
    HW.fLarmor = SeqLoop.dataLoop(Loop-1).fCenter - SeqLoop.dataLoop(Loop-1).CorrectPhase.fOffset(end);
    mySave.HW.B0 = HW.B0;  % save B0 to mySave
    mySave.HW.fLarmor = HW.fLarmor;  % save fLarmor
  end

  if strcmp(Seq.LoopName{Seq.LoopNameCount}, 'normal') && Loop==1
    % Seq.plotSeqTR = 1:3;  % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
    % Seq.plotSeq = 1:3;  % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
  else
    if ~Seq.LoopSeqPlot
      Seq.plotSeqTR = [];  % Plot sequence, all tReps are starting at origin, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
      Seq.plotSeq = [];  % Plot sequence on real timeline, plot RF, AQ and Grad (1==x, 2==y, 3==z, 0 no Grad)
    end
  end

  Seq.AQSlice(1).phaseCycleSteps = numel(Seq.AQSlice(1).phaseCycleAngles);

  %% construct pulse sequence
  % Programming notes:
  % The excitation pulses with slice selection are in a separate tRep.
  % The refocusing pulses (with inversion slice selection) share the same tRep
  % with the following acquisition window including read and phase encoding.
  % The duration of the last tRep in each echo train is adjusted to achieve the
  % set RepetitionTime.

  % number of images acquired per echo train
  numImages = Seq.AQSlice(1).nImages * (Seq.AQSlice(1).oddEvenEchoes+1) * Seq.AQSlice(1).phaseCycleSteps;
  % number of k-lines per image
  numKLines = prod(Seq.AQSlice(1).nPhase) * prod(Seq.AQSlice(1).PhaseOS);
  % assign number to each k-line that will be acquired in this experiment:
  % nImages x odd-even x Turbo-Block x echo train x phase-cycling
  kLinesDim = reshape(1:numImages*numKLines, ...
    Seq.AQSlice(1).nImages, (Seq.AQSlice(1).oddEvenEchoes+1), ...
    Seq.AQSlice(1).TurboBlocksAll, Seq.AQSlice(1).TurboFactorAll, Seq.AQSlice(1).phaseCycleSteps);
  if any(strcmp(Seq.LoopName{Seq.LoopNameCount}, {'normal', 'B0map_tEcho1', 'B0map_tEcho2'}))
    % reshape such that Turbo-Blocks and echo trains share the same dimension
    kLinesDimTemp = reshape(kLinesDim, ...
      Seq.AQSlice(1).nImages, (Seq.AQSlice(1).oddEvenEchoes+1), ...
      numKLines, Seq.AQSlice(1).phaseCycleSteps);
    % permute the order of the k-lines and reshape to the original dimensions
    kLinesDimTemp = kLinesDimTemp(:,:,[Seq.AQSlice(1).kLineOrder, setdiff(1:numKLines, Seq.AQSlice(1).kLineOrder)],:);
    kLinesDimTemp = kLinesDimTemp(:,:,1:numel(Seq.AQSlice(1).kLineOrder),:);
    % consider reduced size if some k-lines are skipped
    szRed = [size(kLinesDim, 1), size(kLinesDim, 2), size(kLinesDim, 3), size(kLinesDim, 4), size(kLinesDim, 5)];
    szRed(4) = szRed(4) / (numel(kLinesDim)/numel(kLinesDimTemp));
    if szRed(4) < 1
      szRed(3) = szRed(3) * szRed(4);
      szRed(4) = 1;
    end
    if any(mod(szRed,1) ~= 0)
      error('PD:sequence_Echo_Standard:incompatibleReduction', ...
        'Number of acquired k-lines incompatible with number of turbo blocks.');
    end
    kLinesDimTemp = reshape(kLinesDimTemp, szRed);
  else
    kLinesDimTemp = kLinesDim;
    % k-line order only applies for imaging loops
    Seq.AQSlice(1).kLineOrder = 1:numel(kLinesDimTemp);
  end

  % order in which k-lines are to be acquired
  % k-lines in echo train x echo trains
  Seq.kLines = reshape(permute(kLinesDimTemp, [2 4 1 5 3]), Seq.AQSlice(1).TurboFactor*(Seq.AQSlice(1).oddEvenEchoes+1)*Seq.AQSlice(1).nImages, []);
  % order of those k-lines for image reconstruction:
  % kLines x nImages x corrections
  Seq.kLinesImages = reshape(permute(kLinesDim, [3 4 1 2 5]), numKLines, Seq.AQSlice(1).nImages, []);
  clear kLinesDimTemp kLinesDim kLinesDim numImages
  % tReps with an excitation pulse
  Seq.P90tReps = cumsum(Seq.SteadyState_PreShots180 + Seq.SteadyState_PostShots180 + ...
    ones(1, size(Seq.kLines, 2) + Seq.SteadyState_PreShots90 + Seq.SteadyState_PostShots90) + ...
    size(Seq.kLines, 1)) - ...
    size(Seq.kLines, 1) - Seq.SteadyState_PostShots180 - Seq.SteadyState_PreShots180;
  % tReps with refocusing pre-shots (i.e. without actual acquisition at the start of the echo train)
  Seq.PreShotstReps = repmat(Seq.P90tReps, Seq.SteadyState_PreShots180, 1) + ...
    repmat((1:Seq.SteadyState_PreShots180).', 1, size(Seq.P90tReps,2));
  % tReps with refocusing pulses
  Seq.P180tReps = repmat(Seq.P90tReps, size(Seq.kLines,1), 1) + ...
    repmat((1:size(Seq.kLines,1)).', 1, size(Seq.P90tReps,2)) + Seq.SteadyState_PreShots180;
  % tReps which contain the k-lines in above ordering
  Seq.kLinestReps = repmat(Seq.P90tReps(1,Seq.SteadyState_PreShots90+1:end-Seq.SteadyState_PostShots90), size(Seq.kLines,1), 1) + ...
    repmat((1:size(Seq.kLines,1)).', 1, size(Seq.kLines,2)) + Seq.SteadyState_PreShots180;
  % tReps with inversion pulses without acquisition
  Seq.P180tRepsPrePost = setdiff(Seq.P180tReps(:), Seq.kLinestReps(:));
  % tReps with refocusing pre-shots (i.e. without actual acquisition at the end of the echo train)
  Seq.PostShotstReps = repmat(Seq.P90tReps, Seq.SteadyState_PostShots180, 1) + ...
    repmat((1:Seq.SteadyState_PostShots180).', 1, size(Seq.P90tReps,2)) + Seq.SteadyState_PreShots180 + size(Seq.kLines,1);
  % Matrix with the tReps grouped into turbo blocks:
  % tReps in Turbo train x Turbo blocks
  Seq.tRepTurboBlock = [Seq.P90tReps; Seq.PreShotstReps; Seq.P180tReps; Seq.PostShotstReps];

  % cycle phase of inversion pulses between blocks by +/-90 degrees
  Seq.tRepPhase180 = zeros(Seq.AQSlice(1).phaseCycleSteps, size(Seq.kLines, 2)/Seq.AQSlice(1).phaseCycleSteps) + 90;
  Seq.tRepPhase180(:,1:2:end) = -90;
  % additionally cycle phase of inversion pulses in requested steps
  Seq.tRepPhase180 = mod(bsxfun(@plus, Seq.tRepPhase180, Seq.AQSlice(1).phaseCycleAngles), 360);
  Seq.tRepPhase180 = Seq.tRepPhase180(:).';
  % cycle phase of excitation pulses between blocks by +/- Seq.AQSlice(1).excitationPhase
  % FIXME: Correction of phase error of excitation pulse?
  Seq.tRepPhase90 = zeros(Seq.AQSlice(1).phaseCycleSteps, size(Seq.kLines, 2)/Seq.AQSlice(1).phaseCycleSteps) - Seq.AQSlice(1).excitationPhase;
  Seq.tRepPhase90(:,1:2:end) = Seq.AQSlice(1).excitationPhase;
  % repeat for all phase cycle steps
  Seq.tRepPhase90 = mod(repmat(Seq.tRepPhase90, Seq.AQSlice(1).phaseCycleSteps/size(Seq.tRepPhase90,1), 1), 360);
  Seq.tRepPhase90 = Seq.tRepPhase90(:).';

  % create (more than) enough blocks for pre- and post-shots
  numPrePostMax = max([Seq.SteadyState_PreShots90, Seq.SteadyState_PostShots90]);
  % inversion pulse pre-/post-shots
  tRepPhase180PreShot = zeros(Seq.AQSlice(1).phaseCycleSteps, 2*numPrePostMax) + 90;
  tRepPhase180PreShot(:,1:2:end) = -90;
  tRepPhase180PreShot = mod(bsxfun(@plus, tRepPhase180PreShot, Seq.AQSlice(1).phaseCycleAngles), 360);
  tRepPhase180PreShot = tRepPhase180PreShot(:).';
  % excitation pulse pre-/post-shots
  tRepPhase90PreShot = zeros(Seq.AQSlice(1).phaseCycleSteps, 2*numPrePostMax) - Seq.AQSlice(1).excitationPhase;
  tRepPhase90PreShot(:,1:2:end) = Seq.AQSlice(1).excitationPhase;
  tRepPhase90PreShot = mod(repmat(tRepPhase90PreShot, Seq.AQSlice(1).phaseCycleSteps/size(Seq.tRepPhase90,1), 1), 360);
  tRepPhase90PreShot = tRepPhase90PreShot(:).';

  % combine pre- and post-shots with actual measurement
  phaseCycleBlock = [tRepPhase180PreShot(end-Seq.SteadyState_PreShots90+1:end), Seq.tRepPhase180, tRepPhase180PreShot(1:Seq.SteadyState_PostShots90)];
  phase90CycleBlock = [tRepPhase90PreShot(end-Seq.SteadyState_PreShots90+1:end), Seq.tRepPhase90, tRepPhase90PreShot(1:Seq.SteadyState_PostShots90)];
  if Seq.AQSlice(1).readoutPhaseInversion
    % read out with phase of inversion pulse
    Seq.tRepTurboBlockPhaseRead = repmat(phaseCycleBlock, size(Seq.tRepTurboBlock,1), 1);
  else
    % read out with phase of excitation pulse flipped by inversion pulse
    Seq.tRepTurboBlockPhaseRead = ...
      repmat(mod(phase90CycleBlock - 2*(phaseCycleBlock - phase90CycleBlock), 360), ...
      size(Seq.tRepTurboBlock,1), 1);
  end
  % phase increment for readout in each turbo block
  Seq.tRepTurboBlockPhaseRead(2:end,:) = ...
    mod(bsxfun(@plus, Seq.tRepTurboBlockPhaseRead(2:end,:), ...
    (0:size(Seq.tRepTurboBlock,1)-2).' * Seq.AQSlice(1).readOutPhaseIncrement), 360);
  % phase increment for inversion pulses in each turbo block
  Seq.AQSlice(1).inversionPhaseIncrement=repmat(Seq.AQSlice(1).inversionPhaseIncrement,...
      ceil((size(Seq.tRepTurboBlock,1)-1)./size(Seq.AQSlice(1).inversionPhaseIncrement,1)),1);
  Seq.AQSlice(1).inversionPhaseIncrement(size(Seq.tRepTurboBlock,1):end,:)=[];
  Seq.tRepTurboBlockPhaseInversion = [phase90CycleBlock;  ... % phase of excitation pulses
    mod(bsxfun(@plus, phaseCycleBlock, ...
        [0;cumsum(Seq.AQSlice(1).inversionPhaseIncrement(1:end-1,1))]), 360)];  % phase of inversion pulses (0:size(Seq.tRepTurboBlock,1)-2).' .* Seq.AQSlice(1).inversionPhaseIncrement), 360)];  % phase of inversion pulses
  Seq.tRepTurboBlockPhaseExcitation = Seq.tRepTurboBlockPhaseInversion;  % only first line is ever used

  % Phase cycle p90 and AQ instead of p180
  Seq.tRepTurboBlockPhaseRead = mod(bsxfun(@minus, Seq.tRepTurboBlockPhaseRead, phaseCycleBlock), 360);
  Seq.tRepTurboBlockPhaseExcitation = mod(bsxfun(@minus, Seq.tRepTurboBlockPhaseExcitation, phaseCycleBlock), 360);
  Seq.tRepTurboBlockPhaseInversion = mod(bsxfun(@minus, Seq.tRepTurboBlockPhaseInversion, phaseCycleBlock), 360);

  % Generate tRep times
  Seq.tRep = Seq.tEcho*ones(1,numel(Seq.P90tReps)+numel(Seq.PreShotstReps)+numel(Seq.P180tReps)+numel(Seq.PostShotstReps));
  Seq.tRep(Seq.P90tReps(2:end)-1) = Seq.tEcho + Seq.AQSlice(1).TurboBreak;
  if isemptyfield(Seq, 'tRepEnd')
    Seq.tRepEnd = Seq.tEcho + max(Seq.tEcho, HW.RX(Seq.AQSlice(1).iDevice).ClampCoil.Enable * (1/2/Seq.AQSlice(1).HzPerPixMin + HW.RX(Seq.AQSlice(1).iDevice).ClampCoil.tPostset + 1e-6));
  end
  Seq.tRep(end) = Seq.tRepEnd;

  % assign order of tReps (as acquired) for image reconstruction
  [~, b] = sort(Seq.kLines(:));
  Seq.AQSlice(1).UsetRep = Seq.kLinestReps(b);

  % optionally cycle gradient sign between loops
  Seq.AQSlice(1).GradSign = Seq.AQSlice(1).GradSign(mod(Loop-1, numel(Seq.AQSlice(1).GradSign))+1);
  Seq.AQSlice(1).ReadGradSign = Seq.AQSlice(1).ReadGradSign(mod(Loop-1, numel(Seq.AQSlice(1).ReadGradSign))+1);
  Seq.AQSlice(1).SliceGradSign = Seq.AQSlice(1).SliceGradSign(mod(Loop-1, numel(Seq.AQSlice(1).SliceGradSign))+1);
  Seq.AQSlice(1).PhaseGradSign = Seq.AQSlice(1).PhaseGradSign(mod(Loop-1, numel(Seq.AQSlice(1).PhaseGradSign))+1);


  %% set up structures for get_ReadParameter (i.e. acquisition windows with read encoding)
  % Readout at tEcho
  Seq.Read(1).HzPerPixelMin = Seq.AQSlice(1).HzPerPixMin;
  Seq.Read(1).CenterOfReadout = Seq.tEcho/2 + Seq.tReadoutOffset;
  Seq.Read(1).UseCoordinate = Seq.AQSlice(1).ReadCoordinate;
  Seq.Read(1).GradTimeDelay = Seq.AQSlice(1).ReadTimeDelay;
  Seq.Read(1).GradTimeIntegralOffset = Seq.AQSlice(1).ReadGradTimeIntegralOffset;
  Seq.Read(1).UseAtRepetitionTime = Seq.kLinestReps;
  Seq.Read(1).UseAtRepetitionTimeDephase = Seq.Read(1).UseAtRepetitionTime;
  % optionally, move dephase pulse of first echo in train to front
  Seq.Read(1).UseAtRepetitionTimeDephase(1,:) = ...
    Seq.Read(1).UseAtRepetitionTimeDephase(1,:) - Seq.DephaseBefore180*(1+Seq.SteadyState_PreShots180);
  Seq.Read(1).Phase = Seq.tRepTurboBlockPhaseRead(Seq.Read(1).UseAtRepetitionTime) + Seq.AQSlice(1).readOutPhase;
  Seq.Read(1).PhaseIncrement = 0;
  Seq.Read(1).distance = Seq.AQSlice(1).Center2OriginImage(3);  % FIXME: Should this use Seq.AQSlice(1).ReadCoordinate instead?
  Seq.Read(1).GradSign = Seq.AQSlice(1).ReadGradSign * Seq.AQSlice(1).GradSign;
  % If the readout gradients are omitted in the 180 degrees pre-shots, toggle the sign
  % of the dephase gradient according to the number of inversion pulses until
  % the first acquired echo for Seq.DephaseBefore180.
  Seq.Read(1).GradDephaseSign = -Seq.AQSlice(1).ReadGradSign * Seq.AQSlice(1).GradSign * ...
    (-1)^(Seq.DephaseBefore180*(1+Seq.SteadyState_PreShots180*(~Seq.SteadyState_PreShots180ReadGrad)));
  Seq.Read(1).GradRephaseSign = -Seq.AQSlice(1).ReadGradSign * Seq.AQSlice(1).GradSign;
  Seq.Read(1).GradDephaseLengthFactor = Seq.AQSlice(1).DephaseLengthFactor;
  Seq.Read(1).GradRephaseLengthFactor = Seq.AQSlice(1).RephaseLengthFactor;
  if Seq.AQSlice(1).dualNuclearImage
    Seq.Read(1).GammaX = Seq.AQSlice(1).GammaX;
  end

  % dummy Readouts at tEcho of 180 degrees pre and post shots
  Seq.Read(2).HzPerPixelMin = Seq.Read(1).HzPerPixelMin;
  Seq.Read(2).CenterOfReadout = Seq.Read(1).CenterOfReadout;
  Seq.Read(2).UseCoordinate = Seq.Read(1).UseCoordinate;
  Seq.Read(2).GradTimeDelay = Seq.Read(1).GradTimeDelay;
  Seq.Read(2).GradTimeIntegralOffset = Seq.Read(1).GradTimeIntegralOffset;
  Seq.Read(2).UseAtRepetitionTime = sort([Seq.PreShotstReps(:);Seq.PostShotstReps(:)]);
  Seq.Read(2).UseAtRepetitionTimeDephase = sort([Seq.PreShotstReps(:);Seq.PostShotstReps(:)]);
  Seq.Read(2).Phase = 0;
  Seq.Read(2).PhaseIncrement = 0;
  Seq.Read(2).distance = Seq.Read(1).distance;
  Seq.Read(2).GradSign = Seq.Read(1).GradSign;
  Seq.Read(2).GradDephaseSign = -Seq.AQSlice(1).ReadGradSign * Seq.AQSlice(1).GradSign ;
  Seq.Read(2).GradRephaseSign = Seq.Read(1).GradRephaseSign;
  Seq.Read(2).GradDephaseLengthFactor = Seq.Read(1).GradDephaseLengthFactor;
  Seq.Read(2).GradRephaseLengthFactor = Seq.Read(1).GradRephaseLengthFactor;
  if Seq.AQSlice(1).dualNuclearImage
    Seq.Read(2).GammaX = Seq.AQSlice(1).GammaX;
  end

  % dummy Readouts at 90 degrees pre and post shots
  Seq.Read(3).HzPerPixelMin = Seq.Read(1).HzPerPixelMin;
  Seq.Read(3).CenterOfReadout = Seq.Read(1).CenterOfReadout;
  Seq.Read(3).UseCoordinate = Seq.Read(1).UseCoordinate;
  Seq.Read(3).GradTimeDelay = Seq.Read(1).GradTimeDelay;
  Seq.Read(3).GradTimeIntegralOffset = Seq.Read(1).GradTimeIntegralOffset;
  % Seq.Read(3).UseAtRepetitionTime = setxor(Seq.AQSlice(1).UsetRep(:), Seq.P180tReps(:));
  Seq.Read(3).UseAtRepetitionTime = Seq.P180tRepsPrePost(:);
  Seq.Read(3).UseAtRepetitionTimeDephase = Seq.P180tRepsPrePost(:);
  Seq.Read(3).Phase = 0;
  Seq.Read(3).PhaseIncrement = 0;
  Seq.Read(3).distance = Seq.Read(1).distance;
  Seq.Read(3).GradSign = Seq.Read(1).GradSign;
  Seq.Read(3).GradDephaseSign = Seq.Read(1).GradDephaseSign;
  Seq.Read(3).GradRephaseSign = Seq.Read(1).GradRephaseSign;
  Seq.Read(3).GradDephaseLengthFactor = Seq.Read(1).GradDephaseLengthFactor;
  Seq.Read(3).GradRephaseLengthFactor = Seq.Read(1).GradRephaseLengthFactor;
  if Seq.AQSlice(1).dualNuclearImage
    Seq.Read(3).GammaX = Seq.AQSlice(1).GammaX;
  end


  % generate readout timings
  Seq = get_ReadParameter(Seq, HW);

  Seq.AQSlice(1).AcquisitionTime = Seq.Read(1).AcquisitionTime;
  Seq.AQSlice(1).AcquisitionFrequency = Seq.Read(1).AcquisitionFrequency;
  Seq.Read(2).GradTimeDelayOffset = Seq.Read(1).GradTimeDelayOffset;
  Seq.Read(3).GradTimeDelayOffset = Seq.Read(1).GradTimeDelayOffset;

  % multiplicate the length of the spoil at end
  Spoil.CenterOfDephaseEnd = Seq.Read(1).CenterOfRephase + ...
    (Seq.Read(1).GradRephaseLength-Seq.Read(1).tRamp) .* (Seq.AQSlice(1).SpoilDephaseLengthFactorEnd-1) / 2 + ...
    Seq.AQSlice(1).SpoilDephaseTimeOffset;
  Spoil.GradDephaseLengthEnd = (Seq.Read(1).GradRephaseLength-Seq.Read(1).tRamp) .* Seq.AQSlice(1).SpoilDephaseLengthFactorEnd + Seq.Read(1).tRamp;
  if Seq.AQSlice(1).SpoilDephaseLengthFactorEnd == 0
    % set minimum gradient length
    Spoil.GradDephaseLengthEnd = 2*Seq.Read(1).tRamp + 2/HW.RX(Seq.AQSlice(1).iDevice).fSample;
  end
   % multiplicate the length of the spoil and rephase pulse to achieve shorter tEcho with reduced grad amplitudes
  Spoil.CenterOfDephase = Seq.Read(1).CenterOfRephase + ...
    (Seq.Read(1).GradRephaseLength-Seq.Read(1).tRamp) .* (Seq.AQSlice(1).SpoilDephaseLengthFactor-1) / 2 + ...
    Seq.AQSlice(1).SpoilDephaseTimeOffset;
  Spoil.GradDephaseLength = (Seq.Read(1).GradRephaseLength-Seq.Read(1).tRamp) .* Seq.AQSlice(1).SpoilDephaseLengthFactor + Seq.Read(1).tRamp;
  if Seq.AQSlice(1).SpoilDephaseLengthFactor == 0
    % set minimum gradient length
    Spoil.GradDephaseLength = 2*Seq.Read(1).tRamp + 2/HW.RX(Seq.AQSlice(1).iDevice).fSample;
  end
  % center of spoiler in tRep with excitation pulse
  Spoil.CenterOfDephaseExcitation = Spoil.CenterOfDephase;
  Spoil.CenterOfRephase = Seq.Read(1).CenterOfDephase - ...
    (Seq.Read(1).GradDephaseLength-Seq.Read(1).tRamp) .* (Seq.AQSlice(1).SpoilRephaseLengthFactor-1) / 2 - ...
    Seq.AQSlice(1).SpoilRephaseTimeOffset;
  Spoil.GradRephaseLength = (Seq.Read(1).GradDephaseLength-Seq.Read(1).tRamp) .* Seq.AQSlice(1).SpoilRephaseLengthFactor + Seq.Read(1).tRamp;
  defGradDephaseLength = Seq.Read(1).GradDephaseLength;
  if Seq.AQSlice(1).RephaseTimeOffset ~= 0
    Seq.Read(1).CenterOfRephase = Seq.Read(1).CenterOfRephase + Seq.AQSlice(1).RephaseTimeOffset;
    Seq.Read(2).CenterOfRephase = Seq.Read(1).CenterOfRephase;
    Seq.Read(3).CenterOfRephase = Seq.Read(1).CenterOfRephase;
  end
  if Seq.AQSlice(1).DephaseTimeOffset ~= 0
    Seq.Read(1).CenterOfDephase = Seq.Read(1).CenterOfDephase - Seq.AQSlice(1).DephaseTimeOffset;
    Seq.Read(2).CenterOfDephase = Seq.Read(1).CenterOfDephase;
    Seq.Read(3).CenterOfDephase = Seq.Read(1).CenterOfDephase;
  end
  if Seq.AQSlice(1).RephaseTimeOffset~=0 || Seq.AQSlice(1).DephaseTimeOffset~=0
      % update readout timings
      Seq = get_ReadParameter(Seq, HW);
  end

  %% set up structures for get_SliceParameter (i.e. rf pulses with optional slice selection)
  % 90 degrees slice excitation
  Seq.Slice(1).Pulse.Function = Seq.AQSlice(1).excitationPulse;
  Seq.Slice(1).Pulse.FlipAngleComposite = Seq.AQSlice(1).excitationFlipAngleComposite;
  Seq.Slice(1).Pulse.FlipPhaseComposite = Seq.AQSlice(1).excitationFlipPhaseComposite;
  Seq.Slice(1).Pulse.FlipAngle = Seq.AQSlice(1).excitationFlipAngle + ...
    Seq.AQSlice(1).excitationFlipAngleIncrement*(Seq.Loop-1);
  % Seq.Slice(1).Pulse.MaxNumberOfSegments = 501;
  Seq.Slice(1).Thickness = Seq.AQSlice(1).thickness;
  Seq.Slice(1).CenterOfPulse = Seq.tEcho/2;
  Seq.Slice(1).UseCoordinate = Seq.AQSlice(1).SliceCoordinate;
  Seq.Slice(1).GradTimeDelay = Seq.AQSlice(1).SliceTimeDelay;
  Seq.Slice(1).UseAtRepetitionTime = Seq.P90tReps;
  Seq.Slice(1).Pulse.Phase = Seq.tRepTurboBlockPhaseExcitation(Seq.Slice(1).UseAtRepetitionTime);
  Seq.Slice(1).Pulse.PhaseIncrement = Seq.AQSlice(1).excitationPhaseIncrement;
  Seq.Slice(1).CenterOfRephase = [];
  Seq.Slice(1).GradDephaseLength = Spoil.GradRephaseLength;  % Seq.Read(1).GradDephaseLength;
  Seq.Slice(1).GradRephaseLength = ...
    (defGradDephaseLength-Seq.Read(1).tRamp) .* Seq.AQSlice(1).SliceRephaseLengthFactor + Seq.Read(1).tRamp;
  Seq.Slice(1).MaxGradAmp = Seq.MaxGradAmpSlice;
  Seq.Slice(1).distance = Seq.AQSlice(1).Center2OriginImage(1);
  Seq.Slice(1).GradTimeIntegralRephaseOffset = Seq.AQSlice(1).SliceGradTimeIntegralRephaseOffset;
  Seq.Slice(1).GradSign = Seq.AQSlice(1).SliceGradSign*Seq.AQSlice(1).GradSign;
  Seq.Slice(1).GradDephaseSign = -Seq.AQSlice(1).SliceGradSign * Seq.AQSlice(1).GradSign;
  Seq.Slice(1).GradRephaseSign = -Seq.AQSlice(1).SliceGradSign * Seq.AQSlice(1).GradSign;
  Seq.Slice(1).iDevice = Seq.AQSlice(1).iDevice;
  if isemptyfield(Seq.Slice(1), 'Gamma')
    if isemptyfield(Seq, 'AQSlice') ...
        || isemptyfield(Seq.AQSlice(1), 'Gamma')
      Seq.Slice(1).Gamma = HW.GammaDef;
    else
      Seq.Slice(1).Gamma = Seq.AQSlice(1).Gamma;
    end
  end
  if Seq.AQSlice(1).dualNuclearImage
    Seq.Slice(1).GammaX = Seq.AQSlice(1).GammaX;
    Seq.Slice(1).Pulse.FlipAngleX = Seq.Slice(1).Pulse.FlipAngle;
  elseif abs(Seq.Slice(1).Gamma - Seq.CorrectPhaseFrequencyGamma) > 5*eps(Seq.Slice(1).Gamma)
    % use dual frequency pulses to excite signal for both nuclei (tolerance for
    % comparison of floating point numbers was chosen arbitrarily)
    Seq.Slice(1).GammaX = Seq.CorrectPhaseFrequencyGamma;
    Seq.Slice(1).Pulse.FlipAngleX = Seq.CorrectPhaseFlipAngle;
  end

  Seq = get_SliceParameter(Seq, HW);  % generate slice parameters

  Seq.Slice(1).CenterOfRephase = Seq.Slice(1).CenterOfRephase + Seq.AQSlice(1).SliceRephaseTimeOffset;

  if Seq.CorrectPhase > 0
    % add acquisition window after acquisition pulses for frequency tracking
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
    Seq.AQSlice(2).nPhase(2) = numel(Seq.P90tReps);           % number of pixels in phase(2) direction
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
      % integer oversampling factor in read direction
      Seq.AQSlice(2).ReadOS = max(1, ...
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
    Seq.AQSlice(2).PhaseOS(1) = 1;                            % oversampling factor in phase(1) direction: 1...
    Seq.AQSlice(2).PhaseOS(2) = 1;                            % oversampling factor in phase(2) direction: 1...
    Seq.AQSlice(2).PhaseOS(3) = 1;                            % oversampling factor in phase(3) direction: 1...
    Seq.AQSlice(2).UsetRep = reshape(Seq.Slice(1).UseAtRepetitionTime, [], 1);  % directly after each excitation pulse
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

    % read out for frequency tracking
    Seq.Read(5).useAQSlice = 2;
    Seq.Read(5).HzPerPixelMin = Seq.AQSlice(2).HzPerPixMin;
    if Seq.AQSlice(1).thickness < 1000
      Seq.CorrectPhaseAQtOffset = ...
        max(Seq.Slice(1).CenterOfRephase + Seq.Slice(1).GradRephaseLength/2 ...
            + HW.Grad(Seq.AQSlice(1).iDevice).tEC - Seq.Slice(1).CenterOfPulse, ...
            Seq.CorrectPhaseAQtOffset);
    else
      Seq.CorrectPhaseAQtOffset = ...
        max(Seq.Slice(1).Pulse.CenterOffset + Seq.Slice(1).Pulse.MaxLength/2 ...
            + get_DeadTimeTX2RX(HW, Seq.AQSlice(2).HzPerPixMin*Seq.AQSlice(2).nRead*Seq.AQSlice(2).ReadOS, Seq.AQSlice(2).iDevice), ...
            Seq.CorrectPhaseAQtOffset);
    end
    Seq.Read(5).CenterOfReadout = Seq.Slice(1).CenterOfPulse + Seq.CorrectPhaseAQtOffset ...
      + 0.5/Seq.Read(5).HzPerPixelMin ...
      + (mod(Seq.AQSlice(2).nRead*Seq.AQSlice(2).ReadOS, 2)~=1)/Seq.Read(5).HzPerPixelMin/Seq.AQSlice(2).nRead/2;
    Seq.Read(5).PhaseIncrement = Seq.Read(1).PhaseIncrement;
    Seq.Read(5).Phase = Seq.tRepTurboBlockPhaseRead(1,:) + Seq.AQSlice(1).readOutPhase;
    % Seq.Read(5).Phase = Seq.Slice(1).Pulse.Phase(:) + 180;
    Seq.Read(5).UseAtRepetitionTime = Seq.AQSlice(2).UsetRep;
    if abs(Seq.Slice(1).Gamma - Seq.CorrectPhaseFrequencyGamma) > 5*eps(Seq.Slice(1).Gamma)
      Seq.Read(5).GammaX = Seq.CorrectPhaseFrequencyGamma;
    end

    % generate readout timings
    Seq = get_ReadParameter(Seq, HW);

    % adjust position of spoiler gradient slot in tRep with excitation if needed
    Spoil.CenterOfDephaseExcitation = Spoil.CenterOfDephaseExcitation + max(0, ...
      ...  % earliest start of gradient pulse to not collide with frequency tracking window
      (Seq.Read(5).CenterOfReadout + Seq.Read(5).AcquisitionTime/2 + get_DeadTimeRX2TX(HW, Seq.Read(5).fSample, Seq.AQSlice(2).iDevice)) ...
      ...   % current start of gradient pulse:
      - (Spoil.CenterOfDephase - Spoil.GradDephaseLength/2));
  end

  readSpoilCoordinate = (Seq.AQSlice(1).PhaseCoordinate == Seq.AQSlice(1).ReadCoordinate);
  if Seq.AQSlice(1).sizePhaseSpoil(readSpoilCoordinate) == 0  % set spoiler amp equal to read grad amp
    % Spoil read direction
    if Seq.Read(1).GradTimeIntegral == 0  % no readout gradient
      Seq.AQSlice(1).SpoilFactor(readSpoilCoordinate) = 0;
    else
      Seq.AQSlice(1).SpoilFactor(readSpoilCoordinate) = ...
        (~Seq.DephaseBefore180) * Seq.Read(1).GradTimeIntegral/Seq.Read(1).GradTimeIntegralAQ/2 + ...  % compensate read dephase
        (Seq.Read(1).GradDephaseLength-Seq.Read(1).tRamp)/Seq.Read(1).AcquisitionTime;  % spoil such that effective amplitude equals the read out gradient amplitude
    end

    Seq.AQSlice(1).sizePhaseSpoil(readSpoilCoordinate) = ...
      Seq.AQSlice(1).sizeRead ./ Seq.AQSlice(1).nRead ./ Seq.AQSlice(1).SpoilFactor(readSpoilCoordinate)./2;
  end
  if Seq.AQSlice(1).sizePhaseSpoilEnd(readSpoilCoordinate) == 0
    % Set spoiler amp equal to read grad amp after last read-out window in echo
    % train.
    % Spoil read direction
    if Seq.Read(1).GradTimeIntegral == 0  % no readout gradient
      Seq.AQSlice(1).SpoilFactorEnd(readSpoilCoordinate) = 0;
    else
      Seq.AQSlice(1).SpoilFactorEnd(readSpoilCoordinate) = ...
        (~Seq.DephaseBefore180) * Seq.Read(1).GradTimeIntegral/Seq.Read(1).GradTimeIntegralAQ/2 + ...  % compensate read dephase
        (Seq.Read(1).GradDephaseLength-Seq.Read(1).tRamp)/Seq.Read(1).AcquisitionTime;  % spoil such that effective amplitude equals the read out gradient amplitude
    end

    Seq.AQSlice(1).sizePhaseSpoilEnd(readSpoilCoordinate) = ...
      Seq.AQSlice(1).sizeRead ./ Seq.AQSlice(1).nRead ./ Seq.AQSlice(1).SpoilFactorEnd(readSpoilCoordinate)./2;
  end

  % 180 degrees (slice) inversion
  Seq.Slice(2).Pulse.Function = Seq.AQSlice(1).inversionPulse;
  Seq.Slice(2).Pulse.FlipAngleComposite = Seq.AQSlice(1).inversionFlipAngleComposite;
  Seq.Slice(2).Pulse.FlipPhaseComposite = Seq.AQSlice(1).inversionFlipPhaseComposite;
  Seq.Slice(2).Pulse.FlipAngle = Seq.AQSlice(1).inversionFlipAngle + ...
    Seq.AQSlice(1).inversionFlipAngleIncrement * (Seq.Loop-1);
  % Seq.Slice(2).Pulse.MaxNumberOfSegments= 401;
  Seq.Slice(2).Thickness = Seq.AQSlice(1).thicknessInversion;  % use slice selection on inversion pulses
  Seq.Slice(2).CenterOfPulse = 0;
  Seq.Slice(2).GradDephaseSign = -1;
  Seq.Slice(2).UseCoordinate = Seq.AQSlice(1).SliceCoordinateInvert;
  Seq.Slice(2).UseAtRepetitionTime = sort([Seq.PreShotstReps; Seq.P180tReps; Seq.PostShotstReps]);
  Seq.Slice(2).UseAtRepetitionTimeDephase = Seq.Slice(2).UseAtRepetitionTime - 1;
  Seq.Slice(2).Pulse.Phase = Seq.tRepTurboBlockPhaseInversion(Seq.Slice(2).UseAtRepetitionTime) + ...
    Seq.AQSlice(1).inversionPhase;
  Seq.Slice(2).Pulse.PhaseIncrement = 0;
  Seq.Slice(2).CenterOfDephase = Seq.Read(1).CenterOfRephase;  % Spoil.CenterOfRephase-mean(Seq.Read(1).GradTimeDelayOffset(1));
  Seq.Slice(2).GradDephaseLength = Seq.Read(1).GradRephaseLength;  % Spoil.GradRephaseLength;
  Seq.Slice(2).CenterOfRephase = Seq.Read(1).CenterOfDephase;  % Spoil.CenterOfDephase-mean(Seq.Read(1).GradTimeDelayOffset(1));
  Seq.Slice(2).GradRephaseLength = Seq.Read(1).GradDephaseLength;  % Spoil.GradDephaseLength;
  Seq.Slice(2).MaxGradAmp = Seq.MaxGradAmpInversion;
  Seq.Slice(2).distance = Seq.AQSlice(1).Center2OriginImage(1);
  Seq.Slice(2).GradSign = Seq.AQSlice(1).SliceGradSign * Seq.AQSlice(1).GradSign;
  Seq.Slice(2).GradDephaseSign = -Seq.AQSlice(1).SliceGradSign * Seq.AQSlice(1).GradSign;
  Seq.Slice(2).GradRephaseSign = -Seq.AQSlice(1).SliceGradSign * Seq.AQSlice(1).GradSign;
  Seq.Slice(2).iDevice = Seq.AQSlice(1).iDevice;
  if Seq.AQSlice(1).dualNuclearImage
    Seq.Slice(2).GammaX = Seq.AQSlice(1).GammaX;
    Seq.Slice(2).Pulse.FlipAngleX = Seq.Slice(2).Pulse.FlipAngle;
  end

  % Move slice dephase of inversion pulse to slice rephase of excitation pulse
  % for first inversion.
  Seq.Slice(3) = Seq.Slice(2);
  Seq.Slice(3).Pulse.Phase = [];
  Seq.Slice(3).UseAtRepetitionTime = Seq.Slice(1).UseAtRepetitionTime;
  Seq.Slice(3).UseAtRepetitionTimeRephase = Seq.Slice(1).UseAtRepetitionTime;
  Seq.Slice(3).UseAtRepetitionTimeDephase = Seq.Slice(1).UseAtRepetitionTime;
  Seq.Slice(3).CenterOfRephase = Seq.Slice(2).CenterOfDephase;
  Seq.Slice(3).CenterOfDephase = Seq.Slice(1).CenterOfRephase;
  Seq.Slice(3).GradRephaseLength = Seq.Slice(2).GradDephaseLength;
  Seq.Slice(3).GradDephaseLength = Seq.Slice(1).GradRephaseLength;
  Seq.Slice(3).GradRephaseSign = -Seq.Slice(2).GradDephaseSign;
  if Seq.AQSlice(1).dualNuclearImage
    Seq.Slice(3).GammaX = Seq.AQSlice(1).GammaX;
    Seq.Slice(3).Pulse.FlipAngleX = Seq.Slice(3).Pulse.FlipAngle;
  end

  if Seq.AQSlice(1).SpoilFactorInversionBlock
    Seq.Slice(2).tEC = 0;
    Seq = get_SliceParameter(Seq, HW);  % generate slice parameters with tEC = 0
    Seq.Slice(2).tEC = max(HW.Grad(Seq.Slice(2).iDevice).tEC, ...
      ((pi/HW.GammaDef/(Seq.AQSlice(1).sizeSpoilInversionBlock/2)) - Seq.Slice(2).GradTimeIntegral/2) / ...
      abs3D(Seq.Slice(2).GradAmp));
    if HW.Grad(Seq.Slice(2).iDevice).tEC > ((pi/HW.GammaDef/(Seq.AQSlice(1).sizeSpoilInversionBlock/2))-Seq.Slice(2).GradTimeIntegral/2)/abs3D(Seq.Slice(2).GradAmp)
      warning('PD:sequence_Spin_Echo:SpoilFactorInversionBlock', ...
        ['Seq.AQSlice(1).SpoilFactorInversionBlock spoiler is too weak, it will be increased by the factor of ', ...
        num2str((Seq.Slice(2).GradTimeIntegral/2 + HW.Grad(Seq.Slice(2).iDevice).tEC*abs3D(Seq.Slice(2).GradAmp)) / ...
        (pi/HW.GammaDef/(Seq.AQSlice(1).sizeSpoilInversionBlock/2)), 3)]);
    end
    % generate slice new parameters
    Seq = get_SliceParameter(Seq, HW);
    Seq.AQSlice(1).inversionGradDephase = 0;
    Seq.AQSlice(1).inversionGradRephase = 0;
  else
    Seq = get_SliceParameter(Seq, HW);  % generate slice parameters
  end

  if ~isfield(Seq.AQSlice(1), 'inversionGradDephase') ...
      || isempty(Seq.AQSlice(1).inversionGradDephase) ...
      || isnan(Seq.AQSlice(1).inversionGradDephase)
    Seq.AQSlice(1).inversionGradDephase = ~Seq.DephaseBefore180;
  end
  % FIXME: Does it make sense to add the dephase pulse but not the rephase pulse
  % (or vice versa)? They should probably either both be present or both be
  % omitted.
  if ~isfield(Seq.AQSlice(1), 'inversionGradRephase') ...
      || isempty(Seq.AQSlice(1).inversionGradRephase) ...
      || isnan(Seq.AQSlice(1).inversionGradRephase)
    Seq.AQSlice(1).inversionGradRephase = Seq.AQSlice(1).inversionGradDephase;
  end


  % move read dephase after excitation pulse for DephaseBefore180
  Seq.Read(4) = Seq.Read(1);
  Seq.Read(4).Phase = [];
  Seq.Read(4).UseAtRepetitionTime = Seq.Read(1).UseAtRepetitionTime(1,:);
  Seq.Read(4).UseAtRepetitionTimeDephase = Seq.Read(1).UseAtRepetitionTimeDephase(1,:);
  Seq.Read(4).UseAtRepetitionTimeRephase = Seq.Read(1).UseAtRepetitionTimeDephase(1,:);
  Seq.Read(4).CenterOfRephase = Seq.Read(1).CenterOfDephase;
  Seq.Read(4).GradRephaseLength = Seq.Read(1).GradDephaseLength;
  Seq.Read(4).GradRephaseSign = -Seq.Read(1).GradDephaseSign;
  if Seq.CorrectPhase > 0
    % read dephase pulse at phase dephase slot
    Seq.Read(4).CenterOfDephase = Spoil.CenterOfDephaseExcitation;
    Seq.Read(4).GradDephaseLength = Spoil.GradDephaseLength;
  else
    % read dephase pulse at slice rephase slot
    Seq.Read(4).CenterOfDephase = Seq.Slice(1).CenterOfRephase;
    Seq.Read(4).GradDephaseLength = Seq.Slice(1).GradRephaseLength;
  end


  %% set up structures for get_PhaseParameter (i.e. phase encoding pulses and spoilers)
  % phase(1) encoding aligned with read dephase and rephase
  % By default, phase(1) is parallel to the slice direction
  Seq.Phase(1).sizePhase = Seq.AQSlice(1).sizePhase(1);
  Seq.Phase(1).nPhase = Seq.AQSlice(1).nPhase(1);
  Seq.Phase(1).PhaseOS = Seq.AQSlice(1).PhaseOS(1);
  if Seq.DephaseBefore180
    Seq.Phase(1).CenterOfDephase = Seq.Read(4).CenterOfDephase;
    Seq.Phase(1).GradDephaseLength = Seq.Read(4).GradDephaseLength;
  else
    Seq.Phase(1).CenterOfDephase = Seq.Read(1).CenterOfDephase;
    Seq.Phase(1).GradDephaseLength = Seq.Read(1).GradDephaseLength;
  end
  Seq.Phase(1).CenterOfRephase = Seq.Read(1).CenterOfRephase;
  Seq.Phase(1).GradRephaseLength = Seq.Read(1).GradRephaseLength;
  Seq.Phase(1).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(1);
  Seq.Phase(1).GradTimeDelay = Seq.AQSlice(1).ReadTimeDelay + Seq.Read(1).GradTimeDelayOffset(1);
  if Seq.DephaseBefore180
    Seq.Phase(1).StepIncrement = (Seq.AQSlice(1).oddEvenEchoes+1);
  else
    Seq.Phase(1).StepIncrement = (Seq.AQSlice(1).oddEvenEchoes+1) * Seq.AQSlice(1).nImages;
  end
  Seq.Phase(1).UseAtRepetitionTime = Seq.AQSlice(1).UsetRep;
  Seq.Phase(1).UseAtRepetitionTimeDephase = Seq.Phase(1).UseAtRepetitionTime - ...
    Seq.DephaseBefore180*(1+Seq.SteadyState_PreShots180);
  Seq.Phase(1).usedkLines = sort(Seq.AQSlice(1).kLineOrder);
  Seq.Phase(1).numKLines = numKLines;
  Seq.Phase(1).kLineIncrement = Seq.Phase(1).StepIncrement;
  Seq.Phase(1).distance = Seq.AQSlice(1).Center2OriginImage(1);
  % Toggle sign of dephase pulse (for Seq.DephaseBefore180) to have the "normal"
  % ordering of phase lines at the first echo with AQ window.
  Seq.Phase(1).GradDephaseSign = 1*Seq.AQSlice(1).GradSign * ...
    (-1)^(Seq.DephaseBefore180*(1+Seq.SteadyState_PreShots180));
  % Toggle sign of rephase pulse (for Seq.DephaseBefore180) to refocus whichever
  % ordering we have at the last echo of the echo train.
  Seq.Phase(1).GradRephaseSign = -1*Seq.AQSlice(1).GradSign * ...
    (-1)^(Seq.DephaseBefore180*(1+size(Seq.P180tReps,1)));

  % phase(2) encoding aligned with read dephase and rephase
  % By default, phase(2) is perpendicular to the read and slice directions
  Seq.Phase(2) = Seq.Phase(1);
  Seq.Phase(2).sizePhase = Seq.AQSlice(1).sizePhase(2);
  Seq.Phase(2).nPhase = Seq.AQSlice(1).nPhase(2);
  Seq.Phase(2).PhaseOS = Seq.AQSlice(1).PhaseOS(2);
  Seq.Phase(2).StepIncrement = Seq.AQSlice(1).nPhase(1) * Seq.AQSlice(1).PhaseOS(1) * ... % phase(1) is parallel to the slice direction.
    Seq.Phase(1).StepIncrement;
  Seq.Phase(2).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(2);
  Seq.Phase(2).distance = Seq.AQSlice(1).Center2OriginImage(2);

  % phase(3) encoding aligned with read dephase and rephase
  % By default, phase(3) is parallel to the read direction
  Seq.Phase(3) = Seq.Phase(1);
  Seq.Phase(3).sizePhase = Seq.AQSlice(1).sizePhase(3);
  Seq.Phase(3).nPhase = Seq.AQSlice(1).nPhase(3);
  Seq.Phase(3).PhaseOS = Seq.AQSlice(1).PhaseOS(3);
  Seq.Phase(3).StepIncrement = Seq.AQSlice(1).nPhase(2) * Seq.AQSlice(1).PhaseOS(2) * ...
    Seq.Phase(2).StepIncrement;
  Seq.Phase(3).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(3);
  Seq.Phase(3).distance = Seq.AQSlice(1).Center2OriginImage(3);

  % spoiler at refocusing pulses aligned with read dephase and rephase
  % This spoiler is in phase(1)/slice direction.
  Seq.Phase(4).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(1);
  Seq.Phase(4).nPhase = 1;
  Seq.Phase(4).PhaseOS = 2;
  Seq.Phase(4).StepOrder = [1,1];
  % dephase signal before inversion
  Seq.Phase(4).CenterOfDephase = Spoil.CenterOfDephase;
  Seq.Phase(4).GradDephaseLength = Spoil.GradDephaseLength;
  Seq.Phase(4).GradDephaseSign = -Seq.AQSlice(1).SliceGradSign * Seq.AQSlice(1).GradSign;
  % rephase signal after inversion
  Seq.Phase(4).CenterOfRephase = Spoil.CenterOfRephase;
  Seq.Phase(4).GradRephaseLength = Spoil.GradRephaseLength;
  Seq.Phase(4).GradRephaseSign = -Seq.AQSlice(1).SliceGradSign * Seq.AQSlice(1).GradSign;
  Seq.Phase(4).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(1);
  Seq.Phase(4).GradTimeDelay = Seq.AQSlice(1).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(1);
  if Seq.SteadyState_PreShots180SpoilerGrad
    Seq.Phase(4).UseAtRepetitionTime = Seq.Slice(2).UseAtRepetitionTime;
  else
    Seq.Phase(4).UseAtRepetitionTime = Seq.P180tReps;
  end
  Seq.Phase(4).UseAtRepetitionTimeDephase = Seq.Phase(4).UseAtRepetitionTime - 1;

  % spoiler at refocusing pulses aligned with read dephase and rephase
  % This spoiler is in phase(2) direction.
  Seq.Phase(5) = Seq.Phase(4);
  Seq.Phase(5).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(2);
  Seq.Phase(5).nPhase = 1;
  Seq.Phase(5).PhaseOS = 2;
  Seq.Phase(5).StepOrder = [1, 1];
  Seq.Phase(5).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(2);
  Seq.Phase(5).GradTimeDelay = Seq.AQSlice(1).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(1);
  Seq.Phase(5).GradDephaseSign = -Seq.AQSlice(1).PhaseGradSign * Seq.AQSlice(1).GradSign;
  Seq.Phase(5).GradRephaseSign = -Seq.AQSlice(1).PhaseGradSign * Seq.AQSlice(1).GradSign;
  Seq.Phase(5).GradTimeIntegralRephaseOffset = -Seq.AQSlice(1).PhaseGradTimeIntegralRephaseOffset;

  % spoiler at refocusing pulses aligned with read dephase and rephase
  % This spoiler is in phase(3)/read direction.
  Seq.Phase(6) = Seq.Phase(4);
  Seq.Phase(6).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(3);
  Seq.Phase(6).nPhase = 1;
  Seq.Phase(6).PhaseOS = 2;
  Seq.Phase(6).StepOrder = [1, 1];
  Seq.Phase(6).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(3);
  Seq.Phase(6).GradTimeDelay = Seq.AQSlice(1).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(1);
  Seq.Phase(6).GradDephaseSign = -Seq.AQSlice(1).ReadGradSign * Seq.AQSlice(1).GradSign;
  Seq.Phase(6).GradRephaseSign = -Seq.AQSlice(1).ReadGradSign * Seq.AQSlice(1).GradSign;

  % spoiler at end of echo train aligned with read rephase (but might be of differing length)
  % This spoiler is in phase(1)/slice direction.
  Seq.Phase(7).sizePhase = Seq.AQSlice(1).sizePhaseSpoilEnd(1);
  Seq.Phase(7).nPhase = 1;
  Seq.Phase(7).PhaseOS = 2;
  Seq.Phase(7).StepOrder = [1, 1];
  Seq.Phase(7).CenterOfRephase = Spoil.CenterOfRephase; % not used
  Seq.Phase(7).CenterOfDephase = Spoil.CenterOfDephaseEnd;
  Seq.Phase(7).GradRephaseLength = Spoil.GradRephaseLength;% not used
  Seq.Phase(7).GradDephaseLength = Spoil.GradDephaseLengthEnd;
  Seq.Phase(7).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(1);
  Seq.Phase(7).GradTimeDelay = Seq.AQSlice(1).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(1);
  Seq.Phase(7).UseAtRepetitionTime = Seq.P90tReps + ...
    (Seq.AQSlice(1).TurboFactor * Seq.AQSlice(1).nImages * (Seq.AQSlice(1).oddEvenEchoes+1) + ...
    Seq.SteadyState_PreShots180 + Seq.SteadyState_PostShots180);
  Seq.Phase(7).UseAtRepetitionTimeRephase = Seq.Phase(7).UseAtRepetitionTime;
  Seq.Phase(7).GradDephaseSign = -Seq.AQSlice(1).SliceGradSign * Seq.AQSlice(1).GradSign;
  Seq.Phase(7).GradRephaseSign = -Seq.AQSlice(1).SliceGradSign * Seq.AQSlice(1).GradSign;

  % spoiler at end of echo train aligned with read rephase
  % This spoiler is in phase(2) direction.
  Seq.Phase(8) = Seq.Phase(7);
  Seq.Phase(8).sizePhase = Seq.AQSlice(1).sizePhaseSpoilEnd(2);
  Seq.Phase(8).nPhase = 1;
  Seq.Phase(8).PhaseOS = 2;
  Seq.Phase(8).StepOrder = [1, 1];
  Seq.Phase(8).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(2);
  Seq.Phase(8).GradTimeDelay = Seq.AQSlice(1).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(1);
  Seq.Phase(8).GradDephaseSign = -Seq.AQSlice(1).PhaseGradSign * Seq.AQSlice(1).GradSign;
  Seq.Phase(8).GradRephaseSign = -Seq.AQSlice(1).PhaseGradSign * Seq.AQSlice(1).GradSign;
  Seq.Phase(8).GradTimeIntegralRephaseOffset = -Seq.AQSlice(1).PhaseGradTimeIntegralRephaseOffset;

  % spoiler at end of echo train aligned with read rephase
  % This spoiler is in phase(3)/read direction.
  Seq.Phase(9) = Seq.Phase(7);
  Seq.Phase(9).sizePhase = Seq.AQSlice(1).sizePhaseSpoilEnd(3);
  Seq.Phase(9).nPhase = 1;
  Seq.Phase(9).PhaseOS = 2;
  Seq.Phase(9).StepOrder = [1, 1];
  Seq.Phase(9).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(3);
  Seq.Phase(9).GradTimeDelay = Seq.AQSlice(1).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(1);
  Seq.Phase(9).GradDephaseSign = -Seq.AQSlice(1).ReadGradSign * Seq.AQSlice(1).GradSign;
  Seq.Phase(9).GradRephaseSign = -Seq.AQSlice(1).ReadGradSign * Seq.AQSlice(1).GradSign;

  % move spoiler from before first inversion pulse with spoiler
  % to before first inversion pulse (aligned with excitation slice rephase)
  % This spoiler is in phase(1)/slice direction.
  Seq.Phase(10).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(1);
  Seq.Phase(10).nPhase = 1;
  Seq.Phase(10).PhaseOS = 2;
  Seq.Phase(10).StepOrder = [1, 1];
  if Seq.CorrectPhase < 1
    Seq.Phase(10).CenterOfDephase = ...
      Seq.Slice(1).CenterOfRephase + Seq.Read(1).GradTimeDelayOffset(1);
    Seq.Phase(10).GradDephaseLength = Seq.Slice(1).GradRephaseLength;
  else
    Seq.Phase(10).CenterOfDephase = Spoil.CenterOfDephaseExcitation;
    Seq.Phase(10).GradDephaseLength = Seq.Phase(4).GradDephaseLength;
  end
  Seq.Phase(10).CenterOfRephase = Seq.Phase(4).CenterOfDephase;
  Seq.Phase(10).GradRephaseLength = Seq.Phase(4).GradDephaseLength;
  Seq.Phase(10).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(1);
  Seq.Phase(10).GradTimeDelay = Seq.AQSlice(1).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(1);
  Seq.Phase(10).UseAtRepetitionTime = Seq.P90tReps;
  Seq.Phase(10).UseAtRepetitionTimeDephase = Seq.P90tReps;
  Seq.Phase(10).UseAtRepetitionTimeRephase = Seq.Phase(4).UseAtRepetitionTimeDephase(1,:);
  Seq.Phase(10).GradDephaseSign = -Seq.AQSlice(1).SliceGradSign * Seq.AQSlice(1).GradSign .* ...
    sign(mod(diff([Seq.P90tReps(1,1); Seq.Phase(4).UseAtRepetitionTime(1,1)], 1, 1), 2) - 0.5);
  Seq.Phase(10).GradRephaseSign = Seq.AQSlice(1).SliceGradSign * Seq.AQSlice(1).GradSign;

  % move spoiler from before first inversion pulse with spoiler
  % to before first inversion pulse (aligned with excitation slice rephase)
  % This spoiler is in phase(2) direction.
  Seq.Phase(11) = Seq.Phase(10);
  Seq.Phase(11).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(2);
  Seq.Phase(11).nPhase = 1;
  Seq.Phase(11).PhaseOS = 2;
  Seq.Phase(11).StepOrder = [1, 1];
  Seq.Phase(11).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(2);
  Seq.Phase(11).GradTimeDelay = Seq.AQSlice(1).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(1);
  Seq.Phase(11).UseAtRepetitionTime = Seq.P90tReps;
  Seq.Phase(11).GradDephaseSign = -Seq.AQSlice(1).PhaseGradSign * Seq.AQSlice(1).GradSign .* ...
    sign(mod(diff([Seq.P90tReps(1,1); Seq.Phase(5).UseAtRepetitionTime(1,1)], 1, 1), 2) - 0.5);
  Seq.Phase(11).GradRephaseSign = Seq.AQSlice(1).PhaseGradSign * Seq.AQSlice(1).GradSign;

  % move spoiler from before first inversion pulse with spoiler
  % to before first inversion pulse (aligned with excitation slice rephase)
  % This spoiler is in phase(3)/read direction.
  Seq.Phase(12) = Seq.Phase(10);
  Seq.Phase(12).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(3);
  Seq.Phase(12).nPhase = 1;
  Seq.Phase(12).PhaseOS = 2;
  Seq.Phase(12).StepOrder = [1, 1];
  Seq.Phase(12).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(3);
  Seq.Phase(12).GradTimeDelay = Seq.AQSlice(1).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(1);
  Seq.Phase(12).UseAtRepetitionTime = Seq.P90tReps;
  Seq.Phase(12).GradDephaseSign = -Seq.AQSlice(1).ReadGradSign * Seq.AQSlice(1).GradSign .* ...
    sign(mod(diff([Seq.P90tReps(1,1); Seq.Phase(6).UseAtRepetitionTime(1,1)], 1, 1), 2) - 0.5);
  Seq.Phase(12).GradRephaseSign = Seq.AQSlice(1).ReadGradSign * Seq.AQSlice(1).GradSign;

  % move spoiler inside tRep with excitation (if Seq.CorrectPhase is true)
  % This spoiler is in phase(1)/slice direction.
  Seq.Phase(13).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(1);
  Seq.Phase(13).nPhase = 1;
  Seq.Phase(13).PhaseOS = 2;
  Seq.Phase(13).StepOrder = [1, 1];
  Seq.Phase(13).CenterOfDephase = Seq.Phase(4).CenterOfDephase;
  Seq.Phase(13).CenterOfRephase = Spoil.CenterOfDephaseExcitation;
  Seq.Phase(13).GradDephaseLength = Seq.Phase(4).GradDephaseLength;
  Seq.Phase(13).GradRephaseLength = Seq.Phase(4).GradDephaseLength;
  Seq.Phase(13).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(1);
  Seq.Phase(13).GradTimeDelay = Seq.AQSlice(1).PhaseTimeDelay + Seq.Read(1).GradTimeDelayOffset(1);
  Seq.Phase(13).UseAtRepetitionTime = Seq.P90tReps;

  % move spoiler inside tRep with excitation (if Seq.CorrectPhase is true)
  % This spoiler is in phase(2) direction.
  Seq.Phase(14) = Seq.Phase(13);
  Seq.Phase(14).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(2);
  Seq.Phase(14).CenterOfDephase = Seq.Phase(5).CenterOfDephase;
  Seq.Phase(14).CenterOfRephase = Spoil.CenterOfDephaseExcitation;
  Seq.Phase(14).GradDephaseLength = Seq.Phase(5).GradDephaseLength;
  Seq.Phase(14).GradRephaseLength = Seq.Phase(5).GradDephaseLength;
  Seq.Phase(14).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(2);
  Seq.Phase(14).GradTimeDelay = Seq.AQSlice(1).PhaseTimeDelay + Seq.Read(2).GradTimeDelayOffset(1);
  Seq.Phase(14).UseAtRepetitionTime = Seq.P90tReps;

  % move spoiler inside tRep with excitation (if Seq.CorrectPhase is true)
  % This spoiler is in phase(3)/read direction.
  Seq.Phase(15) = Seq.Phase(13);
  Seq.Phase(15).sizePhase = Seq.AQSlice(1).sizePhaseSpoil(3);
  Seq.Phase(15).CenterOfDephase = Seq.Phase(6).CenterOfDephase;
  Seq.Phase(15).CenterOfRephase = Spoil.CenterOfDephaseExcitation;
  Seq.Phase(15).GradDephaseLength = Seq.Phase(6).GradDephaseLength;
  Seq.Phase(15).GradRephaseLength = Seq.Phase(6).GradDephaseLength;
  Seq.Phase(15).UseCoordinate = Seq.AQSlice(1).PhaseCoordinate(3);
  Seq.Phase(15).GradTimeDelay = Seq.AQSlice(1).PhaseTimeDelay + Seq.Read(3).GradTimeDelayOffset(1);
  Seq.Phase(15).UseAtRepetitionTime = Seq.P90tReps;


  % apply all phase gradients to the same device
  [Seq.Phase(:).iDevice] = deal(Seq.AQSlice(1).iDevice);

  % generate phase parameters
  Seq = get_PhaseParameter(Seq, HW);

  % Image shift in phase direction
  [~, b] = sort(Seq.kLines(:));
  Seq.Read(1).Phase(b) = reshape(Seq.Read(1).Phase(b), 1, []) + ...
    reshape(Seq.Phase(1).AQPhaseShift + ...
            Seq.Phase(2).AQPhaseShift + ...
            Seq.Phase(3).AQPhaseShift, 1, [])/pi*180;
  % update readout timings
  Seq = get_ReadParameter(Seq, HW);


  if Seq.CorrectPhase > 0
    % check for collision of frequency tracking window
    endOfFreqTrack = Seq.Read(5).CenterOfReadout + Seq.Read(5).AcquisitionTime/2 ...
      + get_DeadTimeRX2TX(HW, Seq.Read(5).fSample, Seq.AQSlice(2).iDevice);
    if Seq.DephaseBefore180 && ...
        endOfFreqTrack > Seq.Phase(1).CenterOfDephase - Seq.Phase(1).GradDephaseLength/2 + Seq.tRep(Seq.Phase(1).UseAtRepetitionTime(1)-1)
      % with dephase gradient pulse
      error('PD:sequence_Spin_Echo:tEchoTooShortForFreqTrack', ...
        ['tEcho too short for selected frequency tracking AQ length ', ...
        '(%.3f ms too long)'], ...
        (endOfFreqTrack - (Seq.Phase(1).CenterOfDephase - Seq.Phase(1).GradDephaseLength/2))*1e3);
    end
    if Seq.Slice(2).Thickness < 1000 && ...
        endOfFreqTrack > Seq.Slice(2).CenterOfDephase - Seq.Slice(2).GradDephaseLength/2
      % with inversion pulse slice dephase gradient
      error('PD:sequence_Spin_Echo:tEchoTooShortForFreqTrack', ...
        ['tEcho too short for selected frequency tracking AQ length ', ...
        '(%.3f ms too long)'], ...
        (endOfFreqTrack - (Seq.Slice(2).CenterOfDephase - Seq.Slice(2).GradDephaseLength/2))*1e3);
    end
    if endOfFreqTrack > ...
        Seq.tRep(Seq.Slice(1).UseAtRepetitionTime(1)) + Seq.Slice(2).CenterOfPulse - Seq.Slice(2).Pulse.MaxLength/2 - Seq.Slice(2).Pulse.CenterOffset
      % with inversion pulse
      error('PD:sequence_Spin_Echo:tEchoTooShortForFreqTrack', ...
        ['tEcho too short for selected frequency tracking AQ length ', ...
        '(%.3f ms too long)'], ...
        (endOfFreqTrack - (Seq.tRep(Seq.Slice(1).UseAtRepetitionTime(1)) + Seq.Slice(2).CenterOfPulse - Seq.Slice(2).Pulse.MaxLength/2 - Seq.Slice(2).Pulse.CenterOffset))*1e3);
    end
    if any(~isinf(Seq.AQSlice(1).sizePhaseSpoil)) && ...
        endOfFreqTrack > ...
        Seq.Phase(13).CenterOfRephase - Seq.Phase(13).GradRephaseLength/2
      % with spoiler at inversion pulse
      error('PD:sequence_Spin_Echo:tEchoTooShortForFreqTrack', ...
        ['tEcho too short for selected frequency tracking AQ length ', ...
        '(%.3f ms too long)'], ...
        (endOfFreqTrack - (Seq.Phase(13).CenterOfDephase - Seq.Phase(13).GradDephaseLength/2))*1e3);
    end

    Seq.AQSlice(2).AcquisitionTime = Seq.Read(5).AcquisitionTime;
    Seq.AQSlice(2).AcquisitionFrequency = Seq.Read(5).AcquisitionFrequency;
    Seq.AQSlice(2).UsetRep = Seq.Read(5).UseAtRepetitionTime;
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

  %% check timing
  if ~Seq.DephaseBefore180
    % timing between inversion pulse and phase dephase slot
    if Seq.AQSlice(1).thicknessInversion < 1
      if Seq.Slice(2).CenterOfPulse + Seq.Slice(2).GradLength/2 ...
          > Seq.Read(1).CenterOfDephase - Seq.Read(1).GradDephaseLength/2
        error('PD:sequence_Spin_Echo:tEchoTooShortForAQ', ...
          'tEcho is %.3f ms too short for selected AQ length (1/HzPerPixelMin)', ...
          ((Seq.Slice(2).CenterOfPulse + Seq.Slice(2).GradLength/2) ...
           - (Seq.Read(1).CenterOfDephase - Seq.Read(1).GradDephaseLength/2))*2*1e3);
      end
    else
      if Seq.Slice(2).CenterOfPulse + Seq.Slice(2).Pulse.MaxLength/2 + ...
          Seq.Slice(2).Pulse.CenterOffset + Seq.Slice(2).tEC ...
          > Seq.Read(1).CenterOfDephase - Seq.Read(1).GradDephaseLength/2
        if any([Seq.Read(1).GradTimeIntegral , [Seq.Phase(1:12).GradTimeIntegral ] ])
          error('PD:sequence_Spin_Echo:tEchoTooShortForAQ', ...
            'tEcho is %.3f ms too short for selected AQ length (1/HzPerPixelMin)', ...
            ((Seq.Slice(2).CenterOfPulse + Seq.Slice(2).Pulse.MaxLength/2 + ...
              Seq.Slice(2).Pulse.CenterOffset + Seq.Slice(2).tEC) ...
             - (Seq.Read(1).CenterOfDephase - Seq.Read(1).GradDephaseLength/2))*2*1e3);
        end
      end
    end

    % timing between excitation pulse and slice rephase slot
    if Seq.AQSlice(1).thickness<1
      if Seq.Slice(1).CenterOfPulse + Seq.Slice(1).GradLength/2 - ...
          Seq.Slice(1).tRamp - 1e-9 ...
          > Seq.Slice(1).CenterOfRephase - Seq.Slice(1).GradRephaseLength/2
        % error('increase thickness or MaxGradAmpSlice, decrease HzPerPixMin, use shorter excitationPulse')
      end
    else
      if Seq.Slice(1).CenterOfPulse + Seq.Slice(1).Pulse.MaxLength/2 + ...
          Seq.Slice(1).Pulse.CenterOffset + Seq.Slice(1).tEC ...
          > Seq.Slice(1).CenterOfRephase - Seq.Slice(1).GradRephaseLength/2
         % error('decrease HzPerPixMin or use shorter excitationPulse')
      end
    end
  end

  % timing between slice rephase slot and inversion pulse slot
  sliceRephaseEnd = -Inf;
  if Seq.DephaseBefore180
    sliceRephaseEnd = max(sliceRephaseEnd, Seq.Phase(1).CenterOfDephase + Seq.Phase(1).GradDephaseLength/2);
  end
  if Seq.AQSlice(1).thickness<1
    sliceRephaseEnd = max(sliceRephaseEnd, Seq.Slice(1).CenterOfRephase + Seq.Slice(1).GradRephaseLength/2);
  end
  if Seq.AQSlice(1).thicknessInversion < 1
    sliceRephaseEnd = max(sliceRephaseEnd, Seq.Slice(2).CenterOfDephase + Seq.Slice(2).GradDephaseLength/2);
  end
  if any(isfinite(Seq.AQSlice(1).sizePhaseSpoil))
    sliceRephaseEnd = max(sliceRephaseEnd, max([Seq.Phase(10:12).CenterOfDephase] + [Seq.Phase(10:12).GradDephaseLength]/2));
  end
  % The different constituents of this gradient slot might cancel each other
  % out. In that case, a collision of the slots is not an issue. It is difficult
  % to generally check for that.
  % As the next best thing, emit a warning and rely on the user acting upon it
  % if necessary.
  if Seq.AQSlice(1).thicknessInversion < 1
    if sliceRephaseEnd > Seq.Slice(2).CenterOfPulse - Seq.Slice(2).GradLength/2 + Seq.tEcho
      warning('PD:sequence_Spin_Echo:SliceRephaseTooLate', ...
        ['The slice rephase slot might be colliding with the inversion pulse slot (%.3f ms). ', ...
        'Consider increasing the echo time or reducing the corresponding length factors.'], ...
        (sliceRephaseEnd - (Seq.Slice(2).CenterOfPulse - Seq.Slice(2).GradLength/2 + Seq.tEcho))*1e3);
    end
  else
    if sliceRephaseEnd + Seq.Slice(1).tEC > Seq.Slice(2).CenterOfPulse + ...
        Seq.Slice(1).Pulse.MaxLength/2 + Seq.Slice(1).Pulse.CenterOffset + Seq.tEcho
      warning('PD:sequence_Spin_Echo:SliceRephaseTooLate', ...
        ['The slice rephase slot might be colliding with the inversion pulse slot (%.3f ms). ', ...
        'Consider increasing the echo time or reducing the corresponding length factors.'], ...
        (sliceRephaseEnd + Seq.Slice(1).tEC - (Seq.Slice(2).CenterOfPulse + ...
        Seq.Slice(1).Pulse.MaxLength/2 + Seq.Slice(1).Pulse.CenterOffset + Seq.tEcho))*1e3);
    end
  end


  %% merge all elements into the complete pulse program
  % add AQ
  AQ = add_AQ(AQ, Seq.Read(1).AQ);

  if Seq.CorrectPhase > 0
    AQ = add_AQ(AQ, Seq.Read(5).AQ);
  end

  % add TX
  if strcmp(Seq.LoopName{Seq.LoopNameCount}, 'CRLoop')
    Seq.Slice(1).TX.Amplitude = Seq.Slice(1).TX.Amplitude.*0;
    Seq.Slice(2).TX.Amplitude = Seq.Slice(2).TX.Amplitude.*0;
  end

  TX = add_TX(TX, Seq.Slice(1).TX);
  TX = add_TX(TX, Seq.Slice(2).TX);

  % Gradients for slice selection
  if Seq.AQSlice(1).thickness < 1
    Grad = add_Grad(Grad, Seq.Slice(1).Grad);
    Grad = add_Grad(Grad, Seq.Slice(1).GradRephase);
  end

  if Seq.AQSlice(1).thicknessInversion<1
    if Seq.AQSlice(1).inversionGradDephase
      Grad = add_Grad(Grad, Seq.Slice(2).GradDephase);
      if Seq.CorrectPhase < 1
        Grad = add_Grad(Grad, Seq.Slice(3).GradRephase);
        Grad = add_Grad(Grad, Seq.Slice(3).GradDephase);
      end
    end
    Grad = add_Grad(Grad, Seq.Slice(2).Grad);
    if Seq.AQSlice(1).inversionGradRephase
      % FIXME: Why is this a different condition to Seq.AQSlice(1).inversionGradDephase
      Grad = add_Grad(Grad, Seq.Slice(2).GradRephase);
    end
  end

  % Gradients
  if Seq.SteadyState_PreShots180SpoilerGradFirst
    % Spoilers to de-phase signal before first inversion pre-shot (without actual acquisition)
    Grad = add_Grad(Grad, Seq.Phase(10).GradDephase);
    Grad = add_Grad(Grad, Seq.Phase(11).GradDephase);
    Grad = add_Grad(Grad, Seq.Phase(12).GradDephase);
  end
  % spoilers to dephase signal before inversion (aligned with read rephase)
  Grad = add_Grad(Grad, Seq.Phase(4).GradDephase);
  Grad = add_Grad(Grad, Seq.Phase(5).GradDephase);
  Grad = add_Grad(Grad, Seq.Phase(6).GradDephase);
  % phase encoding before acquisition (aligned with read dephase)
  Grad = add_Grad(Grad, Seq.Phase(1).GradDephase);
  Grad = add_Grad(Grad, Seq.Phase(2).GradDephase);
  Grad = add_Grad(Grad, Seq.Phase(3).GradDephase);
  % spoilers to rephase signal after inversion (aligned with read dephase)
  Grad = add_Grad(Grad, Seq.Phase(4).GradRephase);
  Grad = add_Grad(Grad, Seq.Phase(5).GradRephase);
  Grad = add_Grad(Grad, Seq.Phase(6).GradRephase);
  % read out gradient
  Grad = add_Grad(Grad, Seq.Read(1).Grad);
  if ~Seq.DephaseBefore180
    % read dephase
    Grad = add_Grad(Grad, Seq.Read(1).GradDephase);
    if Seq.AQSlice(1).RephaseLengthFactor ~= 0
      % read rephase
      Grad = add_Grad(Grad, Seq.Read(1).GradRephase);
    end
  end
  % read out at pre and post shots
  Grad = add_Grad(Grad, Seq.Read(3).Grad);
  if Seq.SteadyState_PreShots180ReadGrad
    % read out gradient blocks at inversion pre-shots (without actual acquisition)
    Grad = add_Grad(Grad, Seq.Read(2).Grad);
    if ~Seq.DephaseBefore180
      Grad = add_Grad(Grad, Seq.Read(2).GradDephase);
      if Seq.AQSlice(1).RephaseLengthFactor ~= 0
        Grad = add_Grad(Grad, Seq.Read(2).GradRephase);
      end
    end
  end
  if ~Seq.DephaseBefore180
    % read dephase and rephase at pre and post shots
    Grad = add_Grad(Grad, Seq.Read(3).GradDephase);
    if Seq.AQSlice(1).RephaseLengthFactor ~= 0
      Grad = add_Grad(Grad, Seq.Read(3).GradRephase);
    end
  else
    Grad = add_Grad(Grad, Seq.Read(4).GradDephase);
  end
  if Seq.AQSlice(1).RephaseLengthFactor ~= 0
    % phase rephase (aligned with read rephase)
    Grad = add_Grad(Grad, Seq.Phase(1).GradRephase);
    Grad = add_Grad(Grad, Seq.Phase(2).GradRephase);
    Grad = add_Grad(Grad, Seq.Phase(3).GradRephase);
  end
  % spoilers to dephase signal at end of echo train
  Grad = add_Grad(Grad, Seq.Phase(7).GradDephase);
  Grad = add_Grad(Grad, Seq.Phase(8).GradDephase);
  Grad = add_Grad(Grad, Seq.Phase(9).GradDephase);
  if Seq.SteadyState_PreShots180SpoilerGradFirst
    % Spoilers to re-phase signal before first inversion pre-shot with acquisition
    Grad = add_Grad(Grad, Seq.Phase(10).GradRephase);
    Grad = add_Grad(Grad, Seq.Phase(11).GradRephase);
    Grad = add_Grad(Grad, Seq.Phase(12).GradRephase);
  end

  if Seq.CorrectPhase > 0 && ...
      (Seq.SteadyState_PreShots180SpoilerGrad && ~Seq.SteadyState_PreShots180SpoilerGradFirst)
    % move spoilers in tReps with excitation pulse
    Grad = add_Grad(Grad, Seq.Phase(13).GradRephase);
    Grad = add_Grad(Grad, Seq.Phase(14).GradRephase);
    Grad = add_Grad(Grad, Seq.Phase(15).GradRephase);
    Grad = add_Grad(Grad, Seq.Phase(13).GradDephase);
    Grad = add_Grad(Grad, Seq.Phase(14).GradDephase);
    Grad = add_Grad(Grad, Seq.Phase(15).GradDephase);
  end

  % % extra Gradients
  % for t=10+1:length(Seq.Phase)
  %   Grad=add_Grad(Grad, Seq.Phase(t).GradDephase);
  %   Grad=add_Grad(Grad, Seq.Phase(t).GradRephase);
  % end
  % for t=3+1:length(Seq.Read)
  %   Grad=add_Grad(Grad, Seq.Read(t).Grad);
  % end
  % for t=2+1:length(Seq.Slice)
  %   Grad=add_Grad(Grad, Seq.Slice(t).GradDephase);
  %   Grad=add_Grad(Grad, Seq.Slice(t).Grad);
  %   Grad=add_Grad(Grad, Seq.Slice(t).GradRephase);
  % end


  if Seq.SingletRep == 1
    %% combine tReps of each echo train into singular logical tReps

    % d90 = Seq.SteadyState_PreShots180+Seq.SteadyState_PostShots180+3;
    n90 = numel(Seq.P90tReps);
    Seq.singletRepStart = cumsum([0, Seq.tRep(1:end-1)]);
    temptRepEnd = cumsum(Seq.tRep);
    temptRepCombine = temptRepEnd([Seq.P90tReps(2:end),end+1]-1) - Seq.singletRepStart(Seq.P90tReps);
    % temptRepCombineStart = cumsum([0,temptRepCombine(1:end-1)]);
    temptRepAdd = zeros(size(Seq.tRep));
    temptRepAdd(Seq.P90tReps) = temptRepAdd(Seq.P90tReps) - [0, temptRepCombine(1:end-1)];
    temptRepAdd = Seq.singletRepStart + cumsum(temptRepAdd);
    Seq.tRep = temptRepCombine;

    Seq.P90tReps = 1:n90;
    Seq.PreShotstReps = 1:double(Seq.SteadyState_PreShots90);  % used anywhere?
    Seq.P180tReps = 1:n90;
    Seq.kLinestReps = repmat(Seq.SteadyState_PreShots90+1:n90-Seq.SteadyState_PostShots90, [size(Seq.kLinestReps,1) 1]);
    Seq.P180tRepsPrePost = setdiff(Seq.P180tReps(:), Seq.kLinestReps(1,:));
    Seq.PostShotstReps = n90-Seq.SteadyState_PostShots90:n90;  % used anywhere?
    Seq.PausetReps = 1:n90;


    % Seq.tRepTurboBlockPhaseExcitation=Seq.tRepTurboBlockPhase(1,:);
    % Seq.tRepTurboBlockPhaseInversion=Seq.tRepTurboBlockPhase(2,:);
    % Seq.tRepTurboBlockPhaseRead=Seq.tRepTurboBlockPhaseRead(2,:);

    % Seq.tRep=(Seq.tEcho*Seq.SteadyState_PreShots180+Seq.SteadyState_PostShots180+2+Seq.AQSlice(1).TurboBreak)*ones(1,numel(Seq.P90tReps));

    % sort pulse program elements by time
    for t = 1:numel(Grad)
      if ~isempty(Grad(t).Time)
        [Grad(t).Time, I] = sort(reshape(bsxfun(@plus, Grad(t).Time, temptRepAdd), [], n90), 1);
        I = bsxfun(@plus, I, 0:size(I,1):size(I,1)*(size(I,2)-1));
        Grad(t).Amp = reshape(Grad(t).Amp, [], n90);
        Grad(t).Amp = Grad(t).Amp(I);
        [row,~] = find(sum(~isnan(Grad(t).Time),2), 1, 'last');
        Grad(t).Time = Grad(t).Time(1:max(row),:);
        Grad(t).Amp = Grad(t).Amp(1:max(row),:);
      end
    end
    for t = 1:numel(TX)
      if ~isempty(TX(t).Start)
        [TX(t).Start, I] = sort(reshape(bsxfun(@plus, TX(t).Start, temptRepAdd), [], n90), 1);
        I = bsxfun(@plus, I, 0:size(I,1):size(I,1)*(size(I,2)-1));
        TX(t).Duration = reshape(TX(t).Duration, [], n90);
        TX(t).Duration = TX(t).Duration(I);
        TX(t).Amplitude = reshape(TX(t).Amplitude, [], n90);
        TX(t).Amplitude = TX(t).Amplitude(I);
        TX(t).Frequency = reshape(TX(t).Frequency, [], n90);
        TX(t).Frequency = TX(t).Frequency(I);
        TX(t).Phase = reshape(TX(t).Phase, [], n90);
        TX(t).Phase = TX(t).Phase(I);
      end
    end
    for t = 1:numel(AQ)
      if ~isempty(AQ(t).Start)
        [AQ(t).Start, I] = sort(reshape(bsxfun(@plus, AQ(t).Start, temptRepAdd), [], n90), 1);
        I = bsxfun(@plus, I, 0:size(I,1):size(I,1)*(size(I,2)-1));
        AQ(t).fSample = reshape(AQ(t).fSample, [], n90);
        AQ(t).fSample = AQ(t).fSample(I);
        AQ(t).Frequency = reshape(AQ(t).Frequency, [], n90);
        AQ(t).Frequency = AQ(t).Frequency(I);
        AQ(t).Phase = reshape(AQ(t).Phase, [], n90);
        AQ(t).Phase = AQ(t).Phase(I);
        AQ(t).nSamples = reshape(AQ(t).nSamples, [], n90);
        AQ(t).nSamples = AQ(t).nSamples(I);
        AQ(t).SamplingFactor = reshape(AQ(t).SamplingFactor, [], n90);
        AQ(t).SamplingFactor = AQ(t).SamplingFactor(I);
        AQ(t).GetData = reshape(AQ(t).GetData, [], n90);
        AQ(t).GetData = any(AQ(t).GetData, 1);
        AQ(t).ResetPhases = reshape(AQ(t).ResetPhases, [], n90);
        AQ(t).ResetPhases = any(AQ(t).ResetPhases, 1);
      end
    end
    if ~isemptyfield(Seq, 'DigitalIO')
      for t = 1:numel(Seq.DigitalIO)
        if all(isnan(Seq.DigitalIO(t).SetTime(:)))
          Seq.DigitalIO(t).SetTime = [];
          Seq.DigitalIO(t).SetValue = [];
          Seq.DigitalIO(t).Repeat = [];
        end
        if ~isempty(Seq.DigitalIO(t).SetTime)
          [Seq.DigitalIO(t).SetTime, I] = sort(reshape(bsxfun(@plus, Seq.DigitalIO(t).SetTime, temptRepAdd), [], n90), 1);
          I = bsxfun(@plus, I, 0:size(I,1):size(I,1)*(size(I,2)-1));
          Seq.DigitalIO(t).SetValue = reshape(Seq.DigitalIO(t).SetValue, [], n90);
          Seq.DigitalIO(t).SetValue = Seq.DigitalIO(t).SetValue(I);
          if ~isemptyfield(Seq.DigitalIO(t), 'Repeat')
            Seq.DigitalIO(t).Repeat = reshape(Seq.DigitalIO(t).Repeat, [], n90);
            Seq.DigitalIO(t).Repeat = all(Seq.DigitalIO(t).Repeat, 1);
          end
        end
      end
    end
  end


  % map AQ windows in tReps to k-lines of images
  Seq.AQSlice(1).UsetRep = zeros(size(Seq.kLinesImages));
  Seq.AQSlice(1).UseAQWindow = ones(size(Seq.kLinesImages));
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
  if Seq.DephaseBefore180
    % In this case, the phase rephase and dephase pulses between echoes are
    % omitted. Thus, every other echo in the echo train has the opposite phase
    % encoding. DephaseBefore180 cannot be combined with TurboFactor. So we know
    % each echo in the train corresponds to the next image.

    % FIXME: Setting DephaseBefore180 with TurboFactor>1 or with nImages>1 or
    % with oddEvenEchoes is prevented with an error above. When does this ever
    % do anything useful?

    % FIXME: This probably needs changes once we implement support for those
    % combinations. For now the only correction we allow in combination with
    % DephaseBefore180 is phaseCycling afaict. That requires to not flip the
    % effective encoding order along the iParts dimension. (That might be
    % different for oddEvenEchoes.)

    Seq.AQSlice(1).UsetRep(:,2:2:end,:) = ...
      flipud(Seq.AQSlice(1).UsetRep(:,2:2:end,:));
    if mod(size(Seq.AQSlice(1).UsetRep, 1), 2) == 0
      % shift such that k-line without encoding stays at center
      Seq.AQSlice(1).UsetRep(:,2:2:end,:) = ...
        Seq.AQSlice(1).UsetRep([end,1:end-1],2:2:end,:);
    end
    if Seq.SingletRep
      Seq.AQSlice(1).UseAQWindow(:,2:2:end,:) = ...
        flipud(Seq.AQSlice(1).UseAQWindow(:,2:2:end,:));
      if mod(size(Seq.AQSlice(1).UseAQWindow, 1), 2) == 0
        % shift such that k-line without encoding stays at center
        Seq.AQSlice(1).UseAQWindow(:,2:2:end,:) = ...
          Seq.AQSlice(1).UseAQWindow([end,1:end-1],2:2:end,:);
      end
    end
  end

  if Seq.SingletRep
    % wrap AQ windows to number of AQ windows

    % all Turbo-Blocks have the same length
    Seq.AQSlice(1).UseAQWindow = mod(Seq.AQSlice(1).UseAQWindow, ...
      (Seq.AQSlice(1).oddEvenEchoes+1)*Seq.AQSlice(1).TurboFactor*Seq.AQSlice(1).nImages);
    Seq.AQSlice(1).UseAQWindow(Seq.AQSlice(1).UseAQWindow==0) = ...
      (Seq.AQSlice(1).oddEvenEchoes+1)*Seq.AQSlice(1).TurboFactor*Seq.AQSlice(1).nImages;

    if Seq.CorrectPhase > 0
      Seq.AQSlice(1).UseAQWindow = Seq.AQSlice(1).UseAQWindow + 1;

      % adjust tReps of frequency tracking windows
      Seq.AQSlice(2).UsetRep = floor(Seq.AQSlice(2).UsetRep / (numel(Seq.singletRepStart) / n90)) + 1;
    end
  end

  if ~any([Seq.PreProcessSequence, ...
           Seq.StartSequence, ...
           Seq.PollPPGfast, ...
           Seq.GetRawData, ...
           Seq.PostProcessSequence])
    % early break if only preparing
    % FIXME: This only works for a single loop. Selecting a loop (apart from
    % type 'normal') doesn't generally work yet.
    SeqLoop = Seq;
    if isobject(HW)
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
    % pre-process sequence to get PPG
    tStartSequence = Seq.StartSequence;
    Seq.StartSequence = 0;
    Seq.PollPPGfast = 0;
    Seq.GetRawData = 0;
    Seq.PostProcessSequence = 0;

    if tStartSequence
      Seq.Reinitialize = double((Seq.LoopNameCount==1) || ~Seq.LoopsBreakExactly);
      [ ~, SeqOut] = set_sequence(HW, Seq, AQ, TX, Grad);  % PreProcessSequence sequence
    else
      [ ~, SeqLoop] = set_sequence(HW, Seq, AQ, TX, Grad); % PreProcessSequence sequence only
      SeqLoop.TimeToNextSequence = SeqLoop.LoopsBreak - (SeqLoop.tRepEnd - SeqLoop.tEcho);
      return;
    end
  else
    SeqOut = Seq;
  end

  %% adjust timing (in loops)
  if ~isempty(SeqOut.LoopsRepetitionTime) && isempty(SeqOut.LoopsBreak)
    init.Seq.LoopsBreak = SeqOut.LoopsRepetitionTime - SeqOut.SequenceTime;
    SeqOut.LoopsBreak = init.Seq.LoopsBreak;
  end

  if SeqOut.LoopsBreakExactly
    if SeqOut.LoopNameCount == 1
      SeqOut.Reinitialize = 1;
      if numel(SeqOut.LoopName) > 1
        SeqOut.TimeToNextSequence = SeqOut.LoopsBreak - (SeqOut.tRepEnd - SeqOut.tEcho);
        SeqOut.TimeFromLastSequence = [];
      else
        SeqOut.TimeToNextSequence = [];
        SeqOut.TimeFromLastSequence = [];
      end
    elseif (SeqOut.LoopNameCount>=2) && (SeqOut.LoopNameCount~=numel(SeqOut.LoopName))
      SeqOut.Reinitialize=0;
      SeqOut.TimeFromLastSequence = SeqOut.LoopsBreak - (SeqOut.tRepEnd - SeqOut.tEcho);
      SeqOut.TimeToNextSequence   = SeqOut.LoopsBreak - (SeqOut.tRepEnd - SeqOut.tEcho);
    elseif SeqOut.LoopNameCount == numel(SeqOut.LoopName)
      SeqOut.Reinitialize=0;
      SeqOut.TimeFromLastSequence = SeqOut.LoopsBreak - (SeqOut.tRepEnd - SeqOut.tEcho);
      SeqOut.TimeToNextSequence = [];
    end
  end

  if isempty(SeqOut.LoopsBreak), tb=1; else tb=SeqOut.LoopsBreak; end
  if strcmp(SeqOut.LoopName{SeqOut.LoopNameCount}, 'normal') && (SeqOut.Loops > 1)
    if isempty(SeqOut.LoopsRepetitionTime)
      fprintf('Time to run remaining loops % 10.1f sec.\n', ...
        SeqOut.SequenceTime * (SeqOut.Loops-Loop+1) + tb*(SeqOut.Loops-Loop+1));
    else
      fprintf('Time to run remaining loops % 10.1f sec.\n', ...
        SeqOut.LoopsRepetitionTime * (SeqOut.Loops-Loop+1) + SeqOut.SequenceTime*double(SeqOut.Loops==Loop));
    end
  end
  clear tb

  %% actual measurement
  SeqOut.PreProcessSequence = 0;
  SeqOut.StartSequence = 1;
  SeqOut.PollPPGfast = 1;
  SeqOut.GetRawData = 1;
  SeqOut.PostProcessSequence = 1;
  [~, SeqOut, data, ~] = set_sequence(HW, SeqOut, AQ, TX, Grad);
  % talker.mySequency.exctractArrayToFile(talker.mySequency.getCommandArray,'test.txt');

  if ~SeqOut.LoopsBreakExactly && ~isempty(SeqOut.LoopsBreak)
    init.Seq.StartSequenceTime = SeqOut.StartSequenceTime + SeqOut.SequenceTime + SeqOut.LoopsBreak;
  end

  if SeqOut.LoopsBreakExactly
    init.Seq.EndTimeFPGA = SeqOut.EndTimeFPGA ...
      + SeqOut.TR_Error ./ [SeqOut.HW.MMRT(:).fSystem] * 4;
    init.Seq.LoopCountStart = SeqOut.LoopCountEnd + SeqOut.TR_Error;
    init.Seq.TimeToNextSequence = SeqOut.TimeFromLastSequence;
    init.Seq.TimeFromLastSequence = SeqOut.TimeToNextSequence;
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


  % re-construct image and k-space from measured data
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

  iAQ = find([SeqOut.AQ(:).Device] == SeqOut.AQSlice(1).iDevice, 1, 'first');
  iData = find([data(:).device] == SeqOut.AQ(iAQ).Device & [data(:).channel] == SeqOut.AQ(iAQ).Channel);  % Might be multiple channels
  iPrimData = iData([data(iData).nucleusX] == false);
  iSecData = iData([data(iData).nucleusX] == true);

  % find k-line without phase encoding (in image indexing)
  iNoPhase = floor(SeqOut.AQSlice(1).nPhase(3)*SeqOut.AQSlice(1).PhaseOS(3)/2)*prod(SeqOut.AQSlice(1).nPhase(1:2).*SeqOut.AQSlice(1).PhaseOS(1:2)) + ...
    floor(SeqOut.AQSlice(1).nPhase(2)*SeqOut.AQSlice(1).PhaseOS(2)/2)*SeqOut.AQSlice(1).nPhase(1)*SeqOut.AQSlice(1).PhaseOS(1) + ...
    floor(SeqOut.AQSlice(1).nPhase(1)*SeqOut.AQSlice(1).PhaseOS(1)/2)+1;
  tRepNoPhase = SeqOut.AQSlice(1).UsetRep(iNoPhase,:);  % in tRep indexing
  % FIXME: Consider case when k-line through k-space center is skipped?
  if any(tRepNoPhase <= numel(Seq.tRep) & tRepNoPhase > 0)
    % find corresponding excitation in tRep indexing
    tRepExcitation = repmat(SeqOut.P90tReps, [size(tRepNoPhase, 2) 1]);
    tRepExcitation = tRepExcitation(sub2ind(size(tRepExcitation), 1:size(tRepExcitation,1), sum(bsxfun(@le, SeqOut.P90tReps, tRepNoPhase.'), 2).')).';

    tRepCum = [0, cumsum(SeqOut.tRep(1:end-1))];
    data(iPrimData).tImageZ = ...
      reshape((tRepCum(tRepNoPhase) + SeqOut.Read(1).CenterOfReadout) ...  % centers of acquisition windows (time_all)
              - (tRepCum(tRepExcitation) + SeqOut.Slice(1).CenterOfPulse), ...  % centers of the respective excitation pulses (time_all)
              [], size(SeqOut.AQSlice(1).UsetRep, 3));

    if SeqOut.SingletRep
      temptRepAdd = reshape(temptRepAdd, [], n90);
      % add time of AQ window in tRep
      % Attention: This assumes that the excitation pulse wasn't moved within the tRep
      data(iPrimData).tImageZ = data(iPrimData).tImageZ + ...
        reshape(temptRepAdd(sub2ind(size(temptRepAdd), SeqOut.AQSlice(1).UseAQWindow(iNoPhase,:), tRepNoPhase)), ...
        [], size(SeqOut.AQSlice(1).UsetRep, 3));
    end
  end


  if Loop == 1
    SeqLoop = SeqOut;
    SeqLoop.data = data(iPrimData);
    SeqLoop.AQSlice(1).raiseFigures = 1;
  end
  SeqLoop.data.StartSequenceTime = SeqOut.StartSequenceTime;
  SeqLoop.data.fCenter = SeqOut.HW.fLarmor;
  if SeqLoop.AQSlice(1).dualNuclearImage
    if Loop == 1
      SeqLoop.dataX = data(iSecData);
    end
    SeqLoop.dataX.fCenter = SeqOut.HW.fLarmorX;
  end


  switch SeqOut.LoopName{SeqOut.LoopNameCount}
    case 'CSRLoopPlus'
      % correct slice rephase
      % in first loop just store the data
      SeqLoop.dataLoop(1) = deal(data(iPrimData));  % FIXME: Why "deal"?

    case 'CSRLoopMinus'
      % correct slice rephase
      SeqLoop.dataLoop(2) = deal(data(iPrimData));  % FIXME: Why "deal"?

      % don't plot unless requested for this correction
      SeqLoop.AQSlice(1).plotB0ppm=0;
      SeqLoop.AQSlice(1).plotB0Hz=0;
      SeqLoop.AQSlice(1).plotB0HzPhase=0;
      SeqLoop.AQSlice(1).plotB0PpmPhase=0;
      if SeqOut.CorrectSliceRephasePlot
        SeqLoop.AQSlice(1).plotImage=1;
        SeqLoop.AQSlice(1).plotkSpace=1;
        SeqLoop.AQSlice(1).plotPhase=1;
      else
        SeqLoop.AQSlice(1).plotImage=0;
        SeqLoop.AQSlice(1).plotkSpace=0;
      end
      [SeqLoop.data, SeqLoop.AQSlice(1)] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));

      if SeqOut.CorrectSliceRephasePlot
        hf = figure(201); clf(hf);
        ax(1) = subplot(2,1,1, 'Parent', hf);
        plot(ax(1), SeqLoop.data.Ticks(1).Read, abs(SeqLoop.dataLoop(1).Image), ...
          SeqLoop.data.Ticks(1).Read, abs(SeqLoop.dataLoop(2).Image));
        ax(2) = subplot(2,1,2, 'Parent', hf);
        plot(ax(2), SeqLoop.data.Ticks(1).Read, unwrap2Dmiddle(angle(SeqLoop.dataLoop(1).Image)), ...
          SeqLoop.data.Ticks(1).Read, unwrap2Dmiddle(angle(SeqLoop.dataLoop(2).Image)));
        ylim(ax(2), [-2*pi, 2*pi]);
        linkaxes(ax, 'x');
      end

      % select ROI
      if ~isemptyfield(SeqLoop.AQSlice(1), 'SliceLimMinMax')
        SeqLoop.AQSlice(1).SliceLim(1) = max(SeqLoop.AQSlice(1).SliceLimMinMax(1), ...
          SeqLoop.AQSlice(1).Center2OriginImage(1) - ...
          min(SeqLoop.AQSlice(1).thickness, SeqLoop.AQSlice(1).thicknessInversion)/2);
        SeqLoop.AQSlice(1).SliceLim(2) = min(SeqLoop.AQSlice(1).SliceLimMinMax(2), ...
          SeqLoop.AQSlice(1).Center2OriginImage(1) + ...
          min(SeqLoop.AQSlice(1).thickness, SeqLoop.AQSlice(1).thicknessInversion)/2);
        SeqLoop.AQSlice(1).SliceLimI(1) = ...
          find(SeqLoop.data.Ticks(1).Read > SeqLoop.AQSlice(1).SliceLim(1), 1, 'first');
        SeqLoop.AQSlice(1).SliceLimI(2) = ...
          find(SeqLoop.data.Ticks(1).Read < SeqLoop.AQSlice(1).SliceLim(2), 1, 'last');
        SeqLoop.AQSlice(1).SliceLimRoi = zeros(size(SeqLoop.dataLoop(1).Image));
        SeqLoop.AQSlice(1).SliceLimRoi(SeqLoop.AQSlice(1).SliceLim(1):SeqLoop.AQSlice(1).SliceLimI(2)) = 1;
      end
      if isemptyfield(SeqLoop.AQSlice(1), 'SliceLimRoi')
        SeqLoop.AQSlice(1).SliceLimRoi = and(abs(SeqLoop.dataLoop(1).Image) > (max(abs(SeqLoop.dataLoop(1).Image))/4), ...
          abs(SeqLoop.dataLoop(2).Image) > (max(abs(SeqLoop.dataLoop(2).Image))/4));
      end

      SeqLoop.AQSlice(1).SliceLimRoi=conv(double(SeqLoop.AQSlice(1).SliceLimRoi),[1;1;1;1],'same')>=3;
      SeqLoop.AQSlice(1).SliceLimDiffRoi=SeqLoop.AQSlice(1).SliceLimRoi(1:end-1);

      % calculate slope difference of image phase in read direction
      SeqLoop.data.DiffTicksSlice = diff(SeqLoop.data.Ticks(1).Read);
      SeqLoop.data.MeanDiffTicksSlice = mean(SeqLoop.data.DiffTicksSlice(SeqLoop.AQSlice(1).SliceLimDiffRoi));
      SeqLoop.dataLoop(1).DiffAngleSlice = diff(unwrap(angle(SeqLoop.dataLoop(1).Image)));
      SeqLoop.dataLoop(2).DiffAngleSlice = diff(unwrap(angle(SeqLoop.dataLoop(2).Image)));
      SeqLoop.dataLoop(1).SlopeAngleSlice = mean(SeqLoop.dataLoop(1).DiffAngleSlice(SeqLoop.AQSlice(1).SliceLimDiffRoi)) ./ SeqLoop.data.MeanDiffTicksSlice;
      SeqLoop.dataLoop(2).SlopeAngleSlice = mean(SeqLoop.dataLoop(2).DiffAngleSlice(SeqLoop.AQSlice(1).SliceLimDiffRoi)) ./ SeqLoop.data.MeanDiffTicksSlice;
      SeqLoop.data.MeanSlopeAngleSlice = mean([SeqLoop.dataLoop(2).SlopeAngleSlice, SeqLoop.dataLoop(1).SlopeAngleSlice]);
      SeqLoop.data.DiffSlopeAngleSlice = diff([SeqLoop.dataLoop(2).SlopeAngleSlice, SeqLoop.dataLoop(1).SlopeAngleSlice]);
      % calculate offset to slice rephase encoder to compensate slope
      SeqLoop.data.SliceReadGradTimeIntegralOffset = ...
        SeqLoop.data.MeanSlopeAngleSlice / SeqLoop.HW.GammaDef ...
        .* (-1).^SeqLoop.SteadyState_PreShots180;
      SeqLoop.data.SliceReadGradTimeIntegralOffset = 0;
      SeqLoop.data.SliceGradTimeIntegralRephaseOffset = ...
        SeqLoop.data.DiffSlopeAngleSlice / SeqLoop.HW.GammaDef ./ 2 ....
        .* (-1).^SeqLoop.SteadyState_PreShots180;
      % SeqLoop.data.SliceGradTimeIntegralRephaseOffset
      if isnan(SeqLoop.data.SliceGradTimeIntegralRephaseOffset)
        warning('isnan(SeqLoop.data.SliceGradTimeIntegralRephaseOffset) set to 0');
        SeqLoop.data.SliceGradTimeIntegralRephaseOffset = 0;
      end
      if SeqOut.CorrectSliceRephasePlot
        hf = figure(202); clf(hf);
        ax(1) = subplot(2,1,1, 'Parent', hf);
        plot(ax(1), SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).SliceLimRoi), abs(SeqLoop.dataLoop(1).Image(SeqLoop.AQSlice(1).SliceLimRoi)), ...
          SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).SliceLimRoi), abs(SeqLoop.dataLoop(2).Image(SeqLoop.AQSlice(1).SliceLimRoi)));
        ax(2) = subplot(2,1,2, 'Parent', hf);
        plot(ax(2), SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).SliceLimRoi), angle(SeqLoop.dataLoop(1).Image(SeqLoop.AQSlice(1).SliceLimRoi)), ...
          SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).SliceLimRoi), angle(SeqLoop.dataLoop(2).Image(SeqLoop.AQSlice(1).SliceLimRoi)));
        ylim(ax(2), [-pi, pi]);
        linkaxes(ax, 'x');
        title(ax(2), ['SliceGradTimeIntegralRephaseOffset = ' num2str(SeqLoop.data.SliceGradTimeIntegralRephaseOffset) ' T s / m']);
      end

    case 'CPRLoopPlus'
      % correct phase rephase
      % in first loop just store the data
      SeqLoop.dataLoop(1) = deal(data(iPrimData));  % FIXME: Why "deal"?

    case 'CPRLoopMinus'
      % correct phase rephase
      SeqLoop.dataLoop(2) = deal(data(iPrimData));  % FIXME: Why "deal"?

      % don't plot unless requested for this correction
      SeqLoop.AQSlice(1).plotB0ppm=0;
      SeqLoop.AQSlice(1).plotB0Hz=0;
      SeqLoop.AQSlice(1).plotB0HzPhase=0;
      SeqLoop.AQSlice(1).plotB0PpmPhase=0;
      if SeqOut.CorrectPhaseRephasePlot
        SeqLoop.AQSlice(1).plotImage=1;
        SeqLoop.AQSlice(1).plotkSpace=1;
        SeqLoop.AQSlice(1).plotPhase=1;
      else
        SeqLoop.AQSlice(1).plotImage=0;
        SeqLoop.AQSlice(1).plotkSpace=0;
      end
      [SeqLoop.data, SeqLoop.AQSlice(1)] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));

      if SeqOut.CorrectPhaseRephasePlot
        hf = figure(201); clf(hf);
        ax(1) = subplot(2,1,1, 'Parent', hf);
        plot(ax(1), SeqLoop.data.Ticks(1).Read, abs(SeqLoop.dataLoop(1).Image), ...
          SeqLoop.data.Ticks(1).Read, abs(SeqLoop.dataLoop(2).Image));
        ax(2) = subplot(2,1,2, 'Parent', hf);
        plot(ax(2), SeqLoop.data.Ticks(1).Read, unwrap2Dmiddle(angle(SeqLoop.dataLoop(1).Image)), ...
          SeqLoop.data.Ticks(1).Read, unwrap2Dmiddle(angle(SeqLoop.dataLoop(2).Image)));
        ylim(ax(2), [-2*pi, 2*pi]);
        linkaxes(ax, 'x');
      end

      % select ROI
      if ~isemptyfield(SeqLoop.AQSlice(1), 'PhaseLimMinMax')
        SeqLoop.AQSlice(1).PhaseLim(1) = ...
          max(SeqLoop.AQSlice(1).PhaseLimMinMax(1), SeqLoop.AQSlice(1).Center2OriginImage(2)-SeqLoop.AQSlice(1).sizePhase(2)/2);
        SeqLoop.AQSlice(1).PhaseLim(2) = ...
          min(SeqLoop.AQSlice(1).PhaseLimMinMax(2), SeqLoop.AQSlice(1).Center2OriginImage(2)+SeqLoop.AQSlice(1).sizePhase(2)/2);

        SeqLoop.AQSlice(1).PhaseLimI(1) = ...
          find(SeqLoop.data.Ticks(1).Read > SeqLoop.AQSlice(1).PhaseLim(1), 1, 'first');
        SeqLoop.AQSlice(1).PhaseLimI(2) = ...
          find(SeqLoop.data.Ticks(1).Read < SeqLoop.AQSlice(1).PhaseLim(2), 1, 'last');
        SeqLoop.AQSlice(1).PhaseLimRoi = zeros(size(SeqLoop.dataLoop(1).Image));
        SeqLoop.AQSlice(1).PhaseLimRoi(SeqLoop.AQSlice(1).PhaseLimI(1):SeqLoop.AQSlice(1).PhaseLimI(2)) = 1;
      end
      if isemptyfield(SeqLoop.AQSlice(1), 'PhaseLimRoi')
        SeqLoop.AQSlice(1).PhaseLimRoi = and(abs(SeqLoop.dataLoop(1).Image) > (max(abs(SeqLoop.dataLoop(1).Image))/4), ...
          abs(SeqLoop.dataLoop(2).Image) > (max(abs(SeqLoop.dataLoop(2).Image))/4));
      end

      SeqLoop.AQSlice(1).PhaseLimRoi=conv(double(SeqLoop.AQSlice(1).PhaseLimRoi),[1;1;1;1],'same')>=3;
      SeqLoop.AQSlice(1).PhaseLimDiffRoi=SeqLoop.AQSlice(1).PhaseLimRoi(1:end-1);

      % calculate slope difference of image phase in read direction
      SeqLoop.data.DiffTicksPhase = diff(SeqLoop.data.Ticks(1).Read);
      SeqLoop.data.MeanDiffTicksPhase = ...
        mean(SeqLoop.data.DiffTicksPhase(SeqLoop.AQSlice(1).PhaseLimDiffRoi));
      SeqLoop.dataLoop(1).DiffAnglePhase = diff(unwrap(angle(SeqLoop.dataLoop(1).Image)));
      SeqLoop.dataLoop(2).DiffAnglePhase = diff(unwrap(angle(SeqLoop.dataLoop(2).Image)));
      SeqLoop.dataLoop(1).SlopeAnglePhase = mean(SeqLoop.dataLoop(1).DiffAnglePhase(SeqLoop.AQSlice(1).PhaseLimDiffRoi)) ./ SeqLoop.data.MeanDiffTicksPhase;
      SeqLoop.dataLoop(2).SlopeAnglePhase = mean(SeqLoop.dataLoop(2).DiffAnglePhase(SeqLoop.AQSlice(1).PhaseLimDiffRoi)) ./ SeqLoop.data.MeanDiffTicksPhase;
      SeqLoop.data.MeanSlopeAnglePhase = mean([SeqLoop.dataLoop(2).SlopeAnglePhase, SeqLoop.dataLoop(1).SlopeAnglePhase]);
      SeqLoop.data.DiffSlopeAnglePhase = diff([SeqLoop.dataLoop(2).SlopeAnglePhase, SeqLoop.dataLoop(1).SlopeAnglePhase]);
      % calculate offset to phase rephase encoder to compensate slope
      SeqLoop.data.PhaseReadGradTimeIntegralOffset = ...
        SeqLoop.data.MeanSlopeAnglePhase / SeqLoop.HW.GammaDef ...
        .* (-1).^SeqLoop.SteadyState_PreShots180;
      SeqLoop.data.PhaseReadGradTimeIntegralOffset = 0;
      SeqLoop.data.PhaseGradTimeIntegralRephaseOffset = ...
        SeqLoop.data.DiffSlopeAnglePhase / SeqLoop.HW.GammaDef./2 ...
        .* (-1).^SeqLoop.SteadyState_PreShots180;
      % SeqLoop.data.PhaseGradTimeIntegralRephaseOffset
      if isnan(SeqLoop.data.PhaseGradTimeIntegralRephaseOffset);
        warning('isnan(SeqLoop.data.PhaseGradTimeIntegralRephaseOffset) set to 0');
        SeqLoop.data.PhaseGradTimeIntegralRephaseOffset=0;
      end
      if SeqOut.CorrectPhaseRephasePlot
        hf = figure(202); clf(hf);
        ax(1) = subplot(2,1,1, 'Parent', hf);
        plot(ax(1), SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).PhaseLimRoi), abs(SeqLoop.dataLoop(1).Image(SeqLoop.AQSlice(1).PhaseLimRoi)), ...
          SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).PhaseLimRoi), abs(SeqLoop.dataLoop(2).Image(SeqLoop.AQSlice(1).PhaseLimRoi)));
        ax(2) = subplot(2,1,2, 'Parent', hf);
        plot(ax(2), SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).PhaseLimRoi), angle(SeqLoop.dataLoop(1).Image(SeqLoop.AQSlice(1).PhaseLimRoi)), ...
          SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).PhaseLimRoi), angle(SeqLoop.dataLoop(2).Image(SeqLoop.AQSlice(1).PhaseLimRoi)));
        ylim(ax(2), [-pi, pi]);
        linkaxes(ax, 'x');
        title(ax(2), ['PhaseGradTimeIntegralRephaseOffset = ' num2str(SeqLoop.data.PhaseGradTimeIntegralRephaseOffset) ' T s / m']);
      end

    case 'CRRLoop'
      % correct read rephase

      % don't plot unless requested for this correction
      SeqLoop.AQSlice(1).plotB0ppm=0;
      SeqLoop.AQSlice(1).plotB0Hz=0;
      SeqLoop.AQSlice(1).plotB0HzPhase=0;
      SeqLoop.AQSlice(1).plotB0PpmPhase=0;
      if SeqOut.CorrectReadRephasePlot
        SeqLoop.AQSlice(1).plotImage=1;
        SeqLoop.AQSlice(1).plotkSpace=1;
        SeqLoop.AQSlice(1).plotPhase=1;
      else
        SeqLoop.AQSlice(1).plotImage=0;
        SeqLoop.AQSlice(1).plotkSpace=0;
      end
      [SeqLoop.data, SeqLoop.AQSlice(1)] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));

      if SeqOut.CorrectReadRephasePlot
        hf = figure(201); clf(hf);
        ax(1) = subplot(2,1,1, 'Parent', hf);
        plot(ax(1), SeqLoop.data.Ticks(1).Read, abs(SeqLoop.data.Image));
        grid(ax(1), 'on');
        title(ax(1), 'Read rephase');
        ax(2) = subplot(2,1,2, 'Parent', hf);
        plot(ax(2), SeqLoop.data.Ticks(1).Read, unwrap2Dmiddle(angle(SeqLoop.data.Image)));
        ylim(ax(2), [-2*pi, 2*pi]);
        grid(ax(2), 'on');
        linkaxes(ax, 'x');
      end

      % select ROI
      if ~isemptyfield(SeqLoop.AQSlice(1), 'ReadLimMinMax')
        SeqLoop.AQSlice(1).ReadLim(1) = ...
          max(SeqLoop.AQSlice(1).ReadLimMinMax(1), SeqLoop.AQSlice(1).Center2OriginImage(3)-SeqLoop.AQSlice(1).sizeRead/2);
        SeqLoop.AQSlice(1).ReadLim(2) = ...
          min(SeqLoop.AQSlice(1).ReadLimMinMax(2), SeqLoop.AQSlice(1).Center2OriginImage(3)+SeqLoop.AQSlice(1).sizeRead/2);

        SeqLoop.AQSlice(1).ReadLimI(1) = ...
          find(SeqLoop.data.Ticks(1).Read > SeqLoop.AQSlice(1).ReadLim(1), 1, 'first');
        SeqLoop.AQSlice(1).ReadLimI(2) = ...
          find(SeqLoop.data.Ticks(1).Read < SeqLoop.AQSlice(1).ReadLim(2), 1, 'last');
        SeqLoop.AQSlice(1).ReadLimRoi = zeros(size(SeqLoop.data.Image));
        SeqLoop.AQSlice(1).ReadLimRoi(SeqLoop.AQSlice(1).ReadLimI(1):SeqLoop.AQSlice(1).ReadLimI(2)) = 1;
      end
      if isemptyfield(SeqLoop.AQSlice(1), 'ReadLimRoi')
        SeqLoop.AQSlice(1).ReadLimRoi = abs(SeqLoop.data.Image)>(max(abs(SeqLoop.data.Image))/4);
      end

      SeqLoop.AQSlice(1).ReadLimRoi = conv(double(SeqLoop.AQSlice(1).ReadLimRoi), [1;1;1;1], 'same') >= 3;
      SeqLoop.AQSlice(1).ReadLimDiffRoi = SeqLoop.AQSlice(1).ReadLimRoi(1:end-1);

      % calculate slope of image phase in read direction
      SeqLoop.data.DiffTicksRead = diff(SeqLoop.data(1).Ticks(1).Read);
      SeqLoop.data.MeanDiffTicksRead = ...
        mean(SeqLoop.data.DiffTicksRead(SeqLoop.AQSlice(1).ReadLimDiffRoi));
      SeqLoop.data.DiffAngleRead = diff(unwrap(angle(SeqLoop.data(1).Image)));
      SeqLoop.data.SlopeAngleRead = ...
        mean(SeqLoop.data.DiffAngleRead(SeqLoop.AQSlice(1).ReadLimDiffRoi)) ./ SeqLoop.data.MeanDiffTicksRead;
      % calculate offset to read encoder to compensate slope
      SeqLoop.data.ReadGradTimeIntegralOffset = SeqLoop.data.SlopeAngleRead / SeqLoop.HW.GammaDef;
      % SeqLoop.data.ReadGradTimeIntegralOffset
      if isnan(SeqLoop.data.ReadGradTimeIntegralOffset);
        warning('isnan(SeqLoop.data.ReadGradTimeIntegralOffset) set to 0');
        SeqLoop.data.ReadGradTimeIntegralOffset=0;
      end
      if SeqOut.CorrectReadRephasePlot
        hf = figure(202); clf(hf);
        ax(1) = subplot(2,1,1, 'Parent', hf);
        plot(ax(1), SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).ReadLimRoi), ...
          abs(SeqLoop.data.Image(SeqLoop.AQSlice(1).ReadLimRoi)));
        grid(ax(1), 'on');
        title(ax(1), 'Read rephase (RoI)');
        ax(2) = subplot(2,1,2, 'Parent', hf);
        plot(ax(2), SeqLoop.data.Ticks(1).Read(SeqLoop.AQSlice(1).ReadLimRoi), ...
          angle(SeqLoop.data.Image(SeqLoop.AQSlice(1).ReadLimRoi)));
        ylim(ax(2), [-pi, pi]);
        grid(ax(2), 'on');
        linkaxes(ax, 'x');
        title(ax(2), ['ReadGradTimeIntegralOffset = ' num2str(SeqLoop.data.ReadGradTimeIntegralOffset) ' T s / m'])
      end

    case 'B0map_tEcho1'
      % measure (and correct for) B0 deviation map
      % in first loop just store the data
      SeqLoop.dataLoop(Loop) = deal(data(iPrimData));  % FIXME: Why "deal"?

    case 'B0map_tEcho2'
      % measure (and correct for) B0 deviation map
      SeqLoop.dataLoop(Loop) = deal(data(iPrimData));  % FIXME: Why "deal"?

      dataB0 = get_B0_map_from_image(HW, SeqLoop, SeqLoop.dataLoop(Loop-1), SeqLoop.dataLoop(Loop));

      % FIXME: Rotate B0 map to magnet coordinate system?

      if SeqOut.CorrectB0Read.WaitForSample
        disp('Waiting for sample exchange...');
        waitfor(msgbox('Please, change the sample and press "OK" to continue.'));
      end

    case 'normal'
      if SeqOut.CorrectPhase && isfield(SeqOut, 'Correct_tfoffset')
        CorrectPhase.tfOffset = SeqOut.Correct_tfoffset;
        CorrectPhase.fOffset = SeqOut.Correct_foffset;
        CorrectPhase.fOffsetStd = SeqOut.Correct_foffsetStd;
      end

      if SeqOut.LoopSaveAllData
        % save complete data structure for each averaging step
        if SeqOut.CorrectB0Read.Use
          [data(iPrimData), SeqOut.AQSlice(1)] = correct_read_B0(data(iPrimData), ...
            SeqOut.AQSlice(1), dataB0, max(sum(SeqOut.Read(1).GradAmp.^2,1).^0.5));
        end
        data(iPrimData).StartSequenceTime = SeqOut.StartSequenceTime;
        data(iPrimData).fCenter = SeqOut.HW.fLarmor;
        if SeqOut.CorrectPhase && exist('CorrectPhase', 'var')
          data(iPrimData).CorrectPhase = CorrectPhase;
        end
        SeqLoop.dataLoop(Loop) = deal(data(iPrimData));  % FIXME: Why "deal"?
      else
        % save only reduced sub-set for each averaging step
        SeqLoop.dataLoop(Loop).Image = single(data(iPrimData).Image);
        SeqLoop.dataLoop(Loop).StartSequenceTime = SeqOut.StartSequenceTime;
        SeqLoop.dataLoop(Loop).fCenter = SeqOut.HW.fLarmor;

        if SeqOut.CorrectPhase && exist('CorrectPhase', 'var')
          SeqLoop.dataLoop(Loop).CorrectPhase = CorrectPhase;
        end
      end

      clear('CorrectPhase');

      if SeqOut.LoopSaveAllSeq
        SeqLoop.SeqLoop(Loop) = SeqOut;
      end

      dataFieldname = {'data', 'dataX'};
      dataIdx = [iPrimData, iSecData];
      for iX = 1:numel(dataIdx)
        % average over loops
        SeqLoop.(dataFieldname{iX}).data = (SeqLoop.(dataFieldname{iX}).data*(Loop-1) + data(dataIdx(iX)).data)/Loop;
        if isfield(data, 'ImageZAll')  % SeqOut.AQSlice(1).nImages > 1
          % FIXME: Can we ever reach here?
          if ~isfield(SeqLoop.(dataFieldname{iX}), 'ImageZAll')
            SeqLoop.(dataFieldname{iX}).ImageZAll = data(dataIdx(iX)).ImageZAll;
          else
            SeqLoop.(dataFieldname{iX}).ImageZAll = (SeqLoop.(dataFieldname{iX}).ImageZAll*(Loop-1)+data(dataIdx(iX)).ImageZAll)/Loop;
          end
        end
        if isfield(data(dataIdx(iX)), 'fft1_dataCut')
          % FIXME: Can we ever reach here?
          if ~isfield(SeqLoop.(dataFieldname{iX}), 'fft1_dataCut')
            SeqLoop.(dataFieldname{iX}).fft1_dataCut = data(dataIdx(iX)).fft1_dataCut;
          else
            SeqLoop.(dataFieldname{iX}).fft1_dataCut = (SeqLoop.(dataFieldname{iX}).fft1_dataCut*(Loop-1) + data(dataIdx(iX)).fft1_dataCut)/Loop;
          end
        end
      end

      if SeqOut.LoopPlot || SeqOut.Loops==Loop  % always plot at last averaging step
        if SeqOut.LoopPlotAverages || (SeqOut.Loops==Loop && SeqOut.LoopPlotLastAverage)
          % plot the (current) averaged image
          try
            SeqLoop.data = get_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
          catch ME
            warning('PD:sequence_Flash:ImageReconstructionError', ...
              'An error occurred during image reconstruction. Trying to continue.\nError:\n%s', ...
              getReport(ME));
            continue;
          end
          if SeqOut.CorrectB0Read.Use
            if ~isfield(SeqOut.AQSlice(1), 'CorrectAmplitude')
              SeqOut.AQSlice(1).CorrectAmplitude = [];
            end
            [SeqLoop.data, SeqOut.AQSlice(1)] = correct_read_B0(SeqLoop.data, ...
              SeqOut.AQSlice(1), dataB0, max(sum(SeqOut.Read(1).GradAmp.^2,1).^0.5));
          end
          SeqLoop.data.RoI = [];   % re-calculate RoI with averaged data
          try
            [SeqLoop.data, SeqLoop.AQSlice(1)] = plot_kSpaceAndImage(SeqLoop.data, SeqLoop.AQSlice(1));
          catch ME
            warning('PD:sequence_Flash:ImageReconstructionError', ...
              'An error occurred while trying to display the results. Trying to continue.\nError:\n%s', ...
              getReport(ME));
            continue;
          end
          SeqLoop.AQSlice(1).raiseFigures = SeqOut.AQSlice(1).raiseFigures;
        elseif SeqOut.LoopPlotAll
          % plot the singular (averaging) step
          % [SeqLoop.data] = get_kSpaceAndImage(SeqLoop.dataLoop(Loop), SeqLoop.AQSlice(1));
          if SeqOut.CorrectB0Read.Use
            if ~isfield(SeqOut.AQSlice(1), 'CorrectAmplitude')
              SeqOut.AQSlice(1).CorrectAmplitude = [];
            end
            [SeqLoop.dataLoop(Loop), SeqOut.AQSlice(1)] = correct_read_B0(SeqLoop.dataLoop(Loop), ...
              SeqOut.AQSlice(1), dataB0, max(sum(SeqOut.Read(1).GradAmp.^2,1).^0.5));
          end
          try
            if SeqOut.Loops == Loop
              % store data from last loop (including RoI)
              [SeqLoop.data, SeqLoop.AQSlice(1)] = plot_kSpaceAndImage(SeqLoop.dataLoop(Loop), SeqLoop.AQSlice(1));
            else
              [~, SeqLoop.AQSlice(1)] = plot_kSpaceAndImage(SeqLoop.dataLoop(Loop), SeqLoop.AQSlice(1));
            end
          catch ME
            warning('PD:sequence_Flash:ImageReconstructionError', ...
              'An error occurred while trying to display the results. Trying to continue.\nError:\n%s', ...
              getReport(ME));
            continue;
          end
          SeqLoop.AQSlice(1).raiseFigures = SeqOut.AQSlice(1).raiseFigures;
        end
      end

    otherwise
      % do nothing

  end



  %%
  % if isfield(SeqLoop.data,'ImageCsi')
  % %%
  % % index=[2,5,3];
  % index=round(0.5.*[size(SeqLoop.data.ImageCsi,2),size(SeqLoop.data.ImageCsi,3),size(SeqLoop.data.ImageCsi,4)]);
  % figure(60)
  % plot(SeqLoop.data.ImageCsiFrequency(:,index(1),index(2),index(3))/SeqLoop.HW.fLarmor*1e6-1e6,abs(SeqLoop.data.ImageCsi(:,index(1),index(2),index(3)))*1e12)%
  % xlim([-100 100])
  % xlabel('Frequency offset in ppm')
  % ylabel('Amplitude in pT')
  %
  % %%
  % % index=[5,4,4];
  % index=round(0.5.*[size(SeqLoop.data.ImageCsiRawZero,2),size(SeqLoop.data.ImageCsiRawZero,3),size(SeqLoop.data.ImageCsiRawZero,4)]);
  % figure(61)
  % plot(repmat(SeqLoop.data.ImageCsiFrequencyZero(:,index(1),index(2),index(3))/SeqLoop.HW.fLarmor*1e6-1e6,1,3), ...
  %   [abs(SeqLoop.data.ImageCsiRawZero(:,index(1),index(2),index(3))), ...
  %   real(SeqLoop.data.ImageCsiRawZero(:,index(1),index(2),index(3))), ...
  %   imag(SeqLoop.data.ImageCsiRawZero(:,index(1),index(2),index(3)))])
  % xlim([-100 100])
  % xlabel('Frequency offset in ppm')
  % ylabel('Amplitude in pT')
  % %%
  % % index=[5,8,4];
  % index=round(0.5.*[size(SeqLoop.data.kSpaceOsRawFrequency,2),size(SeqLoop.data.kSpaceOsRawFrequency,3),size(SeqLoop.data.kSpaceOsRawFrequency,4)]);
  % figure(62)
  % plot(squeeze(SeqLoop.data.time_of_tRep(:,SeqLoop.AQSlice.UseAQWindow,SeqLoop.AQSlice.UsetRep(1)))-SeqLoop.tEcho/2, ...
  %   [abs(SeqLoop.data.kSpaceOsRawFrequency(:,index(1),index(2),index(3))), ...
  %   real(SeqLoop.data.kSpaceOsRawFrequency(:,index(1),index(2),index(3))), ...
  %   imag(SeqLoop.data.kSpaceOsRawFrequency(:,index(1),index(2),index(3)))])
  % xlabel('kSpace Point Signal ')
  % ylabel('Amplitude in pT')
  % %%
  % %%
  % index=[13,25,13];
  % % index=round(0.5.*[size(SeqLoop.data.kSpaceOsRawFrequency,2),size(SeqLoop.data.kSpaceOsRawFrequency,3),size(SeqLoop.data.kSpaceOsRawFrequency,4)]);
  % figure(62)
  % plot([abs(SeqLoop.data.kSpaceOsZ(:,index(1),index(2),index(3))), ...
  %   real(SeqLoop.data.kSpaceOsZ(:,index(1),index(2),index(3))), ...
  %   imag(SeqLoop.data.kSpaceOsZ(:,index(1),index(2),index(3)))])
  % xlabel('kSpace Point Signal ')
  % ylabel('Amplitude in pT')
  % %%
  % % end
end

% SeqLoop.AQSlice(1).plotImageHandle=2000;
% SeqLoop.data.data=SeqLoop.data.data/SeqLoop.Loops;
% [SeqLoop.data]=get_kSpaceAndImage(SeqLoop.data,SeqLoop.AQSlice(1));
% [SeqLoop.data]=plot_kSpaceAndImage(SeqLoop.data,SeqLoop.AQSlice(1));

if SeqOut.CorrectB0Read.Get || SeqOut.CorrectB0Read.Use
  SeqLoop.dataB0 = dataB0;
end

end
