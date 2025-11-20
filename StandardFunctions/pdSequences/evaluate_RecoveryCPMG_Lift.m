function [SeqOut, data] = evaluate_RecoveryCPMG_Lift(HW, SeqOut, data)
%% Evaluate recovery CPMG profiles with sample lift
%
%   evaluate_RecoveryCPMG_Lift(HW, SeqOut, data)
%
% INPUT:
%   HW
%         HW structure or object
%
%   SeqOut
%         sequence structure as returned by "sequence_RecoveryCPMG" with the
%         following additional (optional) fields:
%     doT2iLaplace1D
%           Boolean value that activates calculating and displaying the T2
%           spectrum and profile. (default: (SeqOut.nEcho > 3) )
%     T2iLaplace1D
%           structure with properties for the inverse Laplace transform of the
%           CPMG Echo trains (T2 spectrum) at SeqOut.FitT2AtTau1.
%           See "get_iLaplace1D" for available settings.
%
%   data
%         array of structures containing the measurement data as returned by
%         "sequence_RecoveryCPMG" for each lift position.
%
% OUTPUT:
%   data
%         Same as input structure "data" but with the following added fields
%     T2iLaplace1D
%           array of structures containing the results of the inverse Laplace
%           transform of the CPMG Echo trains at SeqOut.FitT2AtTau1
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% input validation
if nargin ~= 3
  error('PD:evaluate_RecoveryCPMG:nargin', 'Number of input arguments must be 3');
end
if ~isstruct(HW) && ~isa(HW, 'PD.HWClass')
  error('PD:evaluate_RecoveryCPMG:wrongHW', '"HW" must be a structure or PD.HWClass object');
end
if ~isstruct(SeqOut)
  error('PD:evaluate_RecoveryCPMG:wrongSeq', ...
    '"SeqOut" must be a structure as returned by "sequence_RecoveryCPMG".');
end
if ~isstruct(data)
  error('PD:evaluate_RecoveryCPMG:wrongSeq', ...
    '"data" must be a structure as returned by "sequence_RecoveryCPMG" or an array of these structures.');
end

%% default values for fields in SeqOut
if isemptyfield(SeqOut, 'doEvaluation')
  SeqOut.doEvaluation = true;
end
if isemptyfield(SeqOut, 'doDisplay')
  SeqOut.doDisplay = true;
end

if ~isfield(SeqOut, 'fullScaleReference'), SeqOut.fullScaleReference = struct(); end
if isemptyfield(SeqOut.fullScaleReference, 'Amplitude')
  SeqOut.fullScaleReference.Amplitude = 1;
end
if isemptyfield(SeqOut.fullScaleReference, 'CrossSection')
  SeqOut.fullScaleReference.CrossSection = 1;
end
if isemptyfield(SeqOut.fullScaleReference, 'Thickness')
  SeqOut.fullScaleReference.Thickness = 1;
end
if isemptyfield(SeqOut, 'sampleCrossSection')
  SeqOut.sampleCrossSection = 1;
end

if isemptyfield(SeqOut, 'hParentT1Profile')
  SeqOut.hParentT1Profile = 201;
end
if isemptyfield(SeqOut, 'hParentT2Profile')
  SeqOut.hParentT2Profile = 202;
end
if isemptyfield(SeqOut, 'hParentIntensityProfile')
  SeqOut.hParentIntensityProfile = 203;
end
if isemptyfield(SeqOut, 'hParentT2iLaplaceProfile')
  SeqOut.hParentT2iLaplaceProfile = 221;
end

if isemptyfield(SeqOut, 'doT2iLaplace1D'), SeqOut.doT2iLaplace1D = (SeqOut.nEcho > 3); end
if isemptyfield(SeqOut, 'T2iLaplace1D'), SeqOut.T2iLaplace1D = struct(); end

%% evaluation and plot results
iLiftLoop = numel(data);

% FIXME:  For the GUI to work correctly, the lift structure must not be
%         overwritten. What is the correct thing to do for the dummy lift
%         instead?
% if isstruct(HW.Lift) || HW.Lift.isDummy || isemptyfield(data(1), 'lift') % || ~SeqOut.Lift.useLift
%   clear lift
%   if isemptyfield(SeqOut.Lift, 'Displacement'), SeqOut.Lift.Displacement = 1e-3; end
%   pos = num2cell((1:iLiftLoop)*SeqOut.Lift.Displacement); %lift = [data(:).lift];
%   [lift(1:iLiftLoop).position] = deal(pos{:}); lift = num2cell(lift); [data(:).lift] = deal(lift{:});
% end

if SeqOut.doEvaluation
  clear dataL;
  if (SeqOut.nTau1 > 3) && (SeqOut.nEcho > 3)
    % FIXME: This part isn't ready yet for lift slices
    % SeqOut.iLaplace2D.T1Start=SeqOut.Tau1Start*3;
    % SeqOut.iLaplace2D.T1End=SeqOut.Tau1End/2;
    % SeqOut.iLaplace2D.T2Start=SeqOut.tEcho*1;
    % SeqOut.iLaplace2D.T2End=SeqOut.tEcho*SeqOut.nEcho/11;
    % SeqOut.iLaplace2D.IgnoreFirstEcho=5;
    % data.Tau1Time=data.Tau1Time*1;
    % data.EchoTime=data.EchoTime*1;
    % SeqOut.iLaplace2D.T1Start=[];
    % SeqOut.iLaplace2D.T1End=[];
    % SeqOut.iLaplace2D.T2Start=[];
    % SeqOut.iLaplace2D.T2End=[];
    % SeqOut.iLaplace2D.IgnoreFirstEcho=0;
    % SeqOut.iLaplace2D.n_grid=50;
    % SeqOut.iLaplace2D.RingFilter=1;
    % SeqOut.iLaplace2D.LastEchoTrainCorrection=1;
    % SeqOut.iLaplace2D.QualityFactor=10;
    % SeqOut.iLaplace2D.n_tau2Smooth   =  0;
    % SeqOut.iLaplace2D.n_Lcomp     =  101;
    % SeqOut.iLaplace2D.n_Tcomp     =  101;
    [dataL(1), SeqOutiL] = get_iLaplace2D(data(1), SeqOut);
  end

  % SeqOut.intensityAtEcho = 3;
  data(iLiftLoop).T2iLaplace1D = struct();
  if SeqOut.doT2iLaplace1D
    % inverse Laplace 1D
    dataiLaplace1D = data(iLiftLoop);
    SeqiL = SeqOut;
    clear dataiL;
    for t = SeqOut.FitT2AtTau1
      dataiLaplace1D.DataAmplitude = dataiLaplace1D.MeanEchoTau1PhaseCorrected(:,t);
      dataiLaplace1D.DataTime = dataiLaplace1D.EchoTime(:,t);
      if isemptyfield(SeqOut.T2iLaplace1D, 'SpectrumTimeStart')
        SeqiL.iLaplace1D.SpectrumTimeStart = dataiLaplace1D.DataTime(1)/2;%diff(data.DataTime(1:2))*2;
      else
        SeqiL.iLaplace1D.SpectrumTimeStart = SeqOut.T2iLaplace1D.SpectrumTimeStart;
      end
      if isemptyfield(SeqOut.T2iLaplace1D, 'SpectrumTimeEnd')
        SeqiL.iLaplace1D.SpectrumTimeEnd = dataiLaplace1D.DataTime(end)*2;%  max(data.DataTime)*10;
      else
        SeqiL.iLaplace1D.SpectrumTimeEnd = SeqOut.T2iLaplace1D.SpectrumTimeEnd;
      end
      if isemptyfield(SeqOut.T2iLaplace1D, 'SpectrumTimeStartCut')
        SeqiL.iLaplace1D.SpectrumTimeStartCut = SeqiL.iLaplace1D.SpectrumTimeStart*1.4;%diff(data.DataTime(1:2))*1;
      else
        SeqiL.iLaplace1D.SpectrumTimeStartCut = SeqOut.T2iLaplace1D.SpectrumTimeStartCut;
      end
      if isemptyfield(SeqOut.T2iLaplace1D, 'SpectrumTimeEndCut')
        SeqiL.iLaplace1D.SpectrumTimeEndCut = SeqiL.iLaplace1D.SpectrumTimeEnd/1.1;%diff(data.DataTime(1:2))*1;
      else
        SeqiL.iLaplace1D.SpectrumTimeEndCut = SeqOut.T2iLaplace1D.SpectrumTimeEndCut;
      end
      if isemptyfield(SeqOut.T2iLaplace1D, 'nSpectrum')
        SeqiL.iLaplace1D.nSpectrum = 1000;
      else
        SeqiL.iLaplace1D.nSpectrum = SeqOut.T2iLaplace1D.nSpectrum;
      end
      if isemptyfield(SeqOut.T2iLaplace1D, 'QualityFactor')
        SeqiL.iLaplace1D.QualityFactor = 100;
      else
        SeqiL.iLaplace1D.QualityFactor = SeqOut.T2iLaplace1D.QualityFactor;
      end
      dataiL(t) = get_iLaplace1D(dataiLaplace1D, SeqiL);
    end
    data(iLiftLoop).T2iLaplace1D = [dataiL(:).iLaplace1D];
  end

end

if SeqOut.doDisplay
  lift = [data(:).lift];
  if iLiftLoop >= 1
    if SeqOut.hParentT2Profile ~= 0 && SeqOut.nEcho > 3
      % T2 profile
      T2Lift = [data(:).T2];
      if ishghandle(SeqOut.hParentT2Profile, 'figure') || isa(SeqOut.hParentT2Profile, 'double')
        hf = figure(SeqOut.hParentT2Profile);
        set(hf, 'Name', 'T2 Profile');
        clf(hf);
        SeqOut.hParentT2Profile = hf;
      elseif ishghandle(SeqOut.hParentT2Profile, 'uipanel')
        % delete all current children
        delete(get(SeqOut.hParentT2Profile, 'Children'));
      else
        error('"SeqOut.hParentT2Profile" must be a handle to a figure or uipanel.');
      end
      hax = axes(SeqOut.hParentT2Profile);
      plot(hax, [T2Lift(:).tau], [lift(1:numel(T2Lift)).position],'-x');
      title(hax, 'T2 - Singular Exponential Fit');
      ylabel(hax, 'Sample Lift Position in m');
      xlabel(hax, 'T2 in s');
      xlim(hax, [0 SeqOut.T2EstimatedMax]);
      set(hax, 'YDir', 'reverse');
      grid(hax, 'on');
    end

    if SeqOut.nTau1 > 1
      % T1 profile
      T1Lift = [data(:).T1];
      if ishghandle(SeqOut.hParentT1Profile, 'figure') || isa(SeqOut.hParentT1Profile, 'double')
        hf = figure(SeqOut.hParentT1Profile);
        set(hf, 'Name', 'T1 Profile');
        clf(hf);
        SeqOut.hParentT1Profile = hf;
      elseif ishghandle(SeqOut.hParentT1Profile, 'uipanel')
        % delete all current children
        delete(get(SeqOut.hParentT1Profile, 'Children'));
      else
        error('"SeqOut.hParentT1Profile" must be a handle to a figure or uipanel.');
      end
      hax = axes(SeqOut.hParentT1Profile);
      plot(hax, [T1Lift(:).tau], [lift(1:numel(T1Lift)).position],'-x');
      title(hax, 'T1 - Singular Exponential Fit');
      ylabel(hax, 'Sample Lift Position in m');
      xlabel(hax, 'T1 in s');
      set(hax, 'YDir', 'reverse');
      xlim(hax, [0 SeqOut.T1EstimatedMax]);
      grid(hax, 'on');
    end

    if SeqOut.hParentIntensityProfile ~= 0
      % intensity profile of n-th Echo
      meanEchoTau1 = cat(4, data(:).MeanEchoTau1PhaseCorrected);
      if ishghandle(SeqOut.hParentIntensityProfile, 'figure') || isa(SeqOut.hParentIntensityProfile, 'double')
        hf = figure(SeqOut.hParentIntensityProfile);
        set(hf, 'Name', 'Intensity Profile');
        clf(hf);
        SeqOut.hParentIntensityProfile = hf;
      elseif ishghandle(SeqOut.hParentIntensityProfile, 'uipanel')
        % delete all current children
        delete(get(SeqOut.hParentIntensityProfile, 'Children'));
      else
        error('"SeqOut.hParentIntensityProfile" must be a handle to a figure or uipanel.');
      end
      hax = axes(SeqOut.hParentIntensityProfile);
      % FIXME: "real" ist eigentlich richtig, wenn die Phase korrekt
      % geschoben ist. Die Phase kann fuer die Leermessung nicht getrennt
      % korrigiert werden. Deshalb wird die eigentliche Messung (vorerst)
      % auch nicht phasenkorrigiert.
      % Solange das nicht "richtig" implementiert ist, hier erstmal "abs":
      intensityProfile = squeeze(abs(meanEchoTau1(SeqOut.intensityProfileEchoNum,SeqOut.intensityProfileTau1,1,:)));
      plot(hax, intensityProfile, [lift(1:numel(intensityProfile)).position], '-x');
      title(hax, sprintf('Mean Intensity at Echo #%d, Tau_1 #%d', ...
        SeqOut.intensityProfileEchoNum,SeqOut.intensityProfileTau1));
      ylabel(hax, 'Sample Lift Position in m');
      xlabel(hax, 'Amplitude');
      set(hax, 'YDir', 'reverse');
      grid(hax, 'on');
    end
  end

  if SeqOut.doT2iLaplace1D
    % plot profile of T2 spectra
    alliLaplace1D = cat(1, data.T2iLaplace1D);
    dataTime = cat(2, alliLaplace1D(:,1).DataTime);
    dataAmplitude = bsxfun(@times, cat(2, alliLaplace1D(:,1).DataAmplitude), ...
      cat(2, alliLaplace1D(:,1).FullScaleAmplitude));
    if ishghandle(SeqOut.hParentT2iLaplaceProfile, 'figure') || isa(SeqOut.hParentT2iLaplaceProfile, 'double')
      hFigure = figure(SeqOut.hParentT2iLaplaceProfile);
      set(hFigure, 'Name', 'Intensity Profile');
      clf(hFigure);
      SeqOut.hParentT2iLaplaceProfile = hFigure;
      hParent = hFigure;
    elseif ishghandle(SeqOut.hParentT2iLaplaceProfile, 'uipanel')
      % delete all current children
      delete(get(SeqOut.hParentT2iLaplaceProfile, 'Children'));
      hParent = SeqOut.hParentT2iLaplaceProfile;
      hFigure = ancestor(hParent, 'figure');
    else
      error('"SeqOut.hParentT2iLaplaceProfileplace" must be a handle to a figure or uipanel.');
    end
    haxEcho = subplot(3,1,1, 'Parent', SeqOut.hParentT2iLaplaceProfile, 'Tag', 'CPMG_Intensity');
    grid(haxEcho, 'on');
    liftPosition = [lift(1:size(dataTime, 2)).position];
    liftStep = SeqOut.Lift.Displacement; % FIXME: This only works correctly for uniform steps
    liftPosition = liftPosition - liftStep/2;
    liftPosition = [liftPosition, liftPosition(end) + liftStep];
    surface(haxEcho, [dataTime, dataTime(:,end)], ...
      repmat(liftPosition, size(dataTime, 1), 1), ...
      [dataAmplitude, dataAmplitude(:,end)], ...
      'MeshStyle', 'none', ...
      'FaceColor', 'flat');
    set(haxEcho, 'XScale', 'log', 'YDir', 'reverse');
    view(haxEcho, 2);
    title(haxEcho, 'Echo Amplitude in parts water');
    xlabel(haxEcho, 'Echo Ttime in s');
    ylabel(haxEcho, 'Lift Position in m', 'HorizontalAlignment', 'center');
    zlabel(haxEcho, 'Amplitude');

    spectrumTime = cat(2, alliLaplace1D(:,1).SpectrumTime);
    spectrumAmplitude = bsxfun(@times, cat(2, alliLaplace1D(:,1).SpectrumAmplitude), ...
      cat(2, alliLaplace1D(:,1).FullScaleAmplitude));
    haxSpectrum = subplot(3,1,2, 'Parent', SeqOut.hParentT2iLaplaceProfile, 'Tag', 'CPMG_Spectrum');
    surface(haxSpectrum, [spectrumTime, spectrumTime(:,end)], ...
      repmat(liftPosition, size(spectrumTime, 1), 1), ...
      [spectrumAmplitude, spectrumAmplitude(:,end)], ...
      'MeshStyle', 'none', ...
      'FaceColor', 'flat');
    grid(haxSpectrum, 'on');
    set(haxSpectrum, 'XScale', 'log', 'YDir', 'reverse');
    view(haxSpectrum, 2);
    title(haxSpectrum, 'T2 Spectrum');
    xlabel(haxSpectrum, 'T2 in s');
    ylabel(haxSpectrum, 'Lift Position in m', 'HorizontalAlignment', 'center');
    zlabel(haxSpectrum, 'Amplitude  in parts water');

    cumSpectrumAmplitude = cumsum(spectrumAmplitude, 1, 'reverse');
    haxCumSpectrum = subplot(3,1,3, 'Parent', SeqOut.hParentT2iLaplaceProfile, 'Tag', 'CPMG_CumulativeSpectrum');
    surface(haxCumSpectrum, [spectrumTime, spectrumTime(:,end)], ...
      repmat(liftPosition, size(spectrumTime, 1), 1), ...
      [cumSpectrumAmplitude, cumSpectrumAmplitude(:,end)], ...
      'MeshStyle', 'none', ...
      'FaceColor', 'flat');
    grid(haxCumSpectrum, 'on');
    set(haxCumSpectrum, 'XScale', 'log', 'YDir', 'reverse');
    view(haxCumSpectrum, 2);
    title(haxCumSpectrum, 'Cumulative T2 Spectrum');
    xlabel(haxCumSpectrum, 'T2 in s');
    ylabel(haxCumSpectrum, 'Lift Position in m', 'HorizontalAlignment', 'center');
    zlabel(haxCumSpectrum, 'Amplitude in parts water');

    % link rotation
    lp = linkprop([haxEcho haxSpectrum haxCumSpectrum], {'View'});
    setappdata(hParent, 'profile2dLinkView', lp);
    rotate3d(hFigure, 'on');
    linkaxes([haxSpectrum,haxEcho,haxCumSpectrum],'x')
  end
end

end
