function [Spectrum] = get_Spectrum(HW, SeqOut)
%% Get frequency spectrum from measured FID or echo (Spectroscopy)
%
%     [Spectrum] = get_Spectrum(HW, SeqOut)
%
% This function tries to find possible frequency drift and frequency offset
% properties of a series of averaging steps. It applies "best matching"
% corrections and averages the measurement signal after the corrections are
% applied.
%
%
% INPUT:
%
%   HW
%       HW object or structure
%
%   SeqOut
%       Structure with measurement settings and data as returned by, e.g.,
%       sequence_Spectrum. This structure also holds evaluation settings. Among
%       others the following fields can be set:
%
%     loopdata
%         Structure with the measurement data (of potentially several averaging
%         steps) as returned, e.g., by sequence_Spectrum. The field in this
%         structure have the following dimensions:
%         "samples" x "AQ window #" x "tRep" x "averages" x "X-nucleus"
%         This structure must contain the following fields:
%
%       time_all
%           Measurement time for each sample in s.
%       data
%           Complex measurement data at each sample in a.u.
%
%     Spectro
%         Structure with settings for the spectroscopic evalutation. This
%         structure can contain among others the following optional fields:
%
%       plot
%           Boolean value to indicate whether to plot the (intermediate)
%           results. (Default: false)
%
%       FrequencyDrift
%           Span for tested frequency drifts in Hz/s. (Default: 200)
%
%       nFrequencyDrift
%           Number of tested frequency drifts. The lower the number the worse
%           might be the match after the correction. But the higher the number,
%           the more memory is needed. (Default: 200)
%
%       FrequencyMatch
%           Span for tested frequency offsets in ppm. (Default: 10)
%
%       nFrequencyMatch
%           Number of tested frequency offsets. The lower the number the worse
%           might be the match after the correction. But the higher the number,
%           the more memory is needed. (Default: 200)
%
%       ppmStart
%           Start of the resulting spectrum in ppm from the center frequency.
%           (Default: 8)
%
%       ppmStop
%           End of the resulting spectrum in ppm from the center frequency.
%           (Default: -8)
%
%       ZeroFillFactor
%           Factor to the number of acquired samples in the FID which is filled
%           with zeros. This leads to a "smoother" spectrum with "tighter" bins.
%           (Default: 4)
%
%       ReferenceFID
%           Optionally, a reference FID with a sample with one single peak in
%           its spectrum and ideally long T2 can be acquired before the
%           measurement. This reference FID corresponds to an intrinsic line
%           shape of the magnet. Using the reference FID, a de-convolution can
%           be used to reconstruct a spectrum without this intrinsic line
%           widening. The reference FID is stored in the HW object and should be
%           acquired with the function "Find_ReferenceFID". (Default: 0)
%
%       FIDIdealT2Star
%           If ReferenceFID is set to true, this value is the T2* of the "ideal
%           FID" after the deconvolution of a spectrum with one single line. It
%           correlates to the line width of each line in the spectrum after the
%           deconvolution. It is related to the line width in ppm with the
%           following equation:
%             (1e6/(pi*SeqOut.AQ.Frequency(1))) / IdealLineWidthPPM
%           (Default: (1e6/(pi*SeqOut.AQ.Frequency(1))) / 0.4)
%
%       dualNuclear
%           Boolean value to indicate whether the input is a dual-nuclear
%           measurement. In this case, the data of both nuclei must be supplied
%           as described above. The data from the primary  nucleus is used for
%           frequency matching (for averaging). The correction data is applied
%           to both nuclei (taking the relation of HW.GammaDef to HW.GammaX into
%           account).
%           (Default: false)
%
%
% OUTPUT:
%
%   Spectrum
%       A structure with the results of the evaluation. It contains the
%       following fields:
%
%     FID
%         A column vector with the complex-valued time domain data of the
%         averaged data after frequency matching (drift and offset) in a.u..
%         In case of a dual-frequency measurement a matrix with two colums is
%         returned where the second column corresponds to the X nucleus signal.
%
%     Time
%         A column vector with the time corresponding to the FID vector (or
%         matrix) in seconds.
%
%     PhaseOffset
%         Phase at the start of the FID(s) in radians.
%
%     Spectrum
%         A column vector with the complex-valued frequency domain data of the
%         averaged data after frequency matching (drift and offset) in a.u..
%         In case of a dual-frequency measurement, a matrix with two colums is
%         returned where the second column corresponds to the X nucleus signal.
%
%     ppm
%         A column vector with the ppm scale corresponding to the Spectrum
%         vector (or matrix).
%
%     Frequency
%         A column vector with the frequencies corresponding to the FID vector.
%         In case of a dual-frequency measurement, a matrix with two colums is
%         returned where the second column corresponds to the X nucleus
%         spectrum.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2025 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

%% default input
if isemptyfield(SeqOut, 'Spectro'), SeqOut.Spectro = struct(); end

if isemptyfield(SeqOut.Spectro, 'plot'),  SeqOut.Spectro.plot = false;  end
if isemptyfield(SeqOut.Spectro, 'AQChannel'),  SeqOut.Spectro.AQChannel = 1;  end
if isemptyfield(SeqOut.Spectro, 'FrequencyDrift'), SeqOut.Spectro.FrequencyDrift = 200; end
if isemptyfield(SeqOut.Spectro, 'nFrequencyDrift'), SeqOut.Spectro.nFrequencyDrift = 200; end
if isemptyfield(SeqOut.Spectro, 'FrequencyMatch'), SeqOut.Spectro.FrequencyMatch = 10; end
if isemptyfield(SeqOut.Spectro, 'nFrequencyMatch'), SeqOut.Spectro.nFrequencyMatch = 200; end
if isemptyfield(SeqOut.Spectro, 'ppmStart'), SeqOut.Spectro.ppmStart = 8; end
if isemptyfield(SeqOut.Spectro, 'ppmStop'), SeqOut.Spectro.ppmStop = -8; end
if isemptyfield(SeqOut.Spectro, 'ZeroFillFactor'), SeqOut.Spectro.ZeroFillFactor = 4; end
if isemptyfield(SeqOut.Spectro, 'ReferenceFID'), SeqOut.Spectro.ReferenceFID = 0; end
if isemptyfield(SeqOut.Spectro, 'FIDIdealT2Star')
  SeqOut.Spectro.FIDIdealT2Star = (1e6/(pi*SeqOut.AQ.Frequency(1))) / 0.4;
end
if isemptyfield(SeqOut.Spectro, 'dualNuclear')
  if ~isemptyfield(SeqOut, 'dualNuclear')
    SeqOut.Spectro.dualNuclear = SeqOut.dualNuclear;
  else
    SeqOut.Spectro.dualNuclear = false;
  end
end

if isfield(SeqOut.Spectro, 'useSimulate') && SeqOut.Spectro.useSimulate

  %% simulated data
  if isemptyfield(SeqOut, 'Loops')
    if isemptyfield(SeqOut, 'loops')
      SeqOut.Loops = 4;
    else
      % fall back to legacy spelling of the field
      SeqOut.Loops = SeqOut.loops;
    end
  end
  if isemptyfield(SeqOut, 'RepetitionTime'),  SeqOut.RepetitionTime = 10;  end
  if isemptyfield(SeqOut, 'tFID'),  SeqOut.tFID = 3;  end
  if isemptyfield(SeqOut, 'AQ'),  SeqOut.AQ = [];  end
  if isemptyfield(SeqOut.AQ(1), 'fSample'),  [SeqOut.AQ(:).fSample] = deal(20e3);  end
  if isemptyfield(SeqOut.AQ(1), 'nSamples')
    for iAQ = 1:numel(SeqOut.AQ)
      SeqOut.AQ(iAQ).nSamples =  SeqOut.AQ(iAQ).fSample*SeqOut.tFID;
    end
  end
  if isemptyfield(SeqOut.AQ(1), 'Frequency'),  [SeqOut.AQ(:).Frequency] = deal(HW.fLarmor);  end
  if isemptyfield(SeqOut, 'Spectro'),  SeqOut.Spectro = struct();  end
  if isemptyfield(SeqOut.Spectro, 'ReferenceFID'),  SeqOut.Spectro.ReferenceFID = 0;  end
  if isemptyfield(SeqOut.Spectro, 'Loops')
    if isemptyfield(SeqOut.Spectro, 'loops')
      SeqOut.Spectro.Loops = SeqOut.Loops;
    else
      % fall back to legacy spelling of the field
      SeqOut.Spectro.Loops = SeqOut.Spectro.loops;
    end
  end
  if isemptyfield(SeqOut.Spectro.Simulate, 'B0LocalT2'),  SeqOut.Spectro.Simulate.B0LocalT2 = Inf;  end
  SeqOut.Loops = SeqOut.Spectro.Loops;
  tr = 1;
  nSamples = round(SeqOut.tFID*SeqOut.AQ((SeqOut.Spectro.AQChannel)).fSample);
  myTime = (0:nSamples-1).'/SeqOut.AQ((SeqOut.Spectro.AQChannel)).fSample;
  ppm2frequency = 1/SeqOut.Spectro.Simulate.frequency2ppm;
  nSamplesGrid = SeqOut.Spectro.Simulate.nSamplesGrid;
  FIDs = zeros(nSamples, SeqOut.Loops);
  for t = 1:size(FIDs, 2)
    FIDs(:,t)=sum(bsxfun(@times,SeqOut.Spectro.Simulate.Amp.*(randn(1)*SeqOut.Spectro.Simulate.AmpSTD+1),... % amplitude and amplitude nois
                                exp(-myTime*(1./SeqOut.Spectro.Simulate.T2)...                                % T2 decay
                                    -1i.*(SeqOut.Spectro.Simulate.PhaseOffsetRad+SeqOut.Spectro.Simulate.PhaseOffsetRadSTD)... % phase offset of FID
                                    -1i.*myTime.*2.*pi*(SeqOut.Spectro.Simulate.ppm*ppm2frequency...             % chemical shift
                                                        +SeqOut.Spectro.Simulate.PpmPerSec*ppm2frequency*SeqOut.RepetitionTime... %  frequency drift
                                                        +randn(1)*SeqOut.Spectro.Simulate.PpmPerSecSTD*ppm2frequency*SeqOut.RepetitionTime)... %  frequency drift STD
                                    -1i.*myTime.^2*(2*pi*(ones(size(SeqOut.Spectro.Simulate.ppm)).*SeqOut.Spectro.Simulate.PpmPerSec*ppm2frequency+randn(1)*SeqOut.Spectro.Simulate.PpmPerSecSTD*ppm2frequency)./2)))...
                                .*exp(-1i.*2^0.5/SeqOut.Spectro.Simulate.B0LocalT2.^0.5.*diff(myTime(1:2)).^0.5.*cumsum(randn(nSamples,numel(SeqOut.Spectro.Simulate.ppm)),1))...
                  , 2);  % frequency drift during FID
  % for t=1:size(FIDs,2)
  %   FIDs(:,t)=sum(bsxfun(@times,(Amp+randn(1)/100),exp(-myTime*(1./T2)-1i*myTime*2*pi*(fOffset+randn(1)*2)-1i*myTime.^2*(2*pi*(HZperSec+randn(1)*.5)./2))),2);
  %   % FIDs(:,t)=sum(bsxfun(@times,(Amp+randn(1)/100),exp(-myTime*(1./T2)-1i*myTime*2*pi*(fOffset+randn(1)*2)-1i*cumsum(myTime*2*pi*((HZperSec+randn(1)*0).*diff(myTime(1:2))),1))),2);
  % end
  end

  nAQX = 1;  % number of X frequencies
else
  tr = SeqOut.nEchos + 1;
  nAQX = size(SeqOut.loopdata.data, 5);  % number of X frequencies
  myTime = SeqOut.loopdata.time_all(1:SeqOut.AQ(SeqOut.Spectro.AQChannel).nSamples(tr),1,tr,1,1);
  FIDs = reshape(SeqOut.loopdata.data(1:SeqOut.AQ(SeqOut.Spectro.AQChannel).nSamples(tr),1,tr,:), ...
    SeqOut.AQ(SeqOut.Spectro.AQChannel).nSamples(tr), SeqOut.Loops, nAQX);
  nSamplesGrid = round(SeqOut.loopdata.time_all(1,1,tr,1, 1)*SeqOut.AQ(SeqOut.Spectro.AQChannel).fSample(1,tr));
  FIDs = [repmat(FIDs(1,:,:),nSamplesGrid,1,1); FIDs];
  myTime = [repmat(((-nSamplesGrid:-1)+round(myTime(1)*SeqOut.AQ(SeqOut.Spectro.AQChannel).fSample(1,tr))).',1,size(myTime,2))/SeqOut.AQ(SeqOut.Spectro.AQChannel).fSample(1,tr); myTime];
  nSamples = SeqOut.AQ(SeqOut.Spectro.AQChannel).nSamples(tr) + nSamplesGrid;
end
if nAQX > 1
  center_freq = [SeqOut.AQ(SeqOut.Spectro.AQChannel).Frequency(tr), SeqOut.AQ(SeqOut.Spectro.AQChannel).FrequencyX(tr)];
else
  center_freq = SeqOut.AQ(SeqOut.Spectro.AQChannel).Frequency(tr);
end
Frequency = bsxfun(@plus, ...
  get_FFTGrid(SeqOut.AQ(SeqOut.Spectro.AQChannel).fSample(tr), ...
  nSamples*SeqOut.Spectro.ZeroFillFactor), ...
  center_freq);
ppm = -bsxfun(@rdivide, ...
  get_FFTGrid(SeqOut.AQ(SeqOut.Spectro.AQChannel).fSample(tr), ...
              nSamples*SeqOut.Spectro.ZeroFillFactor), ...
  center_freq) * 1e6;
ppmStartI = find(ppm(:,1)<=SeqOut.Spectro.ppmStart, 1, 'first');
ppmStopI = find(ppm(:,1)>=SeqOut.Spectro.ppmStop, 1, 'last');

if SeqOut.Spectro.ReferenceFID
  % linear phase correction of Ref FID and Smoothing and generate Reference FID
  FIDRefRaw = zeros(size(myTime));

  % interpolate to grid of measurement
  FIDRefRaw((myTime>=HW.Spectro.ReferenceTime(1))) = ...
    interp1(HW.Spectro.ReferenceTime, HW.Spectro.ReferenceFid, myTime(myTime>=HW.Spectro.ReferenceTime(1)));
  % extrapolate towards zero on the same grid (if necessary)
  nSamplesRefGrid = find(myTime<=HW.Spectro.ReferenceTime(1), 1, 'last');
  if nSamplesRefGrid > 0
    % number of samples from reference FID used for extrapolation
%     nSamplesExtrap = min(numel(myTime), max(floor(numel(myTime)/50), nSamplesRefGrid*2));
    nSamplesExtrap = min([numel(myTime)-nSamplesRefGrid,... % limit
                          nSamplesRefGrid*2,... % normal
                          find(abs(FIDRefRaw(nSamplesRefGrid))/2>abs(FIDRefRaw(nSamplesRefGrid+(1:nSamplesRefGrid*2),:)),1,'first')]); % bad T2*
    % amplitude - extrapolate quadratically
%     [p, S, mu] = ...
%       polyfit(myTime(nSamplesRefGrid+(1:nSamplesExtrap)), ...
%               abs(FIDRefRaw(nSamplesRefGrid+(1:nSamplesExtrap))), 2);
%     [p,S,mu]=polyfit([-flipud(myTime(nSamplesRefGrid+(1:nSamplesRefGrid*2)));myTime(nSamplesRefGrid+(1:nSamplesRefGrid*2))],...
%                      [ flipud(abs(FIDRefRaw(nSamplesRefGrid+(1:nSamplesRefGrid*2))));abs(FIDRefRaw(nSamplesRefGrid+(1:nSamplesRefGrid*2)))],4); % extrapolate to horizontal
    [p,S,mu]=polyfit([-flipud(myTime(nSamplesRefGrid+(1:nSamplesExtrap)));myTime(nSamplesRefGrid+(1:nSamplesExtrap))],...
                     [ flipud(abs(FIDRefRaw(nSamplesRefGrid+(1:nSamplesExtrap))));abs(FIDRefRaw(nSamplesRefGrid+(1:nSamplesExtrap)))],4); % extrapolate to horizontal
    a = polyval(p, myTime(1:nSamplesRefGrid), S, mu);
%     AmpRef0=a(1);

    % phase - extrapolate linear
    [p, S, mu] = ...
      polyfit(myTime(nSamplesRefGrid+(1:nSamplesExtrap)), ...
      unwrap(angle(FIDRefRaw(nSamplesRefGrid+(1:nSamplesExtrap))),[],1), 1);
    PhaseRef0=polyval(p, myTime(1), S, mu);
    FIDRefRaw(1:nSamplesRefGrid) = a .* exp(1i*polyval(p, myTime(1:nSamplesRefGrid), S, mu));
    if any(abs(FIDRefRaw(1,:))/3>abs(FIDRefRaw(nSamplesRefGrid*3,:)))
      warning('PD:get_Spectrum:extrapolationFIDRefRawTooFar','extrapolation at the beginning of the FIDRefRaw too far')
    end
  else
    PhaseRef0=0;
  end
end



% extrapolate at the beginning of the FID (if necessary)
FIDsMean=mean(FIDs(nSamplesGrid+(1:nSamplesGrid*2),:),2);
StopI = min([nSamplesGrid*2,find(abs(FIDsMean(1))/2>abs(FIDsMean),1,'first')]);
% Amp0Corr=1.000;
Phase0Corr=0.00;
W=5*StopI*1;
if StopI > 1
  for tt = 1:(size(FIDs,2)*size(FIDs,3))
    % fitexpStartGrid=fit_exp( FIDs(nSamplesGrid+(1:StopI),t).*exp(-1i.*get_MeanPhaseDiffWeighted(FIDs(nSamplesGrid+(1:StopI),t)).*(nSamplesGrid+(1:StopI)).'),myTime(nSamplesGrid+(1:StopI),t),0,0,1,0,0);
    % FIDs(1:nSamplesGrid,t)=fitexpStartGrid.xminSingle(2).*exp(-myTime(1:nSamplesGrid)./fitexpStartGrid.xminSingle(3));
%     [p, S, mu] = polyfit(myTime(nSamplesGrid+(1:StopI)), abs(FIDs(nSamplesGrid+(1:StopI),t)), 2);  % extrapolate
    [p,S,mu]=polyfit([-flipud(myTime(nSamplesGrid+(1:StopI)));myTime(nSamplesGrid+(1:StopI))],...
                     [ flipud(abs(FIDs(nSamplesGrid+(1:StopI),tt)));abs(FIDs(nSamplesGrid+(1:StopI),tt))],4); % extrapolate to horizontal
    a = polyval(p, myTime(1:nSamplesGrid), S, mu);
    if PhaseRef0+Phase0Corr ~= 0
      FIDs(nSamplesGrid+1:end,tt) = FIDs(nSamplesGrid+1:end,tt).* exp(1i*(-PhaseRef0+Phase0Corr));
      [p, S, mu] = polyfit([flipud(-myTime(nSamplesGrid+(1:StopI)));                          repmat(myTime(1),W,1);  myTime(nSamplesGrid+(1:StopI))],...
        [flipud(unwrap(angle(conj(FIDs(nSamplesGrid+(1:StopI),tt))),[],1));zeros(W,1); unwrap(angle(FIDs(nSamplesGrid+(1:StopI),tt)),[],1)], 5);
      FIDs(1:nSamplesGrid,tt) = a .* exp(1i*polyval(p, myTime(1:nSamplesGrid), S, mu));
      FIDs(:,tt) = FIDs(:,tt).* exp(1i*(+PhaseRef0));
    else
      [p, S, mu] = polyfit(myTime(nSamplesGrid+(1:StopI)), unwrap(angle(FIDs(nSamplesGrid+(1:StopI),tt)),[],1), 1);
      FIDs(1:nSamplesGrid,tt) = a .* exp(1i*polyval(p, myTime(1:nSamplesGrid), S, mu));
    end
    % figure(251)
    % subplot(2,1,1)
    % plot(myTime(1:nSamplesGrid*4),(abs(FIDs(1:nSamplesGrid*4,t))))
    % subplot(2,1,2)
    % plot(myTime(1:nSamplesGrid*4),(angle(FIDs(1:nSamplesGrid*4,t))))
    % pause

  end
  if any(abs(FIDs(1,:))/3>abs(FIDs(nSamplesGrid+StopI,:)))
    warning('PD:get_Spectrum:extrapolationFIDsTooFar','extrapolation at the beginning of the FID too far')
  end
end


if SeqOut.Spectro.ReferenceFID
  % correct for potential frequency offset (before smoothing)
  nFilter = round(numel(FIDRefRaw)/500);
  PhaseCorr1 = get_MeanPhaseDiffWeighted(...
    conv(FIDRefRaw(round(numel(FIDRefRaw)./20):round(numel(FIDRefRaw)./5)), ...
         RaisedCosine(nFilter)./sum(RaisedCosine(nFilter)), 'valid'));
  FIDRef = FIDRefRaw .* exp(-1i * PhaseCorr1 .* (1:size(FIDRefRaw(:),1)).');

  % smooth reference FID (iteratively)
  RefFIDIn = FIDRef;
  nFilter = cumsum([3,30,300,300,300]);
  for n = nFilter
    FIDRef = conv(RefFIDIn, RaisedCosine(n)./sum(RaisedCosine(n)), 'same');
    FIDRef(1:n) = RefFIDIn(1:n);
    FIDRef(n+(1:n)) = (RefFIDIn(n+(1:n)).*(n:-1:1).' + FIDRef(n+(1:n)).*(1:n).') ./ (n+1);
    [p, S, mu] = polyfit(myTime(end-n*2:end), RefFIDIn((end-n*2:end)),4);
    FIDRef(end-n+1:end) = polyval(p, myTime(end-n+1:end), S, mu);
    RefFIDIn = FIDRef;
  end
  % correct for potential frequency offset (after smoothing)
  PhaseCorr2 = get_MeanPhaseDiffWeighted(FIDRef(nFilter(end):end));
  FIDRef = FIDRef .* exp(-1i * PhaseCorr2 .* (1:size(FIDRef(:),1)).');


  % find a suitable T2* or choose one and generate ideal FID
  if 0
    startI = 800;
    % fitexp = fit_exp( FIDs(startI:end,1).*exp(-1i.*get_MeanPhaseDiffWeighted(FIDs(startI:end,1)).*(startI:size(FIDs,1)).'),myTime(startI:end,1),1,0,1,0,0);
    fitexp = fit_exp(abs(RefFID(startI:end,1)), myTime(startI:end,1), 1, 0, 0, 0, 0);
    SeqOut.Spectro.FIDIdealT2Star = fitexp.xminSingle(3);  %  T2*
  end

  if isemptyfield(SeqOut.Spectro, 'FIDWindowDuration')
    SeqOut.Spectro.FIDWindowDuration = myTime(end);
  end
  if isemptyfield(SeqOut.Spectro, 'FIDWindowResolutionEnhancedDuration')
    SeqOut.Spectro.FIDWindowResolutionEnhancedDuration = SeqOut.Spectro.FIDWindowDuration;
  end
  cutFidIndexREEnd = find(myTime<=SeqOut.Spectro.FIDWindowResolutionEnhancedDuration, 1, 'last');
  % FID window function
  if isequal(SeqOut.Spectro.FIDWindowResolutionEnhanced, @LorentzToGauss)
    if isemptyfield(SeqOut.Spectro, 'FIDWindowResolutionEnhancedLB')
      SeqOut.Spectro.FIDWindowResolutionEnhancedLB = -2/2*cutFidIndexREEnd/cutFidIndexREEnd;
    end
    if isemptyfield(SeqOut.Spectro, 'FIDWindowResolutionEnhancedWmax')
      SeqOut.Spectro.FIDWindowResolutionEnhancedWmax = 1/6;
    end
    Window_ResolutionEnhanced = ...
      [LorentzToGauss(cutFidIndexREEnd, ...
      SeqOut.Spectro.FIDWindowResolutionEnhancedLB, ...
      SeqOut.Spectro.FIDWindowResolutionEnhancedWmax, ...
      SeqOut.AQ(SeqOut.Spectro.AQChannel).fSample(tr), 1); ...
      zeros(size(FIDs,1)-cutFidIndexREEnd, 1)];
    % Window_ResolutionEnhanced= LorentzToGauss(size(FIDs,1), -2/5*size(FIDs,1)/nSamples(tr), 1/10, SeqOut.AQ(SeqOut.Spectro.AQChannel).fSample(tr), 1);
  else
    Window_ResolutionEnhanced = ...
      [SeqOut.Spectro.FIDWindowResolutionEnhanced(cutFidIndexREEnd); ...
      zeros(size(FIDs,1)-cutFidIndexREEnd,1)];
  end

  % mytime(find(abs(FIDRef)<max(abs(FIDRef)),1,'last'))
  % abs(FIDRef(find(abs(FIDRef)<max(abs(FIDRef)),1,'last')))
  FidIdealStartValue = abs(FIDRef(1)); % start value
  FidIdeal = FidIdealStartValue .* exp(-myTime./SeqOut.Spectro.FIDIdealT2Star) .* Window_ResolutionEnhanced;

  % generate Fid Correction window
  FidCorr = FidIdeal ./ FIDRef;
  if SeqOut.Spectro.plot
    hf1 = figure(1000);
    clf(hf1);
    set(hf1, 'Name', 'FidRef FidIdeal and correction');
    axf(1) = subplot(3,1,1, 'Parent', hf1);
    plot(axf(1), ...
      myTime, abs(FIDRefRaw)/HW.RX.AmplitudeUnitScale, ...
      myTime, abs(FIDRef)/HW.RX.AmplitudeUnitScale, ...
      myTime, abs(FidIdeal)/HW.RX.AmplitudeUnitScale)
    ylabel(axf(1), {'amplitude', ['in ' HW.RX.AmplitudeUnit]});
    legend(axf(1), 'FIDRefRaw', 'FidRef', 'FidIdeal');
    axf(2) = subplot(3,1,2, 'Parent', hf1);
    plot(axf(2), ...
      myTime, angle(FIDRefRaw .* exp(-1i * (PhaseCorr1+PhaseCorr2) .* (1:size(FIDRef(:),1)).')), ...
      myTime, angle(FIDRef), ...
      myTime, angle(FidIdeal), ...
      myTime, angle(FidCorr))
    xlabel(axf(2), 'time in s')
    ylabel(axf(2), {'phase', 'in rad'});
    legend(axf(2), 'FIDRefRaw', 'FidRef', 'FidIdeal', 'FidCorr');
    axf(3) = subplot(3,1,3, 'Parent', hf1);
    plot(axf(3), myTime, abs(FidCorr))
    legend(axf(3), 'FidCorr');
    ylabel(axf(3), {'correction', 'factor'});
    linkaxes(axf(1:3), 'x');
  end
else
  FidCorr = ones(size(FIDs,1), 1);
  FidIdealStartValue = mean(abs(FIDs(1,:))); % start value
  FidIdeal = FidIdealStartValue .* exp(-myTime./SeqOut.Spectro.FIDIdealT2Star);
end

FidIdealIntegral = sum(abs(FidIdeal));

cutFidIndexEnd = find(myTime<=SeqOut.Spectro.FIDWindowDuration, 1, 'last');
w = reshape([SeqOut.Spectro.FIDWindow(cutFidIndexEnd);
             zeros(size(FIDs,1)-cutFidIndexEnd,1)] .* FidCorr, ...
            size(FIDs,1), 1);

% figure(345); plot(Window_ResolutionEnhanced);

if nAQX > 1
  % adjust phase in reference FID to center frequency of X nucleus
  w(:,1,2) = w .* ...
    exp(1i * ...
    (SeqOut.AQ(SeqOut.Spectro.AQChannel).FrequencyX(tr)/SeqOut.AQ(SeqOut.Spectro.AQChannel).Frequency(tr)) * ...
    unwrap(angle(w)) .* myTime);
end

FIDsw = bsxfun(@times, FIDs, w);

if SeqOut.Spectro.plot
  for iAQX = 1:nAQX
    hf2 = figure(1000+iAQX);
    clf(hf2);
    set(hf2, 'Name', 'FID + Window Function');
    ax2(1) = subplot(4,1,1, 'Parent', hf2);
    plot(ax2(1), ...
      repmat(myTime,1,size(FIDs,2)), abs(FIDs(:,:,iAQX)/HW.RX.AmplitudeUnitScale), ...
      myTime, abs(w(:,1,iAQX)/max(abs(w(:,1,iAQX)))*max(max(abs(FIDs(:,:,iAQX)))))/HW.RX.AmplitudeUnitScale, 'k--');
    ylabel(ax2(1), {'Raw FIDs', 'amplitude', ['in ' HW.RX.AmplitudeUnit]});
    ax2(2) = subplot(4,1,2, 'Parent', hf2);
    plot(ax2(2), ...
      repmat(myTime,1,size(FIDs,2)), angle(FIDs(:,:,iAQX)), ...
      myTime, angle(w(:,1,iAQX)), 'k--');
    ylabel(ax2(2), {'Raw FIDs', 'phase', 'in rad'});
    ax2(3) = subplot(4,1,3, 'Parent', hf2);
    plot(ax2(3), ...
      repmat(myTime,1,size(FIDs,2)), abs(FIDsw(:,:,iAQX))/HW.RX.AmplitudeUnitScale, ...
      myTime, FidIdeal/HW.RX.AmplitudeUnitScale, 'r--');
    ylabel(ax2(3), {'window FIDs', 'amplitude', ['in ' HW.RX.AmplitudeUnit]});
    ax2(4) = subplot(4,1,4, 'Parent', hf2);
    plot(ax2(4), repmat(myTime,1,size(FIDs,2)), angle(FIDsw(:,:,iAQX)));
    ylabel(ax2(4), {'window FIDs', 'phase', 'in rad'});
    xlabel(ax2(4), 'time in s');
    set(ax2, 'XGrid', 'on', 'YGrid', 'on');
    linkaxes(ax2, 'x');
  end
end

%% check for "most probable" frequency drift
% This is done by applying a range of different frequency drifts to every
% acquisition window and selecting the "best" one

% apply a range of different frequency drifts to signal
if SeqOut.Spectro.nFrequencyDrift == 1
  Omega1 = 0;
else
  Omega1 = linspace(-SeqOut.Spectro.FrequencyDrift/2, ...
    SeqOut.Spectro.FrequencyDrift/2, SeqOut.Spectro.nFrequencyDrift);
end
FIDsw_corrected = zeros(size(FIDsw), 'like', FIDsw);
Hf1 = zeros(1, SeqOut.Loops, numel(Omega1));
Omega1MinI = zeros(1, SeqOut.Loops);

% estimate size of correlation array
memCorrArray = numel(myTime) * SeqOut.Loops * SeqOut.Spectro.nFrequencyDrift * 8 * 2;
[userview, systemview] = memory();
numAvgInGroup = max(1, floor(SeqOut.Loops ...
  / max(1, memCorrArray ...
           / min(userview.MaxPossibleArrayBytes, systemview.PhysicalMemory.Available) * 4)));
startAvg = 1:numAvgInGroup:SeqOut.Loops;
endAvg = [numAvgInGroup:numAvgInGroup:SeqOut.Loops, SeqOut.Loops];

% loop over groups to limit memory consumption
for iGroup = 1:numel(startAvg)
  % only use primary gamma signal for frequency drift detection
  FIDswf1 = bsxfun(@times, FIDsw(:,startAvg(iGroup):endAvg(iGroup),1), ...
    reshape(exp(-1i*2*pi*(myTime.^2)*Omega1), [numel(myTime), 1, SeqOut.Spectro.nFrequencyDrift]));
  % calculate Fourier transform taking all frequency drifts into account
  Sf1 = fftshift(fft(FIDswf1, [], 1), 1);
  clear FIDswf1
  Sf1 = Sf1(round(ppmStartI/SeqOut.Spectro.ZeroFillFactor):round(ppmStopI/SeqOut.Spectro.ZeroFillFactor),:,:);
  % normalize all power spectra
  Pf1 = abs(Sf1).^2 ./ repmat(sum(abs(Sf1(:,:,1)).^2,1), [size(Sf1,1), 1, size(Sf1,3)]);
  clear Sf1
  % select "thinnest"/"most concentrated" power spectrum from range of tested frequency drifts
  Hf1(1,startAvg(iGroup):endAvg(iGroup),:) = -sum(Pf1.*log(Pf1),1);
  clear Pf1
  [~, Omega1MinI(startAvg(iGroup):endAvg(iGroup))] = ...
    min(Hf1(1,startAvg(iGroup):endAvg(iGroup),:), [], 3);
  Omega1Min = Omega1(Omega1MinI(startAvg(iGroup):endAvg(iGroup)));
  if any([Omega1MinI(startAvg(iGroup):endAvg(iGroup))==1, ...
          Omega1MinI(startAvg(iGroup):endAvg(iGroup))==endAvg(iGroup)])
    warning('Frequency drift out of correction range, increase Seq.Spectro.FrequencyDrift');
  end

  % apply frequency drift correction (to all gamma signals)
  % Omega1Min = repmat(Omega1Min, 1, 1, nAQX);
  if nAQX > 1
    Omega1Min(:,:,nAQX) = Omega1Min(:,:,1) / HW.GammaDef * HW.GammaX;
  end
  FIDsw_corrected(:,startAvg(iGroup):endAvg(iGroup),:) = ...
    FIDsw(:,startAvg(iGroup):endAvg(iGroup),:) .* ...
    reshape(exp(-1i .* (myTime.^2) * reshape(Omega1Min, 1, [])./2), size(FIDsw, 1), [], nAQX);
end

Hf1 = permute(Hf1, [3, 2, 1]);
Spectra = fftshift(fft(FIDsw_corrected, size(FIDsw,1)*SeqOut.Spectro.ZeroFillFactor, 1), 1);
SpectraRAWw = fftshift(fft(FIDsw, size(FIDsw,1)*SeqOut.Spectro.ZeroFillFactor, 1), 1);
SpectraRAW = fftshift(fft(FIDs, size(FIDsw,1)*SeqOut.Spectro.ZeroFillFactor, 1), 1);
Ir = (ppmStartI:ppmStopI).';
if SeqOut.Spectro.plot
  for iAQX = 1:nAQX
    hf3 = figure(1010+iAQX);
    set(hf3, 'Name', 'Spectrum - Frequency Drift');
    ax3(1) = subplot(4,1,1, 'Parent', hf3);
    plot(ax3(1), ppm(Ir,iAQX), abs(SpectraRAW(Ir,:,iAQX))/FidIdealIntegral);
    ylabel(ax3(1), ['Spectrum ', HW.RX.AmplitudeName]);
    grid(ax3(1), 'on');

    ax3(2) = subplot(4,1,2, 'Parent', hf3);
    plot(ax3(2), ppm(Ir,iAQX), abs(SpectraRAWw(Ir,:,iAQX))/FidIdealIntegral);
    ylabel(ax3(2), 'Spectrum with Window Function');
    grid(ax3(2), 'on');

    ax3(3) = subplot(4,1,3, 'Parent', hf3);
    plot(ax3(3), ppm(Ir,iAQX), abs(Spectra(Ir,:,iAQX))/FidIdealIntegral);
    ylabel(ax3(3), 'Spectrum Frequency Drift Corr.');
    xlabel(ax3(3), 'ppm');
    xlim(ax3(3), [SeqOut.Spectro.ppmStop, SeqOut.Spectro.ppmStart]);
    grid(ax3(3), 'on');

    set(ax3(1:3), 'XDir', 'reverse');
    linkaxes(ax3(1:3), 'x');

    ax3(4) = subplot(4,1,4, 'Parent', hf3);
    plot(ax3(4), repmat(Omega1.',1,size(Hf1,2)), Hf1);
    hold(ax3(4), 'on');
    plot(ax3(4), [Omega1(Omega1MinI); NaN(size(Omega1(Omega1MinI)))], ...
                 [Hf1(Omega1MinI+(0:size(Hf1,1):(size(Hf1,1)*(size(Hf1,2)-1)))); NaN(size(Omega1(Omega1MinI)))], 'x');
    hold(ax3(4), 'off');
    ylabel(ax3(4), 'Frequency Drift');
    xlabel(ax3(4), 'Hz/s');
    grid(ax3(4), 'on');
  end
  drawnow();
end

%% find "best matching" frequency offset
% This is done by shifting the spectra of each acquisition (potentially
% zero-filled) and selecting the combination with the most overlap.

if any(-SeqOut.Spectro.FrequencyMatch/2 < min(ppm, [], 1)) ...
    || any(SeqOut.Spectro.FrequencyMatch/2 > max(ppm, [], 1))
  warning('PD:get_Spectrum:narrowBandwidth', ...
    'The bandwidth of the signal limits the frequency match to %.1f ppm.', ...
    max(max(abs(ppm), [], 1)));
end

% select part of spectrum that is used for "fit"
Ir = (ppmStartI:ppmStopI).';

% select "steps" in which the spectra are shifted for the check
Omega0 = unique(round(linspace(-SeqOut.Spectro.FrequencyMatch/2, SeqOut.Spectro.FrequencyMatch/2, SeqOut.Spectro.nFrequencyMatch)/diff(ppm(1:2))));
Is = bsxfun(@plus, Ir, Omega0);  % matrix with all indices for "convolution" check
% first "un-shifted" spectrum is the reference
SpectraR = Spectra(Ir,1,1);

D = zeros(1, SeqOut.Loops, numel(Omega0));
Omega0Min = zeros(1, SeqOut.Loops);
SpectraMin = zeros([numel(Ir), SeqOut.Loops, nAQX], 'like', FIDsw);

% estimate size of correlation array
memCorrArray = numel(Ir) * SeqOut.Loops * SeqOut.Spectro.nFrequencyDrift * 8 * 2;
[userview, systemview] = memory();
numAvgInGroup = max(1, floor(SeqOut.Loops ...
  / max(1, memCorrArray ...
           / min(userview.MaxPossibleArrayBytes, systemview.PhysicalMemory.Available) * 4)));
startAvg = 1:numAvgInGroup:SeqOut.Loops;
endAvg = [numAvgInGroup:numAvgInGroup:SeqOut.Loops, SeqOut.Loops];

% loop over groups to limit memory consumption
for iGroup = 1:numel(startAvg)
  numInGroup = endAvg(iGroup) - startAvg(iGroup) + 1;
  % normalize power spectra and repeat for "convolution" check
  Pr = repmat(abs(SpectraR).^2./repmat(sum(abs(SpectraR).^2,1), [size(SpectraR,1),1]), [1, numInGroup, size(Is,2)]);
  % matrix with all "shifted" spectra for "convolution" check
  Is2 = Is;
  % zero-fill if shift exceeds limits of spectrum
  Is2(Is<1 | Is>size(Spectra, 1)) = 1;
  SpectraS = Spectra(Is2(:),startAvg(iGroup):endAvg(iGroup),1);
  SpectraS(Is(:)<1 | Is(:)>size(Spectra, 1),:) = 0;
  SpectraS = permute(reshape(SpectraS, [size(Is,1),size(Is,2),numInGroup]), [1,3,2]);
  % SpectraS = permute(reshape(Spectra(Is2(:),:,1), [size(Is,1),size(Is,2),size(Spectra,2)]), [1,3,2]);
  % normalized power spectra
  Ps = abs(SpectraS).^2./repmat(sum(abs(SpectraS(:,:,1)).^2,1), [size(SpectraS,1),1,size(SpectraS,3)]);
  % calculate determinant with all selected shifts
  D(1,startAvg(iGroup):endAvg(iGroup),:) = (2-sum((Ps.*Pr).^0.5,1)).^0.5;  % convn!
  clear Ps Pr SpectraS

  % get minimum determinant -> best overlap
  [~, Omega0MinI] = min(D(1,startAvg(iGroup):endAvg(iGroup),:), [], 3);
  if any([Omega0MinI==1, Omega0MinI==numel(Omega0)])
    warning('Frequency match out of correction range, increase Seq.Spectro.FrequencyMatch');
  end
  Omega0Min(startAvg(iGroup):endAvg(iGroup)) = Omega0(Omega0MinI);
  Imin = ...
    bsxfun(@plus, Ir, ...
           Omega0Min(startAvg(iGroup):endAvg(iGroup)) ...
           + (0:size(Spectra,1):size(Spectra,1)*size(Omega0Min(startAvg(iGroup):endAvg(iGroup)),2)-1)) ...
    + size(Spectra, 1) * (startAvg(iGroup)-1);
  SpectraMin(:,startAvg(iGroup):endAvg(iGroup),1) = reshape(Spectra(Imin(:)), size(Imin));
end

SpectraOmega0 = zeros(size(Spectra));
SpectraOmega0(:,:,1) = fftshift(fft(FIDsw_corrected(:,:,1) .* ...
  exp(-2i*pi * bsxfun(@times, (1:size(FIDsw_corrected,1)).' / size(FIDsw_corrected,1), Omega0Min/SeqOut.Spectro.ZeroFillFactor)), ...
  size(FIDsw, 1) * SeqOut.Spectro.ZeroFillFactor, 1), 1);
SpectraMin(:,:,1) = SpectraOmega0(Ir,:,1);
if nAQX > 1
  SpectraOmega0(:,:,2) = fftshift(fft(FIDsw_corrected(:,:,2) .* ...
    exp(-2i*pi * bsxfun(@times, (1:size(FIDsw_corrected,1)).' / size(FIDsw_corrected,1), Omega0Min/SeqOut.Spectro.ZeroFillFactor) / HW.GammaDef * HW.GammaX), ...
    size(FIDsw, 1) * SeqOut.Spectro.ZeroFillFactor, 1), 1);
  SpectraMin(:,:,2) = SpectraOmega0(Ir,:,2);
end
FIDswC = ifft(ifftshift(SpectraOmega0,1), [], 1);


if SeqOut.Spectro.plot

  % figure(1041); plot(ppm(:,1), real(SpectraOmega0(:,:,1)));
  % figure(1042); plot(ppm(:,2), real(SpectraOmega0(:,:,2)));

  hf4 = figure(1021); clf(hf4);
  set(hf4, 'Name', 'Frequency Match');
  ax4 = axes('Parent', hf4);
  plot(ax4, Omega0*abs(diff(ppm(1:2,1))), squeeze(D));
  ylabel(ax4, 'frequency match');
  xlabel(ax4, 'ppm');
  grid(ax4, 'on');

  for iAQX = 1:nAQX
    hf5 = figure(1030+iAQX);
    set(hf5, 'Name', 'Spectrum - Frequency Matched');
    ax5(1) = subplot(4,1,1, 'Parent', hf5);
    fid_matched = reshape(mean(FIDswC(1:end/SeqOut.Spectro.ZeroFillFactor,:,:), 2), [], nAQX);
    plot(ax5(1), ...
      myTime, abs(fid_matched(:,iAQX)), ...
      myTime, real(fid_matched(:,iAQX)), ...
      myTime, imag(fid_matched(:,iAQX)));
    xlabel(ax5(1), 'time in s');
    ylabel(ax5(1), 'amplitude');
    legend(ax5(1), 'abs', 'real', 'imag');
    grid(ax5(1), 'on');

    ax5(2) = subplot(4,1,2, 'Parent', hf5);
    plot(ax5(2), repmat(ppm(Ir,iAQX), [1,size(SpectraMin,2)+1]), ...
                 [abs(SpectraMin(:,:,iAQX)),mean(abs(SpectraMin(:,:,iAQX)),2)]/FidIdealIntegral);
    ylabel(ax5(2), 'absolute spectrum');
    grid(ax5(2), 'on');

    ax5(3) = subplot(4,1,3, 'Parent', hf5);
    plot(ax5(3), repmat(ppm(Ir,iAQX), [1,size(SpectraMin,2)+1]), ...
                 [angle(SpectraMin(:,:,iAQX)),angle(mean(SpectraMin(:,:,iAQX),2))]/FidIdealIntegral);
    ylabel(ax5(3), 'spectrum phase');
    grid(ax5(3), 'on');

    ax5(4) = subplot(4,1,4, 'Parent', hf5);
    plot(ax5(4), ppm(Ir,iAQX), [mean(abs(SpectraMin(:,:,iAQX)),2), mean(real(SpectraMin(:,:,iAQX)),2), mean(imag(SpectraMin(:,:,iAQX)),2)]/FidIdealIntegral);
    ylabel(ax5(4), 'spectrum');
    xlabel(ax5(4), 'ppm');
    grid(ax5(4), 'on');

    set(ax5(2:4), 'XDir', 'reverse');
    linkaxes(ax5(2:4), 'x');
  end
end

PhaseOffset = angle(mean(FIDswC(1,:),2));
SpectraMinCP = bsxfun(@times, SpectraMin, exp(-1i.*PhaseOffset));

Spectrum.FID = reshape(mean(FIDswC(1:end/SeqOut.Spectro.ZeroFillFactor,:,:), 2), [], nAQX);
Spectrum.PhaseOffset = PhaseOffset;
Spectrum.Time = myTime;
Spectrum.Spectrum = reshape(mean(SpectraMinCP ,2), [], nAQX)/FidIdealIntegral;
Spectrum.FidIdealIntegral = FidIdealIntegral;
Spectrum.ppm = ppm(Ir,:);
Spectrum.Frequency = Frequency(Ir,:);

end
