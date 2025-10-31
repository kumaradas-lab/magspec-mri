function [Spectrum] = get_Spectrum(HW, SeqOut)
%% Get frequency spectrum from measured FID or echo (Spectroscopy)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2020 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

if isemptyfield(SeqOut.Spectro, 'plot'),  SeqOut.Spectro.plot = 0;  end
if isemptyfield(SeqOut.Spectro, 'AQChannel'),  SeqOut.Spectro.AQChannel = 1;  end

if isfield(SeqOut.Spectro, 'Simulate') && SeqOut.Spectro.useSimulate

  %% simulated data
  if isemptyfield(SeqOut, 'loops'),  SeqOut.loops = 4;  end
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
  if isemptyfield(SeqOut.Spectro, 'loops'),  SeqOut.Spectro.loops = SeqOut.loops;  end
  if isemptyfield(SeqOut.Spectro.Simulate, 'B0LocalT2'),  SeqOut.Spectro.Simulate.B0LocalT2 = Inf;  end
  SeqOut.loops = SeqOut.Spectro.loops;
  tr = 1;
  nSamples = round(SeqOut.tFID*SeqOut.AQ((SeqOut.Spectro.AQChannel)).fSample);
  myTime = (0:nSamples-1).'/SeqOut.AQ((SeqOut.Spectro.AQChannel)).fSample;
  ppm2frequency = 1/SeqOut.Spectro.Simulate.frequency2ppm;
  nSamplesGrid = SeqOut.Spectro.Simulate.nSamplesGrid;
  FIDs = zeros(nSamples, SeqOut.loops);
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

else
  tr = SeqOut.nEchos + 1;
  myTime = SeqOut.loopdata.time_all(1:SeqOut.AQ(SeqOut.Spectro.AQChannel).nSamples(tr),1,tr,1);
  FIDs = squeeze(SeqOut.loopdata.data(1:SeqOut.AQ(SeqOut.Spectro.AQChannel).nSamples(tr),1,tr,:));
  nSamplesGrid = round(SeqOut.loopdata.time_all(1,1,tr,1)*SeqOut.AQ(SeqOut.Spectro.AQChannel).fSample(1,tr));
  FIDs = [repmat(FIDs(1,:),nSamplesGrid,1); FIDs];
  myTime = [repmat(((-nSamplesGrid:-1)+round(myTime(1)*SeqOut.AQ(SeqOut.Spectro.AQChannel).fSample(1,tr))).',1,size(myTime,2))/SeqOut.AQ(SeqOut.Spectro.AQChannel).fSample(1,tr); myTime];
  nSamples = SeqOut.AQ(SeqOut.Spectro.AQChannel).nSamples(tr) + nSamplesGrid;
end
Frequency = get_FFTGrid(SeqOut.AQ(SeqOut.Spectro.AQChannel).fSample(tr),nSamples*SeqOut.Spectro.ZeroFillFactor)+SeqOut.AQ(SeqOut.Spectro.AQChannel).Frequency(tr);
ppm = -(get_FFTGrid(SeqOut.AQ(SeqOut.Spectro.AQChannel).fSample(tr),nSamples*SeqOut.Spectro.ZeroFillFactor)./SeqOut.AQ(SeqOut.Spectro.AQChannel).Frequency(tr))*1e6;
ppmStartI = find(ppm<SeqOut.Spectro.ppmStart, 1, 'first');
ppmStopI = find(ppm>SeqOut.Spectro.ppmStop, 1, 'last');


% extrapolate at the beginning of the FID
StopI = nSamplesGrid*2;
for t = 1:size(FIDs,2)
  % fitexpStartGrid=fit_exp( FIDs(nSamplesGrid+(1:StopI),t).*exp(-1i.*get_MeanPhaseDiffWeighted(FIDs(nSamplesGrid+(1:StopI),t)).*(nSamplesGrid+(1:StopI)).'),myTime(nSamplesGrid+(1:StopI),t),0,0,1,0,0);
  % FIDs(1:nSamplesGrid,t)=fitexpStartGrid.xminSingle(2).*exp(-myTime(1:nSamplesGrid)./fitexpStartGrid.xminSingle(3));
  [p, S, mu] = polyfit(myTime(nSamplesGrid+(1:StopI)), abs(FIDs(nSamplesGrid+(1:StopI),t)), 2);  % extrapolate
  % [p,S,mu]=polyfit([-flipud(myTime(nSamplesGrid+(1:StopI)));myTime(nSamplesGrid+(1:StopI))],...
  %                  [ flipud(abs(FIDs(nSamplesGrid+(1:StopI),t)));abs(FIDs(nSamplesGrid+(1:StopI),t))],4); % extrapolate to horizontal
  a = polyval(p, myTime(1:StopI), S, mu);
  [p, S, mu] = polyfit(myTime(nSamplesGrid+(1:StopI)), unwrap(angle(FIDs(nSamplesGrid+(1:StopI),t)),[],1), 2);
  FIDs(1:StopI,t) = a .* exp(1i*polyval(p, myTime(1:StopI), S, mu));
  % figure(251)
  % subplot(2,1,1)
  % plot(myTime(1:nSamplesGrid*4),(abs(FIDs(1:nSamplesGrid*4,t))))
  % subplot(2,1,2)
  % plot(myTime(1:nSamplesGrid*4),(angle(FIDs(1:nSamplesGrid*4,t))))
  % pause
end

if SeqOut.Spectro.ReferenceFID
  % linear phase correction of Ref FID and Smoothing and generate Reference FID
  FIDRefRaw=zeros(size(myTime));
  FIDRefRaw((myTime>=HW.Spectro.ReferenceTime(1))) = interp1(HW.Spectro.ReferenceTime, HW.Spectro.ReferenceFid, myTime(myTime>=HW.Spectro.ReferenceTime(1)));
  nSamplesRefGrid=find(myTime<=HW.Spectro.ReferenceTime(1),1,'last');
  [p,S,mu]=polyfit(myTime(nSamplesRefGrid+(1:nSamplesRefGrid*2)),abs(FIDRefRaw(nSamplesRefGrid+(1:nSamplesRefGrid*2))),2); % extrapolate
  % [p,S,mu]=polyfit([-flipud(myTime(nSamplesRefGrid+(1:nSamplesRefGrid*2)));myTime(nSamplesRefGrid+(1:nSamplesRefGrid*2))],...
  %                  [ flipud(abs(FIDRefRaw(nSamplesRefGrid+(1:nSamplesRefGrid*2))));abs(FIDRefRaw(nSamplesRefGrid+(1:nSamplesRefGrid*2)))],4); % extrapolate to horizontal
  a=polyval(p,myTime(1:nSamplesRefGrid*2),S,mu);
  [p,S,mu]=polyfit(myTime(nSamplesRefGrid+(1:nSamplesRefGrid*2)),unwrap(angle(FIDRefRaw(nSamplesRefGrid+(1:nSamplesRefGrid*2))),[],1),2);
  FIDRefRaw(1:nSamplesRefGrid*2)=a.*exp(1i*polyval(p,myTime(1:nSamplesRefGrid*2),S,mu));

  nFilter=round(numel(FIDRefRaw)/500);
  PhaseCorr1=get_MeanPhaseDiffWeighted(conv(FIDRefRaw(round(numel(FIDRefRaw)./20):round(numel(FIDRefRaw)./5)),RaisedCosine(nFilter)./sum(RaisedCosine(nFilter)),'valid'));
  FIDRef=FIDRefRaw.*exp(-1i*PhaseCorr1.*(1:size(FIDRefRaw(:),1)).');

  RefFIDIn=FIDRef;
  nFilter=cumsum([3,30,300,300,300]);
  for n=nFilter
    FIDRef=conv(RefFIDIn,RaisedCosine(n)./sum(RaisedCosine(n)),'same');
    FIDRef(1:n)=RefFIDIn(1:n);
    FIDRef(n+(1:n))=(RefFIDIn(n+(1:n)).*(n:-1:1).'+FIDRef(n+(1:n)).*(1:n).')./(n+1);
    [p,S,mu]=polyfit(myTime(end-n*2:end),RefFIDIn((end-n*2:end)),4);
    FIDRef(end-n+1:end)=polyval(p,myTime(end-n+1:end),S,mu);
    RefFIDIn=FIDRef;
  end
  PhaseCorr2=get_MeanPhaseDiffWeighted(FIDRef(nFilter(end):end));
  FIDRef=FIDRef.*exp(-1i*PhaseCorr2.*(1:size(FIDRef(:),1)).');


  % find a suitable T2* or choose one and generate ideal FID
  if 0
    startI = 800;
    % fitexp = fit_exp( FIDs(startI:end,1).*exp(-1i.*get_MeanPhaseDiffWeighted(FIDs(startI:end,1)).*(startI:size(FIDs,1)).'),myTime(startI:end,1),1,0,1,0,0);
    fitexp = fit_exp( abs(RefFID(startI:end,1)),myTime(startI:end,1),1,0,0,0,0);
    SeqOut.Spectro.FIDIdealT2Star = fitexp.xminSingle(3);  %  T2*
  end
  % mytime(find(abs(FIDRef)<max(abs(FIDRef)),1,'last'))
  % abs(FIDRef(find(abs(FIDRef)<max(abs(FIDRef)),1,'last')))
  FidIdealStartValue=abs(FIDRef(1)); % start value
  FidIdeal=FidIdealStartValue.*exp(-myTime./SeqOut.Spectro.FIDIdealT2Star);

  %generate Fid Correction window
  FidCorr=FidIdeal./FIDRef;
  if SeqOut.Spectro.plot
    hf1 = figure(1);
    clf(hf1);
    set(hf1, 'Name', 'FidRef FidIdeal and correction');
    axf(1)=subplot(3,1,1,'parent',hf1);
    plot(axf(1),myTime,abs(FIDRefRaw)/HW.RX.AmplitudeUnitScale,myTime,abs(FIDRef)/HW.RX.AmplitudeUnitScale,myTime,abs(FidIdeal)/HW.RX.AmplitudeUnitScale)
    ylabel(axf(1),{'amplitude',['in ' HW.RX.AmplitudeUnit]});
    legend(axf(1),'FIDRefRaw','FidRef','FidIdeal')
    axf(2)=subplot(3,1,2,'parent',hf1);
    plot(axf(2),myTime,angle(FIDRefRaw.*exp(-1i*(PhaseCorr1+PhaseCorr2).*(1:size(FIDRef(:),1)).')),myTime,angle(FIDRef),myTime,angle(FidIdeal),myTime,angle(FidCorr))
    xlabel(axf(2),'time in s')
    ylabel(axf(2),{'phase','in rad'});
    legend(axf(2),'FIDRefRaw','FidRef','FidIdeal','FidCorr')
    axf(3)=subplot(3,1,3,'parent',hf1);
    plot(axf(3),myTime,abs(FidCorr))
    legend(axf(3),'FidCorr')
    ylabel(axf(3),{'correction','factor'});
    linkaxes(axf(1:3),'x');
  end
else
  FidCorr = ones(size(FIDs,1), 1);
  FidIdealStartValue = mean(abs(FIDs(1,:))); % start value
  FidIdeal = FidIdealStartValue .* exp(-myTime./SeqOut.Spectro.FIDIdealT2Star);
end

FidIdealIntegral = sum(abs(FidIdeal));

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

cutFidIndexEnd = find(myTime<=SeqOut.Spectro.FIDWindowDuration, 1, 'last');
w = reshape([SeqOut.Spectro.FIDWindow(cutFidIndexEnd);
             zeros(size(FIDs,1)-cutFidIndexEnd,1)] .* Window_ResolutionEnhanced.*FidCorr, ...
            size(FIDs,1), 1);

FIDsw = FIDs .* repmat(w, 1, size(FIDs,2));

if SeqOut.Spectro.plot
  hf2 = figure(2);
  clf(hf2);
  set(hf2, 'Name', 'FID + Window Function');
  ax2(1) = subplot(4,1,1, 'Parent', hf2);
  plot(ax2(1), repmat(myTime,1,size(FIDs,2)), abs(FIDs/HW.RX.AmplitudeUnitScale), ...
    myTime, abs(w/max(abs(w(:)))*max(abs(FIDs(:))))/HW.RX.AmplitudeUnitScale, 'k--');
  ylabel(ax2(1), {'Raw FIDs', 'amplitude', ['in ' HW.RX.AmplitudeUnit]});
  ax2(2) = subplot(4,1,2, 'Parent', hf2);
  plot(ax2(2), repmat(myTime,1,size(FIDs,2)), angle(FIDs),myTime,angle(w), 'k--');
  ylabel(ax2(2), {'Raw FIDs', 'phase', 'in rad'});
  ax2(3) = subplot(4,1,3, 'Parent', hf2);
  plot(ax2(3), repmat(myTime,1,size(FIDs,2)), abs(FIDsw)/HW.RX.AmplitudeUnitScale, ...
    myTime, FidIdeal/HW.RX.AmplitudeUnitScale, 'r--');
  ylabel(ax2(3), {'window FIDs', 'amplitude', ['in ' HW.RX.AmplitudeUnit]});
  ax2(4) = subplot(4,1,4, 'Parent', hf2);
  plot(ax2(4), repmat(myTime,1,size(FIDs,2)), angle(FIDsw));
  ylabel(ax2(4), {'window FIDs', 'phase', 'in rad'});
  xlabel(ax2(4), 'time in s');
  set(ax2, 'XGrid', 'on', 'YGrid', 'on');
  linkaxes(ax2, 'x');
end

Omega1 = linspace(-SeqOut.Spectro.FrequencyDrift/2,...
  SeqOut.Spectro.FrequencyDrift/2, SeqOut.Spectro.nFrequencyDrift);
FIDswf1 = bsxfun(@times, FIDsw, ...
  reshape(exp(-1i*2*pi*(myTime.^2)*Omega1), [numel(myTime), 1, SeqOut.Spectro.nFrequencyDrift]));
Sf1 = fftshift(fft(FIDswf1, [], 1), 1);
clear FIDswf1
Sf1 = Sf1(round(ppmStartI/SeqOut.Spectro.ZeroFillFactor):round(ppmStopI/SeqOut.Spectro.ZeroFillFactor),:,:);
Pf1 = abs(Sf1).^2 ./ repmat(sum(abs(Sf1(:,:,1)).^2,1), [size(Sf1,1), 1, size(Sf1,3)]);
clear Sf1
Hf1 = -sum(Pf1.*log(Pf1),1);
clear Pf1
[~, Omega1MinI] = min(Hf1, [], 3);
Hf1 = permute(Hf1, [3, 2, 1]);
Omega1Min = Omega1(Omega1MinI);
if any([Omega1MinI==1, Omega1MinI==numel(Omega1)])
  warning('Frequency drift out of correction range, increase Seq.Spectro.FrequencyDrift');
end

% Frequeny drift

Spectra=fftshift(fft(FIDsw.*exp(-1i.*(myTime.^2)*Omega1Min./2),size(FIDsw,1)*SeqOut.Spectro.ZeroFillFactor,1),1);
SpectraRAWw=fftshift(fft(FIDsw,size(FIDsw,1)*SeqOut.Spectro.ZeroFillFactor,1),1);
SpectraRAW=fftshift(fft(FIDs,size(FIDsw,1)*SeqOut.Spectro.ZeroFillFactor,1),1);
Ir=(ppmStartI:ppmStopI).';
if SeqOut.Spectro.plot
  hf3 = figure(3);
  set(hf3, 'Name', 'Spectrum - Frequency Drift');
  ax3(1) = subplot(4,1,1, 'Parent', hf3);
  plot(ax3(1),ppm(Ir), abs(SpectraRAW(Ir,:))/FidIdealIntegral);
  ylabel(ax3(1), ['Spectrum ', HW.RX.AmplitudeName]);
  ax3(2) = subplot(4,1,2, 'Parent', hf3);
  plot(ax3(2), ppm(Ir), abs(SpectraRAWw(Ir,:))/FidIdealIntegral);
  ylabel(ax3(2), 'Spectrum with Window Function');
  ax3(3) = subplot(4,1,3, 'Parent', hf3);
  plot(ax3(3), ppm(Ir), abs(Spectra(Ir,:))/FidIdealIntegral);
  ylabel(ax3(3), 'Spectrum Frequency Drift Corr.');
  xlabel(ax3(3), 'ppm');
  xlim(ax3(3), [SeqOut.Spectro.ppmStop, SeqOut.Spectro.ppmStart]);
  set(ax3(1:3), 'XDir', 'reverse');
  linkaxes(ax3(1:3), 'x');
  ax3(4) = subplot(4,1,4, 'Parent', hf3);
  plot(ax3(4), repmat(Omega1.',1,size(Hf1,2)), Hf1);
  hold(ax3(4), 'all');
  plot(ax3(4), [Omega1(Omega1MinI); NaN(size(Omega1(Omega1MinI)))], ...
               [Hf1(Omega1MinI+(0:size(Hf1,1):(size(Hf1,1)*(size(Hf1,2)-1)))); NaN(size(Omega1(Omega1MinI)))], 'x');
  hold(ax3(4), 'off');
  ylabel(ax3(4), 'Frequency Drift');
  xlabel(ax3(4), 'Hz/s');
end

% Frequency Match

Ir = (ppmStartI:ppmStopI).';
Omega0 = round(linspace(-SeqOut.Spectro.FrequencyMatch/2, SeqOut.Spectro.FrequencyMatch/2, SeqOut.Spectro.nFrequencyMatch)/diff(ppm(1:2)));
Is = bsxfun(@plus, Ir, Omega0);
SpectraR = Spectra(Ir,end);
Pr = repmat(abs(SpectraR).^2./repmat(sum(abs(SpectraR).^2,1), [size(SpectraR,1),1]), [1,size(Spectra,2),size(Is,2)]);
SpectraS = permute(reshape(Spectra(Is(:),:), [size(Is,1),size(Is,2),size(Spectra,2)]), [1,3,2]);
Ps = abs(SpectraS).^2./repmat(sum(abs(SpectraS(:,:,1)).^2,1), [size(SpectraS,1),1,size(SpectraS,3)]);
D = (2-sum((Ps.*Pr).^0.5,1)).^0.5;  % convn!
clear Ps Pr SpectraS SpectraR

[~, Omega0MinI] = min((D), [], 3);
if any([Omega0MinI==1, Omega0MinI==numel(Omega0)])
  warning('Frequency match out of correction range, increase Seq.Spectro.FrequencyMatch');
end
Omega0Min = Omega0(Omega0MinI);
Imin = bsxfun(@plus, Ir, Omega0Min+(0:size(Spectra,1):size(Spectra,1)*size(Omega0Min,2)-1));
SpectraMin = reshape(Spectra(Imin(:)),size(Imin));

Imin = bsxfun(@plus, (1:size(Spectra,1)).', Omega0Min);
Imin(Imin<1) = Imin(Imin<1) + size(Spectra,1);
Imin(Imin>size(Spectra,1)) = Imin(Imin>size(Spectra,1)) - size(Spectra,1);
Imin = bsxfun(@plus, Imin, (0:size(Spectra,1):size(Spectra,1)*size(Omega0Min,2)-1));
SpectraOmega0 = reshape(Spectra(Imin(:)), size(Spectra));
FIDswC = ifft(ifftshift(SpectraOmega0,1), [], 1);

if SeqOut.Spectro.plot
  hf4 = figure(4);
  set(hf4, 'Name', 'Spectrum - Frequency Match');
  ax4(1) = subplot(4,1,1, 'Parent', hf4);
  plot(ax4(1),Omega0*abs(diff(ppm(1:2))),squeeze(D));
  ylabel(ax4(1), 'frequency match');
  ax4(2) = subplot(4,1,2, 'Parent', hf4);
  plot(ax4(2), repmat(ppm(Ir), [1,size(SpectraMin,2)+1]), ...
               [abs(SpectraMin),mean(abs(SpectraMin),2)]/FidIdealIntegral);
  ylabel(ax4(2), 'absolute spectrum');
  ax4(3) = subplot(4,1,3, 'Parent', hf4);
  plot(ax4(3), repmat(ppm(Ir), [1,size(SpectraMin,2)+1]), ...
               [angle(SpectraMin),angle(mean(SpectraMin,2))]/FidIdealIntegral);
  ylabel(ax4(3), 'spectrum phase');
  ax4(4) = subplot(4,1,4, 'Parent', hf4);
  plot(ax4(4), ppm(Ir), [mean(abs(SpectraMin),2), mean(real(SpectraMin),2), mean(imag(SpectraMin),2)]/FidIdealIntegral);
  ylabel(ax4(4), 'spectrum');
  xlabel(ax4(4), 'ppm');
  set(ax4(2:4), 'XDir', 'reverse');
  linkaxes(ax4(2:4), 'x');
end

PhaseOffset = angle(mean(FIDswC(1,:),2));
SpectraMinCP = bsxfun(@times, SpectraMin, exp(-1i.*PhaseOffset));

Spectrum.FID = mean(FIDswC, 2);
Spectrum.PhaseOffset = PhaseOffset;
Spectrum.Time = myTime;
Spectrum.Spectrum = mean(SpectraMinCP,2)/FidIdealIntegral;
Spectrum.FidIdealIntegral = FidIdealIntegral;
Spectrum.ppm = ppm(Ir);
Spectrum.Frequency = Frequency(Ir);

end
