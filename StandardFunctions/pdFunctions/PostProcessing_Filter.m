function [Seq, data_1D, data_S_AQs_TRs] = PostProcessing_Filter(nargout_set_sequence, Seq, raw_data_AQs, num_averages, Channel, iDevice)
%% Sort raw data into structures and do some post-processing (incl. Frequency Filter)
%
%   [Seq, data_1D, data_S_AQs_TRs] = PostProcessing_Filter(nargout_set_sequence, Seq, raw_data_AQs, num_averages, Channel)
%
% This function uses zeroFill_image as a "band-pass" filter for the measured
% data.
%
% See also: get_data, get_data_1d
%
% ------------------------------------------------------------------------------
% (C) Copyright 2020-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


iAQ = find([Seq.AQ(:).Channel] == Channel & [Seq.AQ(:).Device] == iDevice, 1, 'first');

if isempty(raw_data_AQs)
  % create arrays with NaN
  if all(Seq.AQ(iAQ).fSample(Seq.AQ(iAQ).nSamples(:)>0) == Seq.HW.RX(iDevice).fSample)
    % un-mixed raw data
    raw_data_AQs = NaN(max(Seq.AQ(iAQ).nSamples(:)), sum(Seq.AQ(iAQ).nSamples(:)>0));
  else
    % raw data from mixer
    raw_data_AQs = NaN(max(Seq.AQ(iAQ).nSamples(:)), sum(Seq.AQ(iAQ).nSamples(:)>0)) + 1i*NaN;
  end
end

if nargout_set_sequence >= 3 % calculate data (samples x AQ windows x tReps)
  [Seq, data_S_AQs_TRs] = get_data(Seq, raw_data_AQs, num_averages, Channel, iDevice);

  % Filter data
  if isemptyfield(Seq, 'PostProcessing'), Seq.PostProcessing = struct(); end
  if isemptyfield(Seq.PostProcessing, 'BW'), Seq.PostProcessing.BW = HW.fLarmor/8; end  % frequency bandwidth in Hz
  if isemptyfield(Seq.PostProcessing, 'offset'), Seq.PostProcessing.offset = 0; end  % frequency offset in Hz
  if isemptyfield(Seq.PostProcessing, 'onlyValid'), Seq.PostProcessing.onlyValid = 0; end % remove re-sampled samples at border that are probably invalid
  
  iAQ_dev = find([Seq.AQ(:).Device] == iDevice, 1, 'first');
  for itRep = 1:size(Seq.AQ(iAQ_dev).nSamples, 2)
    for iAQ = 1:size(Seq.AQ(iAQ_dev).nSamples, 1)
      if isnan(Seq.AQ(iAQ_dev).Start(iAQ,itRep))
        continue;
      end
      removeSamples = Seq.PostProcessing.onlyValid*ceil(Seq.AQ(iAQ_dev).fSample(iAQ,itRep)/Seq.PostProcessing.BW);
      
      DataOfEcho = data_S_AQs_TRs.data(1:Seq.AQ(iAQ_dev).nSamples(iAQ,itRep),iAQ,itRep);
      TimeOfEcho = data_S_AQs_TRs.time_all(1:Seq.AQ(iAQ_dev).nSamples(iAQ,itRep),iAQ,itRep);
      % FIXME: Which is the reference sample with fixed phase?
      % DataOfEchoShift = DataOfEcho .* exp(1i*2*pi*Seq.PostProcessing.offset * (TimeOfEcho-Seq.tEcho));
      DataOfEchoShift = DataOfEcho .* exp(1i*2*pi*Seq.PostProcessing.offset * TimeOfEcho);
      DataOfEchoBW = zeroFill_image(DataOfEchoShift, size(DataOfEchoShift), ...
        [2*Seq.PostProcessing.BW/Seq.AQ(iAQ_dev).fSample(iAQ,itRep), Inf, Inf]);
      if 0
        % replace the invalid samples by complex NaN
        data_S_AQs_TRs.data(1:removeSamples,iAQ,itRep) = NaN+1i*NaN;
        data_S_AQs_TRs.data((removeSamples+1):(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-removeSamples),iAQ,itRep) = ...
          DataOfEchoBW((removeSamples+1):(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-removeSamples),iAQ,itRep);
        data_S_AQs_TRs.data(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-removeSamples:end,iAQ,itRep) = NaN+1i*NaN;
      else
        % remove invalid samples and move valid samples to front
        data_S_AQs_TRs.data(1:(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-2*removeSamples),iAQ,itRep) = ...
          DataOfEchoBW((removeSamples+1):(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-removeSamples));
        data_S_AQs_TRs.data(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-2*removeSamples+1:end,iAQ,itRep) = NaN+1i*NaN;
        data_S_AQs_TRs.time_all(1:(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-2*removeSamples),iAQ,itRep) = ...
          TimeOfEcho((removeSamples+1):(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-removeSamples));
        data_S_AQs_TRs.time_all(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-2*removeSamples+1:end,iAQ,itRep) = NaN;
        data_S_AQs_TRs.time_of_tRep(1:(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-2*removeSamples),iAQ,itRep) = ...
          data_S_AQs_TRs.time_of_tRep((removeSamples+1):(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-removeSamples),iAQ,itRep);
        data_S_AQs_TRs.time_of_tRep(Seq.AQ(iAQ_dev).nSamples(iAQ,itRep)-2*removeSamples+1:end,iAQ,itRep) = NaN;
      end
    end
  end
%   % reduce size of resulting array
%   allNaN = all(all(isnan(data_S_AQs_TRs.data), 2), 3);
%   firstAllNaN = find(~allNaN, 1, 'last')+1;
%   data_S_AQs_TRs.data(firstAllNaN:end,:,:) = [];
%   data_S_AQs_TRs.time_all(firstAllNaN:end,:,:) = [];
%   data_S_AQs_TRs.time_of_tRep(firstAllNaN:end,:,:) = [];
else
  data_S_AQs_TRs = struct('time_of_tRep', NaN, 'time_all', NaN, ...
  'data', NaN, 'averages', [], 'f_fft1_data', [], 'fft1_data', [], ...
  'cic_N', NaN, 'cic_M', NaN, 'cic_corr', [], ...
  'Amplitude2Norm', [], 'Amplitude2Uin', [], 'WindowPhaseOffset', []);
end

if nargout_set_sequence >= 4 % also calculate data_1D (all AQ windows separated by NaN)
  [Seq, data_1D] = get_data_1D(Seq, data_S_AQs_TRs, Channel, iDevice);
else
  data_1D = struct('time_of_tRep', NaN, 'time_all', NaN, ...
  'data', NaN, 'AqFrequency', NaN);
end

end

%%
% if nargout >=3;
%     aiE=0;
%     if nargout>=4;
%         data_1D.time_of_tRep=zeros(sum(Seq.AQ.nSamples(:))+sum(Seq.AQ.nSamples(:)>0),1);
%         data_1D.time_all=data_1D.time_of_tRep;
%         data_1D.data=data_1D.time_of_tRep;
%         data_1D.AqFrequency=data_1D.time_of_tRep;
%     end
%     data_S_AQs_TRs.time_of_tRep=nan(max(Seq.AQ.nSamples(:)),size(Seq.AQ.Start,1),size(Seq.AQ.Start,2));
%     data_S_AQs_TRs.time_all=data_S_AQs_TRs.time_of_tRep;
%     if Seq.CalculateFFTOfData
%         data_S_AQs_TRs.f_fft1_data=data_S_AQs_TRs.time_of_tRep;
%         data_S_AQs_TRs.fft1_data=data_S_AQs_TRs.time_of_tRep;
%     end
%     data_S_AQs_TRs.data=data_S_AQs_TRs.time_of_tRep;
%     if Seq.CalculateFFTOfData;  data_S_AQs_TRs.cic_corr=CIC_corr(Seq.AQ, Seq.HW); end;
%     data_S_AQs_TRs.Amplitude2Norm=1./Seq.AQ.Norm2Amplitude;
%     data_S_AQs_TRs.Amplitude2Uin=Seq.AQ.Amplitude2Uin;
%
%     data_S_AQs_TRs.WindowPhaseOffset=nan(size(Seq.AQ.Start,1),size(Seq.AQ.Start,2));
%     if Seq.CorrectAQWindowPhase
%         Seq.AQ.WindowPhaseCorrection.PhaseIncrementfLarmor=int64(Seq.HW.fLarmor/Seq.HW.RX.fSample*2^Seq.HW.RX.DdsPicBits);
%         Seq.AQ.WindowPhaseCorrection.PhaseIncrementOffset=int64(0);
%         Seq.AQ.WindowPhaseCorrection.LastStart=int64(0);
%     end
%
%     nAQ=0;
%     for TR=1:size(Seq.AQ.Start,2)
%         AQs=find(~isnan(Seq.AQ.Start(:,TR)));
%             if Seq.CorrectAQWindowPhase
%                 if or(Seq.AQ.WindowPhaseCorrection.ResetPhases(TR),TR==1);
%                     Seq.AQ.WindowPhaseCorrection.PhaseOffsetDDS=int64(0);
%                     Seq.AQ.WindowPhaseCorrection.PhaseIncrementOffset=int64(Seq.AQ.Frequency(1,TR)/Seq.HW.RX.fSample*2^Seq.HW.RX.DdsPicBits)-Seq.AQ.WindowPhaseCorrection.PhaseIncrementfLarmor;
%                     Seq.AQ.WindowPhaseCorrection.LastStart=int64(0);
%                 end
%             end
%             for t=AQs.'
%                 nAQ=nAQ+1;
%                 if Seq.CorrectAQWindowPhase
%                     Seq.AQ.WindowPhaseCorrection.NextStart=int64(Seq.AQ.WindowPhaseCorrection.StartTimes(t,TR)*Seq.HW.RX.fSample);
%                     Seq.AQ.WindowPhaseCorrection.PhaseOffsetDDS=mod(Seq.AQ.WindowPhaseCorrection.PhaseOffsetDDS+Seq.AQ.WindowPhaseCorrection.PhaseIncrementOffset*(Seq.AQ.WindowPhaseCorrection.NextStart-Seq.AQ.WindowPhaseCorrection.LastStart),2^Seq.HW.RX.DdsPicBits);
%                     Seq.AQ.WindowPhaseCorrection.LastStart=Seq.AQ.WindowPhaseCorrection.NextStart;
%                     Seq.AQ.WindowPhaseCorrection.PhaseIncrementOffset=mod(int64(Seq.AQ.Frequency(t,TR)/Seq.HW.RX.fSample*2^Seq.HW.RX.DdsPicBits)-Seq.AQ.WindowPhaseCorrection.PhaseIncrementfLarmor,2^Seq.HW.RX.DdsPicBits);
%                     data_S_AQs_TRs.WindowPhaseOffset(t,TR)=mod(2*pi*double(Seq.AQ.WindowPhaseCorrection.PhaseOffsetDDS)/2^Seq.HW.RX.DdsPicBits+pi,2*pi)-pi;
% %                     data_S_AQs_TRs.data(1:Seq.AQ.nSamples(t,TR),t,TR)=raw_data_AQs(1:Seq.AQ.nSamples(t,TR),nAQ).*Seq.AQ.Norm2Amplitude(TR)./Seq.AQ.AmplitudeCalGain(t,TR).*Seq.AQ.rawData2Norm*exp(-1i*data_S_AQs_TRs.WindowPhaseOffset(t,TR));
%                     data_S_AQs_TRs.data(1:Seq.AQ.nSamples(t,TR),t,TR)=raw_data_AQs(1:Seq.AQ.nSamples(t,TR),nAQ).*exp(-1i*data_S_AQs_TRs.WindowPhaseOffset(t,TR));
%                 else
% %                     data_S_AQs_TRs.data(1:Seq.AQ.nSamples(t,TR),t,TR)=raw_data_AQs(1:Seq.AQ.nSamples(t,TR),nAQ).*Seq.AQ.Norm2Amplitude(TR)./Seq.AQ.AmplitudeCalGain(t,TR).*Seq.AQ.rawData2Norm;
%                     data_S_AQs_TRs.data(1:Seq.AQ.nSamples(t,TR),t,TR)=raw_data_AQs(1:Seq.AQ.nSamples(t,TR),nAQ);
%                 end
%                 if Seq.AQ.fSample(t,TR)==Seq.HW.RX.fSample
%                     data_S_AQs_TRs.time_of_tRep(1:Seq.AQ.nSamples(t,TR),t,TR)=linspace(Seq.AQ.Start(t,TR),Seq.AQ.Start(t,TR)+ (Seq.AQ.nSamples(t,TR)-1)/Seq.AQ.fSample(t,TR), Seq.AQ.nSamples(t,TR)).';
%                 else
%                     data_S_AQs_TRs.time_of_tRep(1:Seq.AQ.nSamples(t,TR),t,TR)=linspace(Seq.AQ.Start(t,TR)+0.5/Seq.AQ.fSample(t,TR),Seq.AQ.Start(t,TR)+0.5/Seq.AQ.fSample(t,TR)+ (Seq.AQ.nSamples(t,TR)-1)/Seq.AQ.fSample(t,TR), Seq.AQ.nSamples(t,TR)).';
%                 end
%                 if Seq.CalculateFFTOfData
%                     if rem(Seq.AQ.nSamples(t,TR),2)
%                         data_S_AQs_TRs.f_fft1_data(1:Seq.AQ.nSamples(t,TR),t,TR)=linspace(-Seq.AQ.fSample(t,TR)/2+Seq.AQ.fSample(t,TR)/Seq.AQ.nSamples(t,TR),Seq.AQ.fSample(t,TR)/2-Seq.AQ.fSample(t,TR)/Seq.AQ.nSamples(t,TR), Seq.AQ.nSamples(t,TR)).'+Seq.AQ.Frequency(t,TR);
%                     else
%                         data_S_AQs_TRs.f_fft1_data(1:Seq.AQ.nSamples(t,TR),t,TR)=linspace(-Seq.AQ.fSample(t,TR)/2,Seq.AQ.fSample(t,TR)/2-Seq.AQ.fSample(t,TR)/Seq.AQ.nSamples(t,TR), Seq.AQ.nSamples(t,TR)).'+Seq.AQ.Frequency(t,TR);
%     %                   data_S_AQs_TRs.f_fft1_data(Seq.AQ.nSamples(t,TR):end,t,TR)=nan;%Seq.AQ.fSample(TR)/2-Seq.AQ.fSample(TR)/Seq.AQ.nSamples(t,TR);
%                     end
%                     data_S_AQs_TRs.fft1_data(1:Seq.AQ.nSamples(t,TR),t,TR)=(fftshift(ifft(ifftshift(data_S_AQs_TRs.data(1:Seq.AQ.nSamples(t,TR),t,TR))))).*data_S_AQs_TRs.cic_corr(1:Seq.AQ.nSamples(t,TR),t,TR);
%                 end
%             end
%
%
%             if Seq.CorrectAQWindowPhase
%                 if AQs==0
%                     Seq.AQ.WindowPhaseCorrection.PhaseOffsetDDS=mod(Seq.AQ.WindowPhaseCorrection.PhaseOffsetDDS+Seq.AQ.WindowPhaseCorrection.PhaseIncrementOffset*int64(Seq.tRep(TR)*Seq.HW.RX.fSample),2^Seq.HW.RX.DdsPicBits);
%                 else
%                     Seq.AQ.WindowPhaseCorrection.PhaseOffsetDDS=mod(Seq.AQ.WindowPhaseCorrection.PhaseOffsetDDS+Seq.AQ.WindowPhaseCorrection.PhaseIncrementOffset*(int64(Seq.tRep(TR)*Seq.HW.RX.fSample)-Seq.AQ.WindowPhaseCorrection.LastStart),2^Seq.HW.RX.DdsPicBits);
%                 end
%                 Seq.AQ.WindowPhaseCorrection.LastStart=int64(0);
%             end
% %         end
%     end
%     data_S_AQs_TRs.time_all=data_S_AQs_TRs.time_of_tRep+repmat(reshape(cumsum([0,Seq.tRep(1:end-1)],2),1,1,[]),[size(data_S_AQs_TRs.time_of_tRep,1),size(data_S_AQs_TRs.time_of_tRep,2),1]);
% end
%
% if nargout>=4;
%     nAQ=0;
%     for TR=1:size(Seq.AQ.Start,2)
%         AQs=find(~isnan(Seq.AQ.Start(:,TR)));
%         for t=AQs.'
%             nAQ=nAQ+1;
%             aiE=aiE+Seq.AQ.nSamples(t,TR)+1;
%             aiS=aiE-Seq.AQ.nSamples(t,TR);
%                 data_1D.time_of_tRep(aiS:aiE)=[data_S_AQs_TRs.time_of_tRep(1:Seq.AQ.nSamples(t,TR),t,TR);nan];
%                 data_1D.time_all(aiS:aiE)=[data_S_AQs_TRs.time_all(1:Seq.AQ.nSamples(t,TR),t,TR);nan];
%                 data_1D.data(aiS:aiE)=[data_S_AQs_TRs.data(1:Seq.AQ.nSamples(t,TR),t,TR);nan];
%                 data_1D.AqFrequency(aiS:aiE)=[zeros(Seq.AQ.nSamples(t,TR),1)+Seq.AQ.Frequency(t,TR);nan];
%
%         end
%     end
% end
