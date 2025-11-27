function outputTxt = plotSeq_DataCursorFcn(pointDataTip, eventObj, HW, iDevice)
%% Custom data cursor function for plotSeq
%
% pointDataTip A PointDataTip object (undocumented Matlab and currently not used)
% eventObj     Object containing event data structure
% outputTxt    Data cursor text (string or cell array of strings)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2017-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

pos = get(eventObj, 'Position');
target = get(eventObj, 'Target');

targetTag = get(target, 'Tag');
switch targetTag
  case {'plotSeqLineGrad1', 'plotSeqLineGrad2', 'plotSeqLineGrad3', 'plotSeqLineGrad4'}
    gradNum = str2double(targetTag(end));
    outputTxt{1} = HW.Grad(iDevice).Name{gradNum};
    outputTxt{2} = ['time: ' num2str(pos(1), 6) ' s'];
    outputTxt{3} = ['amp: ' num2str(pos(2), 3) ' ' HW.Grad(iDevice).AmpUnit{gradNum}];

  case {'plotSeqLineTX1', 'plotSeqLineTX2', 'plotSeqLineTX3', 'plotSeqLineTX4', ...
      'plotSeqLineTX1DC', 'plotSeqLineTX2DC', 'plotSeqLineTX3DC', 'plotSeqLineTX4DC'}
    % FIXME: Will there ever be more than 2 TX?
    ud = get(target, 'UserData');
    idx = get(eventObj, 'DataIndex');
    if isfield(ud, 'channel')
      % mixed TX mode is supported
      txNum = ud.channel(min(idx, end));
    else
      txNum = str2double(targetTag(end));
    end
    if strcmp(targetTag(end-1:end), 'DC')
      isAmpDC = true;
      ampDC_str = ' DC';
    else
      isAmpDC = false;
      ampDC_str = '';
    end
    outputTxt{1} = sprintf('%s (channel %d)%s', HW.TX(iDevice).AmplitudeName, txNum, ampDC_str);
    outputTxt{2} = ['time: ' num2str(pos(1), 6) ' s'];
    TXCal = PD.TXMaxDef(HW.TX(iDevice), 'Cal', ud.freq(idx));
    TXCal.Amplitude(txNum) = ud.abs(idx)*HW.TX(iDevice).AmplitudeUnitScale;
    if isAmpDC
      outputTxt{3} = ['value: ' num2str(pos(2), 3) ' ' HW.TX(iDevice).AmplitudeUnit ...
        ' (' char(8793) ' ' num2str(TXCal.PaUoutCalibrated(txNum)*sign(pos(2)), 4) ' V)'];
    else
      outputTxt{3} = ['value: ' num2str(pos(2), 3) ' ' HW.TX(iDevice).AmplitudeUnit ...
        ' (amp: ' num2str(ud.abs(idx), 3) ' ' HW.TX(iDevice).AmplitudeUnit ' ' char(8793) ' ' ...
        num2str(TXCal.PaUoutCalibrated(txNum), 4) ' V)'];
      if ud.freq(idx) < 1e6
        outputTxt{4} = ['frequency: ' num2str(ud.freq(idx)/1e3, '%.6f') ' kHz'];
      else
        outputTxt{4} = ['frequency: ' num2str(ud.freq(idx)/1e6, '%.6f') ' MHz'];
      end
      if abs(ud.abs(idx)) > eps
        outputTxt{5} = ['phase: ' num2str(ud.phase(idx)/pi*180, 4) char(176)];
      end
    end
    if abs(pos(2)) > eps
      time = get(target, 'XData');
      value = get(target, 'YData');
      neighbors = idx+[-1,1];
      neighbors(neighbors < 0) = [];
      neighbors(neighbors > numel(value)) = [];
      neighbors = neighbors(value(neighbors) == pos(2));
      neighbor = neighbors(time(neighbors) ~= pos(1));
      if isscalar(neighbor)
        duration = abs(time(neighbor) - pos(1));
        outputTxt{end+1} = ['duration: ' num2str(duration*1e6, '%.3f') ' ' char(181) 's'];
        if ~isAmpDC
          if isempty(HW.TX(iDevice).PaUout2AmplitudeX) || isequal(HW.fLarmor, HW.fLarmorX)
            useXNucleus = false;
          else
            useXNucleus = abs(ud.freq(idx) - HW.fLarmorX) < abs(ud.freq(idx) - HW.fLarmor);
          end
          if useXNucleus
            flipAngle = duration * ud.abs(idx)*HW.TX(iDevice).AmplitudeUnitScale * HW.GammaX * 360/(2*pi);
          else
            flipAngle = duration * ud.abs(idx)*HW.TX(iDevice).AmplitudeUnitScale * HW.GammaDef * 360/(2*pi);
          end
          outputTxt{end+1} = ['flip angle: ' num2str(flipAngle, 3) char(176)];
        end
      end
    end

  case 'plotSeqLineAQ'
    if numel(pos) < 3
      fSample = pos(2);
    else
      fSample = pos(3);
    end
    outputTxt{1} = ['time: ' num2str(pos(1), 6) ' s'];
    outputTxt{end+1} = ['sampling rate: ' num2str(fSample, 6) ' ' HW.RX(iDevice).fSampleUnit];
    % get number of samples from time
    ud = get(target, 'UserData');
    idx = find(ud.tStart <= pos(1), 1, 'last');
    nSamples = ud.nSamples(idx);
    if nSamples == 1
      outputTxt{end+1} = '1 sample';
    else
      outputTxt{end+1} = [num2str(nSamples, '%d') ' samples'];
    end
    if isfield(ud, 'SamplingFactor')
      SamplingFactor = ud.SamplingFactor(idx);
    else
      SamplingFactor = 1;
    end
    if SamplingFactor ~= 1
      outputTxt{end+1} = ['sampling factor: ' num2str(SamplingFactor, '%d')];
    end
    duration = nSamples/fSample/HW.RX(iDevice).fSampleUnitScale;
    outputTxt{end+1} = ['duration: ' num2str(duration*1e3) ' ms'];
    outputTxt{end+1} = ['BW per pixel: ' num2str(1/duration) ' Hz'];
    if ud.Frequency(idx) < 1e6
      outputTxt{end+1} = ['frequency: ' num2str(ud.Frequency(idx)/1e3, '%.6f') ' kHz'];
    else
      outputTxt{end+1} = ['frequency: ' num2str(ud.Frequency(idx)/1e6, '%.6f') ' MHz'];
    end
    outputTxt{end+1} = ['phase: ' num2str(mod(ud.Phase(idx), 360), 4) char(176)];
    if isfield(ud, 'FrequencyX') && isfield(ud, 'PhaseX') && ...
        isfinite(ud.FrequencyX(idx)) && isfinite(ud.PhaseX(idx))
      if ud.FrequencyX(idx) < 1e6
        outputTxt{end+1} = ['frequency X: ' num2str(ud.FrequencyX(idx)/1e3, '%.6f') ' kHz'];
      else
        outputTxt{end+1} = ['frequency X: ' num2str(ud.FrequencyX(idx)/1e6, '%.6f') ' MHz'];
      end
      outputTxt{end+1} = ['phase X: ' num2str(mod(ud.PhaseX(idx), 360), 4) char(176)];
    end

  case 'plotSeqLineDigiIO'
    if rem(pos(2), 1) > 0.5
      state = 'off';
    else
      state = 'on';
    end
    outputTxt{1} = ['time: ' num2str(pos(1), 6) ' s'];
    outputTxt{2} = sprintf('DigiIO%d: %s', round(pos(2)), state);

  otherwise
    outputTxt = [];

end

end
