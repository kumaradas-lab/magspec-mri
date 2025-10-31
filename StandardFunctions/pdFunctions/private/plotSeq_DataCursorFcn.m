function outputTxt = plotSeq_DataCursorFcn(pointDataTip, eventObj, HW, iDevice)
%% Custom data cursor function for plotSeq
%
% pointDataTip A PointDataTip object (undocumented Matlab and currently not used)
% eventObj     Object containing event data structure
% outputTxt    Data cursor text (string or cell array of strings)
%
% ------------------------------------------------------------------------------
% (C) Copyright 2017-2021 Pure Devices GmbH, Wuerzburg, Germany
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

  case {'plotSeqLineTX1', 'plotSeqLineTX2', 'plotSeqLineTX3', 'plotSeqLineTX4'}
    % FIXME: Will there ever be more than 2 TX?
    ud = get(target, 'UserData');
    idx = get(eventObj, 'DataIndex');
    if isfield(ud, 'channel')
      % mixed TX mode is supported
      txNum = ud.channel(min(idx, end));
    else
      txNum = str2double(targetTag(end));
    end
    outputTxt{1} = sprintf('%s (channel %d)', HW.TX(iDevice).AmplitudeName, txNum);
    outputTxt{2} = ['time: ' num2str(pos(1), 6) ' s'];
    outputTxt{3} = ['value: ' num2str(pos(2), 3) ' ' HW.TX(iDevice).AmplitudeUnit ...
      ' (amp: ' num2str(ud.abs(idx), 3) ' ' HW.TX(iDevice).AmplitudeUnit ' ' char(8793) ' ' ...
      num2str(ud.abs(idx)*HW.TX(iDevice).AmplitudeUnitScale/HW.TX(iDevice).PaUout2Amplitude(txNum), 4) ' V)'];
    outputTxt{4} = ['frequency: ' num2str(ud.freq(idx)/1e6, 6) ' MHz'];
    if abs(ud.abs(idx)) > eps
      outputTxt{5} = ['phase: ' num2str(ud.phase(idx)/pi*180, 4) char(176)];
    end
    if abs(pos(2)) > eps
      time = get(target, 'XData');
      value = get(target, 'YData');
      neighbors = idx+[-1,1];
      neighbors(neighbors < 0) = [];
      neighbors(neighbors > numel(value)) = [];
      neighbors = neighbors(value(neighbors) == pos(2));
      neighbor = neighbors(time(neighbors) ~= pos(1));
      if numel(neighbor) == 1
        duration = abs(time(neighbor) - pos(1));
        outputTxt{end+1} = ['duration: ' num2str(duration*1e6, 3) ' ' char(181) 's'];
        flipAngle = duration * ud.abs(idx)*HW.TX(iDevice).AmplitudeUnitScale * HW.GammaDef * 360/(2*pi);
        outputTxt{end+1} = ['flip angle: ' num2str(flipAngle, 3) char(176)];
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
    SamplingFactor = ud.SamplingFactor(idx);
    if SamplingFactor ~= 1
      outputTxt{end+1} = ['sampling factor: ' num2str(SamplingFactor, '%d')];
    end
    duration = nSamples/fSample/HW.RX(iDevice).fSampleUnitScale;
    outputTxt{end+1} = ['duration: ' num2str(duration*1e3) ' ms'];
    outputTxt{end+1} = ['BW per pixel: ' num2str(1/duration) ' Hz'];
    outputTxt{end+1} = ['frequency: ' num2str(ud.Frequency(idx)/1e6, 6) ' MHz'];
    outputTxt{end+1} = ['phase: ' num2str(mod(ud.Phase(idx), 360), 4) char(176)];

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
