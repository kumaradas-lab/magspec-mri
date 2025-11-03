function HW = Find_ReferenceFID_Press(HW, Seq)
%% Acquire an FID as a reference for spectroscopy
%
% The used sample should preferably have a long T2* and must have only one peak
% (without any chemical shifts).
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018-2020 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

%% input check


%% run measurement
SeqOut = sequence_Spectrum_Press(HW, Seq);
tr = SeqOut.nEchos+1;
HW.Spectro.ReferenceFid = SeqOut.data.data(:,1,tr);
HW.Spectro.ReferenceTime = SeqOut.data.time_all(:,1,tr);
HW.Spectro.ReferenceStartTime = SeqOut.StartSequenceTime;
HW.Spectro.SliceSelect = SeqOut.AQSlice(1);


% write new reference FID to file
if (isfield(HW, 'Spectro') || isa(HW, 'PD.HW')) && ...
    ~isemptyfield(HW.Spectro, 'ReferenceFidPath')
  fprintf('Saving Reference FID to file "%s"\n', HW.Spectro.ReferenceFidPath)
  Spectro.ReferenceFid = HW.Spectro.ReferenceFid;
  Spectro.ReferenceTime = HW.Spectro.ReferenceTime;
  Spectro.ReferenceStartTime = HW.Spectro.ReferenceStartTime;
  Spectro.SliceSelect = HW.Spectro.SliceSelect;
  % Spectro.ReferenceStartTime = HW.Spectro.ReferenceStartTime;
  % Spectro.SliceSelect = HW.Spectro.SliceSelect;
  % Spectro.SlicePulse = HW.Spectro.SlicePulse;
  % Spectro.FlipPulse = HW.Spectro.FlipPulse;
  % Spectro.useSliceSelect = HW.Spectro.useSliceSelect;
  save(HW.Spectro.ReferenceFidPath, '-struct', 'Spectro');
else
  warning('PD:ReferenceFID_Press:NonPermanent', ...
    ['"HW.Spectro.ReferenceFidPath" is not set. Reference FID is ' ...
    'stored only temporarily in HW.Spectro.ReferenceFid.\n' ...
    'It will be lost at next "LoadSystem".']);
end

end
