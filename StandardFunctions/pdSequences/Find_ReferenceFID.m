function [HW, mySave] = Find_ReferenceFID(HW, mySave, Seq, SliceSelect)
%% Acquire an FID as a reference for spectroscopy
%
% The used sample should preferably have a long T2* and must have only one peak
% (without any chemical shifts).
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2021 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

%% input check
if nargin<2, mySave = [];       end
if nargin<3, Seq = [];          end
if nargin<4, SliceSelect = [];  end

%% default values
if isemptyfield(Seq, 'SlicePulse'),  Seq.SlicePulse = HW.FindShim.SlicePulse;  end  % 90 degrees pulse function handle
if isemptyfield(Seq, 'FlipPulse'),  Seq.FlipPulse = HW.FindShim.SlicePulse;  end  % 90 degrees pulse function handle
if isemptyfield(Seq, 'useSliceSelect'),  Seq.useSliceSelect= HW.FindShim.useSliceSelect;  end  % use slice gradient

if isemptyfield(SliceSelect, 'alfa'),      SliceSelect.alfa      = HW.FindShim.SliceSelect.alfa;       end  % around x axis
if isemptyfield(SliceSelect, 'phi'),       SliceSelect.phi       = HW.FindShim.SliceSelect.phi;        end  % around y axis
if isemptyfield(SliceSelect, 'theta'),     SliceSelect.theta     = HW.FindShim.SliceSelect.theta;      end  % around z axis
if isemptyfield(SliceSelect, 'CenterRot'), SliceSelect.CenterRot = HW.FindShim.SliceSelect.CenterRot;  end
if isemptyfield(SliceSelect, 'nRead'),     SliceSelect.nRead     = HW.FindShim.SliceSelect.nRead;      end
if isemptyfield(SliceSelect, 'nPhase'),    SliceSelect.nPhase    = HW.FindShim.SliceSelect.nPhase;     end
if isemptyfield(SliceSelect, 'sizeRead'),  SliceSelect.sizeRead  = HW.FindShim.SliceSelect.sizeRead;   end
if isemptyfield(SliceSelect, 'sizePhase'), SliceSelect.sizePhase = HW.FindShim.SliceSelect.sizePhase;  end
if isemptyfield(SliceSelect, 'thickness'), SliceSelect.thickness = HW.FindShim.SliceSelect.thickness;  end
if isempty(SliceSelect.thickness)
  if Seq.useSliceSelect
    SliceSelect.thickness = 0.008;
  end
end
if ~(Seq.useSliceSelect)
  SliceSelect.thickness = 1e12;
end
Seq.SliceSelect=SliceSelect;

%% run measurement
[HW, mySave, SeqOut, SliceSelect] = sequence_Spectrum(HW, mySave, Seq, SliceSelect);
tr = SeqOut.nEchos + 1;
HW.Spectro.ReferenceFid = SeqOut.loopdata.data(:,1,tr);
HW.Spectro.ReferenceTime = SeqOut.loopdata.time_all(:,1,tr);
HW.Spectro.ReferenceStartTime = SeqOut.StartSequenceTime;
HW.Spectro.SliceSelect = SliceSelect;
HW.Spectro.SlicePulse = SeqOut.SlicePulse;
HW.Spectro.FlipPulse = SeqOut.FlipPulse;
HW.Spectro.useSliceSelect = SeqOut.useSliceSelect;


% write new reference FID to file
if (isfield(HW, 'Spectro') || isa(HW, 'PD.HW')) && ...
    ~isemptyfield(HW.Spectro, 'ReferenceFidPath')
  fprintf('Saving Reference FID to file "%s"\n', HW.Spectro.ReferenceFidPath)
  Spectro.ReferenceFid = HW.Spectro.ReferenceFid;
  Spectro.ReferenceTime = HW.Spectro.ReferenceTime;
  Spectro.ReferenceStartTime = HW.Spectro.ReferenceStartTime;
  Spectro.SliceSelect = HW.Spectro.SliceSelect;
  Spectro.SlicePulse = HW.Spectro.SlicePulse;
  Spectro.FlipPulse = HW.Spectro.FlipPulse;
  Spectro.useSliceSelect = HW.Spectro.useSliceSelect;
  save(HW.Spectro.ReferenceFidPath, '-struct', 'Spectro');
else
  warning(['"HW.Spectro.ReferenceFidPath" is not set. Reference FID is ' ...
    'stored only temporary in HW.Spectro.ReferenceFid.' 10 ...
    'It will be lost at next "LoadSystem".']);
end

end
