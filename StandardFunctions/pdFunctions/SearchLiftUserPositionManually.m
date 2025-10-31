function [SeqOut, data] = SearchLiftUserPositionManually(HW, Seq)
%% Manually determine the "user zero" position
%
%   SearchLiftUserPositionManually(HW, Seq)
%
% Start a repeated CPMG train measurement (using "sequence_RecoveryCPMG") while
% manually moving the sample lift.
% When the user interrupts the measurement the position of the lift can be used
% for referencing the "user zero" position of the Lift
%
% INPUT:
%   HW
%         HW structure or object
%   Seq
%         Sequence settings. See "sequence_RecoveryCPMG". Additional settings
%         are:
%     Lift.SampleHolderOffset
%           Distance between lower edge of reference sample and actual sample
%           in m. The offset is positive if the "user zero" position should be
%           higher than the lower edge of the reference sample. (Default: 0)
%
% OUTPUT:
%   SeqOut and data of the last measurement (see "sequence_RecoveryCPMG")
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% check input
if nargin ~= 2
  error('Correct syntax is: SearchLiftUserPositionManually(HW, Seq)');
end
if ~(isa(HW, 'PD.HW') || ~isemptyfield(HW, 'Lift')) ...
    || ~isa(HW.Lift, 'PD.SampleLift') || ~isvalid(HW.Lift)
  error('You must have a sample lift to use this function');
end

%% default input
if isemptyfield(Seq, 'Lift'), Seq.Lift = struct(); end
if isemptyfield(Seq.Lift, 'SampleHolderOffset'), Seq.Lift.SampleHolderOffset = 0; end

%%
HW.Lift.UseJoystick();

hMsg = msgbox(sprintf(['Move the sample lift to the desired position with the joystick' ...
  ' while the measurement is running.\n' ...
  'Then close this dialog (or interrupt the measurement with Ctrl-C) to finalize ' ...
  'referencing the "user zero" position.']), ...
  'Manual "user zero" position');

% Use onCleanup guard to trigger referencing the user zero when user interrupts
% measurement.
guard = onCleanup(@() SaveLiftUserPosition(HW, Seq, hMsg));

SeqOut = []; data = [];

hf = figure(301);
clf(hf);
set(hf, 'Name', 'Mean signal');
hax = axes(hf);
xlabel(hax, 'Amplitude');
ylabel(hax, 'Absolute position in m');

lineColors = get(hax, 'ColorOrder');

hl_all = line(hax, 'XData', NaN, 'YData', NaN, 'LineStyle', 'none', ...
  'Marker', 'x', 'Color', lineColors(1,:));
hl_last = line(hax, 'XData', NaN, 'YData', NaN, 'LineStyle', 'none', ...
  'Marker', 'x', 'Color', lineColors(2,:), 'MarkerSize', 10);

pos = [];
amp = [];
Seq.doDisplay = false;
while ishghandle(hMsg, 'figure')
  [data, SeqOut] = sequence_RecoveryCPMG(HW, Seq);
  Seq.LoopCountStart = SeqOut.LoopCountEnd;
  Seq.Reinitialize = 0;

  %% plot results
  [SeqOut, data] = evaluate_RecoveryCPMG_Lift(HW, SeqOut, data);

  pos(end+1) = HW.Lift.userZero + (1 - 2*HW.Lift.userSampleCoord) * data.lift.position;
  amp(end+1) = mean(abs(data.SampleEchoTau1(:)));

  set(hl_all, 'XData', amp(1:end-1), 'YData', pos(1:end-1));
  set(hl_last, 'XData', amp(end), 'YData', pos(end));

end

end


function SaveLiftUserPosition(HW, Seq, hMsg)
%% "Cleanup" function that executes when user interrupts measurement

if ishghandle(hMsg, 'figure'), close(hMsg); end

doSet = questdlg('Do you want to set the current position as new "user zero" position?', ...
  'New user reference position', 'Yes', 'No', 'Yes');

if ~strcmp(doSet, 'Yes'), return; end

err = HW.Lift.SetUserZero([], Seq.Lift.SampleHolderOffset);

if err
  errordlg(sprintf('Error %s setting the new "zero" position.', err));
else
  msgbox(sprintf(['New "zero" position successfully set to %.4f m.\n', ...
    'Maximum sample height: %.1f mm'], HW.Lift.userZero, (HW.Lift.userZero - HW.Lift.minPosition)*1e3), ...
    'New "user zero" position');
end

HW.Lift.UseJoystick();

end
