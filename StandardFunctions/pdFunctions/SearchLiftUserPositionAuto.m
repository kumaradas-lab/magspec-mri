function [SeqOut, data] = SearchLiftUserPositionAuto(HW, Seq)
%% Automatically determine the "user zero" position
%
%   [SeqOut, data] = SearchLiftUserPositionAuto(HW, Seq)
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
% (C) Copyright 2018-2023 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% check input
if nargin ~= 2
  error('Correct syntax is: SearchLiftUserPositionAuto(HW, Seq)');
end
if ~(isa(HW, 'PD.HWClass') || ~isemptyfield(HW, 'Lift')) ...
    || ~isa(HW.Lift, 'PD.SampleLift') || ~isvalid(HW.Lift)
  error('You must have a sample lift to use this function');
end

%% default input
if isemptyfield(Seq, 'Lift'), Seq.Lift = struct(); end
if isemptyfield(Seq.Lift, 'SampleHolderOffset'), Seq.Lift.SampleHolderOffset = 0; end
if isemptyfield(Seq.Lift, 'Speed'), Seq.Lift.Speed = HW.Lift.defaultSpeed; end


%% default output
SeqOut = []; data = [];


%% Find start position
HW.Lift.UseJoystick();
Seq.Lift.useLift = false; % Boolean, move lift at end of sequence
Seq.Lift.checkEncoder = false; % Boolean, check whether motor position and encoder position coincide

hMsg = msgbox(sprintf(['Move the sample lift to a position below the desired ' ...
  '"user zero".\n' ...
  'Click ok when ready,']), ...
  'Start position below "user zero"');

hf = figure(301);
clf(hf);
set(hf, 'Name', 'Mean signal');
hax = axes(hf);
xlabel(hax, 'Amplitude');
ylabel(hax, 'Absolute position in m');

lineColors = get(hax, 'ColorOrder');
maxNumColors = 5; % maximum number of colors used in this function
if size(lineColors, 1) < 5
  lineColors = repmat(lineColors, ceil(maxNumColors/size(lineColors, 1)), 1);
end

hl_all = line(hax, 'XData', NaN, 'YData', NaN, 'LineStyle', 'none', ...
  'Marker', 'x', 'Color', lineColors(1,:));
hl_last = line(hax, 'XData', NaN, 'YData', NaN, 'LineStyle', 'none', ...
  'Marker', 'x', 'Color', lineColors(2,:), 'MarkerSize', 10);

posManual = [];
ampManual = [];
Seq.doDisplay = false;
while ishghandle(hMsg, 'figure')
  [data, SeqOut] = sequence_RecoveryCPMG(HW, Seq);
  Seq.LoopCountStart = SeqOut.LoopCountEnd;
  Seq.Reinitialize = 0;

  %% plot results
  [SeqOut, data] = evaluate_RecoveryCPMG_Lift(HW, SeqOut, data);

  posManual(end+1) = HW.Lift.userZero + (1 - 2*HW.Lift.userSampleCoord) * data.lift.position;
  ampManual(end+1) = mean(abs(data.SampleEchoTau1(:)));

  set(hl_all, 'XData', ampManual(1:end-1), 'YData', posManual(1:end-1));
  set(hl_last, 'XData', ampManual(end), 'YData', posManual(end));

end

set(hl_all, 'XData', ampManual, 'YData', posManual);


%% Scan profile coarsely
set(hl_last, 'XData', NaN, 'YData', NaN, 'Marker', '.');

Seq.Lift.useLift = true; % Boolean, move lift at end of sequence
Seq.Lift.TotalDisplacement = Seq.thicknessSlice * 4; % total displacement in m
Seq.Lift.nLoops = 6; % number of steps
Seq.Lift.Displacement = Seq.Lift.TotalDisplacement / (Seq.Lift.nLoops-1); % displacement in m

Seq.plotSeq = 1:4;

HW.Lift.SetSpeed(Seq.Lift.Speed);
HW.Lift.SetRelativePosition(Seq.Lift.Displacement);

HW.Lift.SaveSet(1); % Save to set 1

movementTime = HW.Lift.GetMovementTime();
Seq.tRelax = max(Seq.tRelax, movementTime+0.5);

HW.Lift.UseSetMode();

% activate joystick mode on interrupt
unwind = onCleanup(@() HW.Lift.UseJoystick());

Seq.Lift.checkEncoder = true;

Seq.Reinitialize = 1;
Seq.LoopCountStart = 0;

posScan = NaN(1, Seq.Lift.nLoops);
ampScan = NaN(1, Seq.Lift.nLoops);
for iLiftLoop = 1:Seq.Lift.nLoops
  if iLiftLoop == Seq.Lift.nLoops
    % Don't move after last step
    Seq.Lift.skipLift = true;
  end

  [data, SeqOut] = sequence_RecoveryCPMG(HW, Seq);
  Seq.LoopCountStart = SeqOut.LoopCountEnd;
  Seq.Reinitialize = 0;

  %% plot results
  [SeqOut, data] = evaluate_RecoveryCPMG_Lift(HW, SeqOut, data);

  % convert position from user coordinates to absolute coordinates
  posScan(iLiftLoop) = HW.Lift.userZero + (1 - 2*HW.Lift.userSampleCoord) * data.lift.position;
  ampScan(iLiftLoop) = mean(abs(data.SampleEchoTau1(:)));

  set(hl_last, 'XData', ampScan, 'YData', posScan);

  % FIXME: Detect early when we can stop scanning (i.e. the max is found).
end
pause(Seq.tRelax);


%% Scan with smaller steps around maximum
weightedMeanPos = sum(posScan.*ampScan)/sum(ampScan);
Seq.Lift.useLift = true; % Boolean, move lift at end of sequence
Seq.Lift.skipLift = false;
Seq.Lift.TotalDisplacement = -Seq.thicknessSlice * 2.2; % total displacement in m
Seq.Lift.nLoops = 8; % number of steps
Seq.Lift.Displacement = Seq.Lift.TotalDisplacement / (Seq.Lift.nLoops-1); % displacement in m

HW.Lift.SetRelativePosition(Seq.Lift.Displacement);

HW.Lift.SaveSet(1); % Save to set 1

movementTime = HW.Lift.GetMovementTime();
Seq.tRelax = max(Seq.tRelax, movementTime+0.5);

HW.Lift.MoveAbsolute(weightedMeanPos - Seq.Lift.TotalDisplacement/2, true);
Seq.Reinitialize = 1;
Seq.LoopCountStart = 0;

hl_fine = line(hax, 'XData', NaN, 'YData', NaN, 'LineStyle', 'none', ...
  'Marker', '+', 'Color', lineColors(3,:));

posFine = NaN(1, Seq.Lift.nLoops);
ampFine = NaN(1, Seq.Lift.nLoops);
for iLiftLoop = 1:Seq.Lift.nLoops
  if iLiftLoop == Seq.Lift.nLoops
    % Don't move after last step
    Seq.Lift.skipLift = true;
  end

  [data, SeqOut] = sequence_RecoveryCPMG(HW, Seq);
  Seq.LoopCountStart = SeqOut.LoopCountEnd;
  Seq.Reinitialize = 0;

  %% plot results
  [SeqOut, data] = evaluate_RecoveryCPMG_Lift(HW, SeqOut, data);

  % convert position from user coordinates to absolute coordinates
  posFine(iLiftLoop) = HW.Lift.userZero + (1 - 2*HW.Lift.userSampleCoord) * data.lift.position;
  ampFine(iLiftLoop) = mean(abs(data.SampleEchoTau1(:)));

  set(hl_fine, 'XData', ampFine, 'YData', posFine);
end


%% fit sinc function to profile
startParam = [weightedMeanPos, max(ampFine), 4*Seq.thicknessSlice]; % posMax, ampMax, width
[fitParam, data.fval, data.exitflag, data.output] = ...
  fminsearch(@(x) fitSinc(x, posFine, ampFine), startParam);

posFit = linspace(0, 3*fitParam(3), 100) + fitParam(1) - fitParam(3)*3/2;
line(hax, 'XData', innerSinc(posFit, fitParam(1), fitParam(2), fitParam(3)), ...
  'YData', posFit, 'Color', lineColors(4,:));
line(hax, 'XData', fitParam(2), 'YData', fitParam(1), 'LineStyle', 'none', ...
  'Marker', 'o', 'Color', lineColors(5,:));
line(hax, 'XData', [0 fitParam(2)], 'YData', [fitParam(1) fitParam(1)], ...
  'LineStyle', '-', 'Color', lineColors(5,:));


%% set new value for "user zero"
figure(hf); % bring figure with data to front
doSet = questdlg(sprintf(['Do you want to use %.1f mm as the reference for the ' ...
  'new "user zero" position?'], fitParam(1)*1e3), ...
  'New user reference position', 'Yes', 'No', 'Yes');

if ~strcmp(doSet, 'Yes'), return; end

err = HW.Lift.SetUserZero(fitParam(1), Seq.Lift.SampleHolderOffset);

if err
  errordlg(sprintf('Error %s setting the new "zero" position.', err));
else
  msgbox(sprintf(['New "user zero" position successfully set to %.4f m.\n', ...
    'Maximum sample height: %.1f mm'], HW.Lift.userZero, (HW.Lift.userZero - HW.Lift.minPosition)*1e3), ...
    'New "user zero" position');
end


end


function residual = fitSinc(fitParam, posMeas, ampMeas)
%% Objective function for fitting the inner lobe os a sinc function

sinc = innerSinc(posMeas, fitParam(1), fitParam(2), fitParam(3));
numPointsInRange = sum(~isnan(sinc));
if numPointsInRange < 3
  % This is necessary to not move the data completely out of range of the inner
  % lobe of the sinc function.
  residual = fitParam(2) * (3 - numPointsInRange) * 1e3;
else
  residual = sqrt(sum((ampMeas - sinc).^2, 'omitnan')/numPointsInRange);
end

end


function val = innerSinc(x, pos, amp, width)
%% Inner lobe of sinc function (NaN otherwise)

posScaled = (x-pos) * pi / (width / 2);

val = amp * sin(posScaled)./posScaled;

% FIXME: Is this good enough? Or is the problem with 0/0 emerging for larger
% values already?
val(abs(posScaled) < eps) = amp;
% "kill" the outer lobes
val(abs(posScaled) > pi) = NaN;

end
