function Seq = get_LiftMovement(Seq, HW)
%% Generate sequence of spikes at the DigitalIO for sample lift movement
%
%   Seq = get_LiftMovement(Seq, HW)
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2018 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%%
if ~isfield(Seq, 'Lift'), Seq.Lift = []; end
if isemptyfield(Seq.Lift, 'Displacement'), Seq.Lift.Displacement = HW.Lift.defaultDisplacement; end % Displacement in m
if isemptyfield(Seq.Lift, 'Direction') % Direction of movement: 1 is up, 0 is down
  Seq.Lift.Direction = double(Seq.Lift.Displacement >= 0);
else
  Seq.Lift.Direction = double(~xor(Seq.Lift.Direction, Seq.Lift.Displacement >= 0));
end
Seq.Lift.Displacement = abs(Seq.Lift.Displacement);
if isemptyfield(Seq.Lift, 'Speed'), Seq.Lift.Speed = HW.Lift.defaultSpeed; end % Speed of movement in m/s
if isscalar(Seq.Lift.Speed), Seq.Lift.Speed = Seq.Lift.Speed * ones(size(Seq.Lift.Displacement)); end
if isemptyfield(Seq.Lift, 'Start'), Seq.Lift.Start = zeros(size(Seq.Lift.Displacement)); end % Start time of movement in s
if isscalar(Seq.Lift.Start), Seq.Lift.Start = Seq.Lift.Start * ones(size(Seq.Lift.Displacement)); end
% if isemptyfield(Seq.Lift, 'clockInputChannel'), Seq.Lift.clockInputChannel = HW.Lift.clockChannel; end % number of the port at the controller for the clock signal
if isemptyfield(Seq.Lift, 'clockOutputChannel'), Seq.Lift.clockOutputChannel = 5; end % number of the channel at the DigitalIO port of the drive for the clock signal
% if isemptyfield(Seq.Lift, 'directionInputChannel'), Seq.Lift.directionInputChannel = HW.Lift.directionChannel; end % number of the port at the controller for the direction signal
if isemptyfield(Seq.Lift, 'directionOutputChannel'), Seq.Lift.directionOutputChannel = 6; end % number of the channel at the DigitalIO port of the drive for the direction signal


% number of steps for the distance
numSteps = round(Seq.Lift.Displacement / HW.Lift.steps2meter / HW.Lift.microsteps2steps);
if any(sum(numSteps, 1) > 500/2)
  error('PD:get_LiftMovement:IOJam', ...
    'Number of steps per tRep is exceeded. Reduce lift displacement');
end
% frequency for selected speed
freqStepper = Seq.Lift.Speed / HW.Lift.steps2meter / HW.Lift.microsteps2steps;
if any(freqStepper(:) ~= 0 & freqStepper(:) < 1) || any(freqStepper(:) > 1e6) || ... % controller imposed limits
    any(Seq.Lift.Speed(:) ~= 0 & Seq.Lift.Speed(:) < HW.Lift.minSpeed) || any(Seq.Lift.Speed(:) > HW.Lift.maxSpeed) % user imposed limits
  error('PD:get_LiftMovement:frequency', ...
    'Stepper frequency must be between 1Hz and 1MHz (controller) or speed between %f and %f m/s (user). Check lift speed.', ...
    HW.Lift.minSpeed, HW.Lift.maxSpeed);
end
% debounce time
if 1/freqStepper < 4*HW.Lift.debounceTime
  % FIXME: Is 4 a good enough factor?
  error('PD:get_LiftMovement:debounceTime', ...
    'Debounce time of IO inputs (%d ms) must be shorter than 1/4 of the stepper period (%.2f ms).', ...
    HW.Lift.debounceTime*1e3, 1e3/freqStepper);
end
% check limits
totalDist = Seq.Lift.Displacement(:) * (2*Seq.Lift.Direction(:)-1);
absPos = HW.Lift.GetPosition() + cumsum(totalDist);
if any(absPos > HW.Lift.maxPosition) || any(absPos < HW.Lift.minPosition)
  error('PD:get_LiftMovement:posLimits', ...
    'Movement of lift exceeds position limits. Reduce displacements.');
end

TimeOn = repmat((0:max(numSteps(:))-1).', [1, size(numSteps)]);
selector = zeros(max(numSteps(:))+1, size(numSteps, 1), size(numSteps, 2));
[xx, yy] = meshgrid(1:size(numSteps, 2), 1:size(numSteps, 1));
selector(sub2ind(size(selector), numSteps+1, yy, xx)) = 1; % mark first index that is not used with 1
selector = cumsum(selector(1:max(numSteps(:)),:,:), 1); % use cumsum to mark all unused indices
TimeOn(selector>0) = NaN;
TimeOn = bsxfun(@plus, permute(Seq.Lift.Start, [3 1 2]), ...
  bsxfun(@rdivide, TimeOn, permute(freqStepper, [3 1 2])));
Value = ones(size(TimeOn));
TimeOff = TimeOn + bsxfun(@rdivide, ones(size(TimeOn)), permute(freqStepper, [3 1 2]))/2;
Seq.Lift.SetTime = reshape([TimeOn; TimeOff], [], size(Seq.Lift.Displacement, 2));
Value = reshape([Value; 0*Value], [], size(Seq.Lift.Displacement, 2));
[Seq.Lift.SetTime, iSort] = sort(Seq.Lift.SetTime, 1);
Direction = reshape(repmat(permute(Seq.Lift.Direction, [3 1 2]), 2*max(numSteps(:)), 1, 1), ...
  [], size(Seq.Lift.Displacement, 2));
Seq.Lift.SetValue = Value(iSort) * 2^(Seq.Lift.clockOutputChannel-1) + ...
  Direction * 2^(Seq.Lift.directionOutputChannel-1);

% unset direction channel at end of movement
lastIdx = sub2ind(size(Seq.Lift.SetValue), cumsum(numSteps,1)*2, xx);
Seq.Lift.SetValue(lastIdx) = ...
  Seq.Lift.SetValue(lastIdx) - Seq.Lift.Direction * 2^(Seq.Lift.directionOutputChannel-1);

% reduce matrix sizes
lastUseful = find(any(~isnan(Seq.Lift.SetTime), 2), 1, 'last');
Seq.Lift.SetTime = Seq.Lift.SetTime(1:lastUseful,:);
Seq.Lift.SetValue = Seq.Lift.SetValue(1:lastUseful,:);
Seq.Lift.SetValue(isnan(Seq.Lift.SetTime)) = NaN;

end
