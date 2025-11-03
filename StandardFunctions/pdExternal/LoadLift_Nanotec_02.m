%% Configure and initialize sample lift with S/N 2
if ~exist('PD.SampleLift', 'class')
  return;
end

HW.Lift.SN = 2;

if (~isa(HW, 'PD.HW') && isemptyfield(HW, 'Lift')) || ~isa(HW.Lift, 'PD.SampleLift') || ...
    isempty(HW.Lift) || (isa(HW.Lift, 'PD.SampleLift') && ~isvalid(HW.Lift))
  % only re-create object if it was explicitly deleted before
  if isa(HW, 'PD.HW') || (~isemptyfield(HW, 'Lift') && ~isemptyfield(HW.Lift, 'SN'))
    lift = HW.Lift;
  else
    lift = struct();
  end
  HW.Lift = PD.SampleLift.GetInstance(lift);
end

% define any settings that differ from the default values here
% HW.Lift.SetMicroStepsDivisor(5);
HW.Lift.SetInputPolarityDirect([5 6]);
HW.Lift.startSetChannel = 6;  % This should connect to Digital Output 3
HW.Lift.setBit0Channel = 5;   % This should connect to Digital Output 4
% HW.Lift.directionChannel = 5;
% HW.Lift.clockChannel = 6;
HW.Lift.upperLimitChannel = 1;
HW.Lift.lowerLimitChannel = 2;
HW.Lift.upDirection = 0;


% Set up "AbortMeasurement" function
HW.Function_Abort_Measurement = @(HW) PD.AbortLift();

% Initialize
HW.Lift.ClearErrorMessage();
err = HW.Lift.StartReference();

if err > 0
  if isempty(HW.Lift.errorMessage)
    errorStr = err.MessageStr();
  else
    errorStr = HW.Lift.errorMessage;
  end
  error('PD:SampleLift:Configuration', ...
    'Error %d (%s) during configuration of sample lift: %s', err, err, errorStr);
end

% HW.Lift.UseJoystick();
