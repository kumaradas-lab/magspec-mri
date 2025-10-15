additionalGradientResistorOhms = 20; % 0, 20 or 94 when using RJ45 cable attenuator

global talker;
if(additionalGradientResistorOhms>0)
  % Check if gradient amplifier is connected
  if(talker.myMon.DeviceStatus.COUNT_GRADIENT1_FAULTS~=0)
    error('Please disconnect external gradient amplifier when using additional gradient resistors.');
  end;
  
  % Check gradient resistance 
  if all(HW.Grad.LoadRin<10)
    HW.Grad.LoadRin=HW.Grad.LoadRin+additionalGradientResistorOhms;
  else
    error('Unexpected gradient resistance.');
  end
end