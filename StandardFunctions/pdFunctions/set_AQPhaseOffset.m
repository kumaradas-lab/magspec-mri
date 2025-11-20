function HW = set_AQPhaseOffset(HW, AQPhaseOffset)
% set AQ Phase Offset in DEG
for iDevice = numel(HW.MMRT):-1:1
  if (isa(HW.RX(iDevice), 'PD.RX') || isfield(HW.RX(iDevice), 'Calibration')) ...
      && ~isempty(HW.RX(iDevice).Calibration) ...
      && isfield(HW.RX(iDevice).Calibration, 'Frequency') ...
      && ~isempty(HW.RX(iDevice).Calibration.Frequency)
    HW.RX.Calibration.Phase=HW.RX.Calibration.Phase-AQPhaseOffset;
  else
    HW.RX.Calibration.Frequency=(-2:2)*HW.RX.fSample;
    HW.RX.Calibration.Gain=ones(size(HW.RX.Calibration.Frequency));
    HW.RX.Calibration.Phase=-AQPhaseOffset*ones(size(HW.RX.Calibration.Frequency));
  end
end
end