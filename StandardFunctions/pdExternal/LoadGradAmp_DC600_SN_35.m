HW.Grad.ExtGradSN = 35;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[0.001402637,0.005817095,-0.003808113,-0.002240312]; % Offset Current 15-Mar-2021 16:13:37
 
HW.Grad.PaUin2PaIout=([0.333951,0.3382715,0.3286292,0.3303183]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio% Input Voltage to output Current ratio 15-Mar-2021 16:13:37

HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=50e-6;                                    % minimum ramp time in s
HW.Grad.tEC=50e-6;                                      % Setting time;

switch HW.UserName
  case 'magnet_01_probe_1H'
    % Magnet 1: 1H Wechselprobenkopf
    % Magnet 165
    HW.Grad.SystemTimeDelay(1:3) = [1.68496e-05   1.9252e-05   1.6668e-05]; % Time delay of grad amp
  case 'magnet_01_probe_31P'
    % Magnet 1: 31P Wechselprobenkopf
    % Magnet 165
    HW.Grad.SystemTimeDelay(1:3) = [1.68496e-05   1.9252e-05   1.6668e-05]; % Time delay of grad amp
  case 'magnet_02'
    % Magnet 2: 10mm Magnet KEIN Wechselprobenkopf
    % Magnet 166
    HW.Grad.SystemTimeDelay(1:3) = [2.27648e-05  2.64584e-05  2.35296e-05]; % Magnet #166
end

HW.Grad.MaxAmpSlice=0.1;                                % Max Grad Amp ?

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Gradient and temperature ok of DC600
