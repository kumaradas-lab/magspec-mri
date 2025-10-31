HW.Grad.ExtGradSN = 24;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[-0.00318512,-0.001497788,0.00189265,0.001343368]; % Offset Current 05-Nov-2018 14:53:29

HW.Grad.PaUin2PaIout=([0.3293668,0.331638,0.3347297,0.3344647]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio% Input Voltage to output Current ratio 05-Nov-2018 14:53:29

HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=50e-6;                                    % minimum ramp time in s
HW.Grad.tEC=50e-6;                                      % Setting time;
HW.Grad.SystemTimeDelay(1:3) = [50e-6, 50e-6, 50e-6]; % Time delay of grad amp
HW.Grad.MaxAmpSlice=0.1;                                % Max Grad Amp ?

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Gradient and temperature ok of DC600
