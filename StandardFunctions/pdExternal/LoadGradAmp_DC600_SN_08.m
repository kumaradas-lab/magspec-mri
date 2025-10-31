HW.Grad.ExtGradSN = 8;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[0.002781935,0.0007123647,-0.004950025,-0.002564188]; % Offset Current 11-Nov-2016 16:11:09

HW.Grad.PaUin2PaIout=([0.3359287,0.3340262,0.3302193,0.3319968]-HW.Grad.PaOffsetI)./1; % Input Voltage to output Current ratio% Input Voltage to output Current ratio 11-Nov-2016 16:12:48

HW.Grad.PaPmaxInt=[100,100,100,100];                    % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Output impedance

HW.Grad.tRamp=50e-6;                                    % minimum ramp time
HW.Grad.tEC=50e-6;                                      % Setting time;
HW.Grad.SystemTimeDelay(1:3) = [1.7648e-05   1.9408e-05    1.756e-05]; % Time delay of grad amp
HW.Grad.MaxAmpSlice=0.1;                                % Max Grad Amp ?

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Gradient and temperature ok of DC600
