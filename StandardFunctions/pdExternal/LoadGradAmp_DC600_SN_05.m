HW.Grad.ExtGradSN = 5;

HW.Grad.PowerDown=0;
HW.Grad.PaEnable=1;

% Grad Arrays
HW.Grad.PaCurrentControlled=[1,1,1,1];                  % If current controlled set 1, if voltage controlled set 0

HW.Grad.PaRin=[24e3,24e3,24e3,24e3];                    % INA137 input impedance

HW.Grad.PaOffsetU=[0,0,0,0];                            % Offset voltage
HW.Grad.PaOffsetI=[-0.0021933,0.003974,0.00016173,0.0040016]; % Offset current

HW.Grad.PaUin2PaUout=[0,0,0,0];                         % Voltage Gain
HW.Grad.PaUin2PaIout=([0.33108,0.33731,0.33135,0.33866]-HW.Grad.PaOffsetI)./1; % Input voltage to output current ratio

HW.Grad.PaPmaxInt=[80,80,80,80];                        % Maximum internal power dissipation

HW.Grad.PaRout=[15000,15000,15000,15000];               % Ouput impedance

HW.Grad.tRamp=18e-6;                                    % minimum ramp time
HW.Grad.tEC=50e-6;                                      % Setting time;
HW.Grad.SystemTimeDelay(1:3) = [1.7784e-05   2.0488e-05   2.1112e-05]; % Time delay of grad amp

HW.Grad.Status1=1;                                      % Power supply ok of DC600
HW.Grad.Status2=1;                                      % Graient and temperatur ok of DC600

% Piezo
HW.Grad.LoadRin(4)=100;  %x y z B0 resistance parallel to Piezo
% name for Plot
HW.Grad.Name(4)={'Piezo Current'};
HW.Grad.AmpUnit(4)={'A'};
HW.Grad.AmpUnitScale(4)=1;
% efficence
HW.Grad.LoadIin2Amp(4)=1; % 1 A per Amper
