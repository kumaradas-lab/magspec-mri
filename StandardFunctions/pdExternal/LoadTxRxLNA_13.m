%% Settings for LNA S/N 13

HW.RX.LnaSN = 13;

% Resonant coil 10 mm
HW.TX.Uout2PaUout(2)=HW.TX.Uout2PaUout(2)*10^(-0.2/20);
HW.TX.ChannelDef=2;         % default TX rf channel
HW.TX.Def.Uout(2)=3.8;      % default TX Uout in V

HW.RX2TXdeadTime=4e-6;      % Totzeit des Empfängers vor dem Senden
HW.TX2RXdeadTime=50e-6;     % 50µs (langsame Pindioden im LNA 12) Totzeit des Empfängers nach dem Senden 
HW.TX.BlankOffset=3e-6;      % Offset des Blankingsignals vor dem TX Puls
HW.TX.BlankPostset=1e-6;    % Zeitzusatz des Blankingsignals nach dem TX Puls

switch 4
  case 1
    HW.RX.VGAGainDef=HW.RX.VGAGainMax/1; % max Empfangsverstärkung -0 dB
    HW.RX.LNAGain=10^(((-55.34)-(-80))/20); % 24.66 dB gain @ 24.35 MHz F=1.3 dB
  case 4
    HW.RX.VGAGainDef=HW.RX.VGAGainMax/4; % max Empfangsverstärkung -12 dB
    HW.RX.LNAGain=10^(((-55.39)-(-80))/20); % 26.0 dB gain @ 24.35 MHz F=1.4 dB
  case 10
    HW.RX.VGAGainDef=HW.RX.VGAGainMax/10; % max Empfangsverstärkung -20 dB
    HW.RX.LNAGain=10^(((-55.36)-(-80))/20); % 24.64 dB gain @ 24.35 MHz F=1.6 dB
  case 20
    HW.RX.VGAGainDef=HW.RX.VGAGainMax/20; % max Empfangsverstärkung -20 dB
    HW.RX.LNAGain=10^(((-55.30)-(-80))/20); % 24.70 dB gain @ 24.35 MHz F=2.2 dB
  case 40
    HW.RX.VGAGainDef=HW.RX.VGAGainMax/40; % max Empfangsverstärkung -20 dB
    HW.RX.LNAGain=10^(((-55.27)-(-80))/20); % 24.73 dB gain @ 24.35 MHz F=4 dB
  case 100
    HW.RX.VGAGainDef=HW.RX.VGAGainMax/100; % max Empfangsverstärkung -20 dB
    HW.RX.LNAGain=10^(((-55.46)-(-80))/20); % 24.54 dB gain @ 24.35 MHz F=9.4 dB
end
HW.TX.Max.PaUout(2)=min(HW.TX.Max.PaUout(2),100);

%% TRx switch during transmit at Tx2
HW.TX.BlankOffsetAQ=1e-6;        % Blank of receiver before TX pulse
HW.TX.BlankPostsetAQ=4e-6;       % Blank of receiver after TX pulse
HW.TX.BlankAQ=1;                 % Switch TRx to 50 Ohm Resistor during TX pulse, to avoid saturation.
