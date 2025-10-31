classdef TX < handleHidden
%% TX  Class storing properties for transmission to coil
%
% constructor:
%   TXobj = PD.TX(HW)
%
% Input:
%   HW      object of class PD.HW
%
% TX properties:
%   fSample             sample frequency of HF-DAC
%   DdsPicBits          bitwidth of DDS PIC register
%   DdsPofBits          bitwidth of DDS POF register
%   DacBits             bitwidth of TX-DAC
%   ExtRFSN             serial number of external RF amplifier
%
%   n                   number of TX channels
%
%   Latenz              latency of TX channels in seconds
%   BlankPostsetLatenz  latency of blank (falling edge?) in seconds
%   BlankOffsetLatenz   latency of blank (rising edge?) in seconds
%
%   Dac2Norm            voltage step when increasing value at digital
%                       output by one step (in Volts/step). This value might
%                       be unreliable for calibrated amplifiers. Consider
%                       using "get_TX_Amplitude" instead.
%
%   ChannelDef          index of default TX channel
%   BlankOffset         offset of blanking signal before TX pulse in
%                       seconds
%   BlankPostset        postset of blanking signal after TX pulse in
%                       seconds
%   Rout                output impedance of HF amplifier (???)
%
%   PaUout2Amplitude    conversion factor from voltage amplitude at
%                       inlet of coil to B1 field amplitude in T/V
%   PaUout2AmplitudePath path to .m-file with those factors
%
%   AmplitudeName       string with name of the amplitude (for e.g.
%                       plots)
%   AmplitudeUnit       string with unit of the amplitude (e.g. 'nT')
%   AmplitudeUnitScale  [double] conversion factor for unit
%                       (signal in Tesla = AmplitudeUnitScale × signal
%                       in AmplitudeUnit)
%   DampCoil            structure with setting for coil damping
%
%
% TX methods:
%   dBm2Amp
%
%
% See also:
%   HW_Standard
%
% ----------------------------------------------------------------------------
% (C) Copyright 2016-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ----------------------------------------------------------------------------

end


%#function LoadCoil
%#function PD.HW
%#function PD.TXMaxDef
%#function handleHidden

