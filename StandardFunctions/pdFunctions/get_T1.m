function [t1, T1, data, SeqOut] = get_T1(HW, Seq)
% [t1, T1, data, SeqOut] = get_T1(HW, Seq)
% measures T1 value using sequence_Recovery.m
%
% This function is obsolete and will be removed in a future version of
% OpenMatlab. Please, call "sequence_Recovery" directly.

% % Create standard parameters if missing
% if nargin == 1
%   Seq.Plot=0;
% end
% if isemptyfield(Seq, 'tRelax')
%   Seq.tRelax = 0;
% end
% if isemptyfield(Seq, 'Flip')
%   Seq.Flip = 5;
% end
% if isemptyfield(Seq, 'nFids')
%   Seq.nFids = 50;
% end
% if isemptyfield(Seq, 'FlipPulse')
%   Seq.FlipPulse = @Pulse_Rect;
% end
% if isemptyfield(Seq, 'InvertPulse')
%   Seq.InvertPulse = @Pulse_Rect_Composite180;
% end
% if isemptyfield(Seq, 'Plot')
%   Seq.Plot = 0;
% end
% if isemptyfield(Seq, 'ConsoleOut')
%   Seq.ConsoleOut  = 0;
% end
% if isemptyfield(Seq, 'tFlip')
%   Seq.tFlip = 10e-3;
% end
% if isemptyfield(Seq, 'tFlipStart')
%   Seq.tFlipStart =  eq.tFlip;
% end
% if isemptyfield(Seq, 'tFlipEnd')
%   Seq.tFlipEnd = Seq.tFlip*Seq.nFids;
% end
% if isemptyfield(Seq, 'tFlipLog')
%   Seq.tFlipLog = 0;
% end

warning('PD:obsolete_function', ['This function is obsolete and will be removed ' ...
  'in a future version of OpenMatlab. Please, call "sequence_Recovery" directly.']);

% Start measurement by calling the pre-made sequence
[t1, T1, data, SeqOut] = sequence_Recovery(HW, Seq);


%% -----------------------------------------------------------------------------
% (C) Copyright 2012-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------
