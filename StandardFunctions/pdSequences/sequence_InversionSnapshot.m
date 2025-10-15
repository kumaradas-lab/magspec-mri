function [t1, T1, data, SeqOut] = sequence_InversionSnapshot(HW, Seq)
%% Sequence using a inversion snapshot method to measure the T1 time
%
%   [t1, T1, data, SeqOut] = sequence_InversionSnapshot(HW, Seq)
%
% This sequence sets default parameters for an inversion snapshot and re-routes
% to "sequence_Recovery".
% This function is deprecated and might be removed in a future version of
% OpenMatlab. Consider using "sequence_Recovery" directly.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2012-2018 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%%
Seq = set_EmptyField(Seq, 'recovery', 'Inversion');
Seq = set_EmptyField(Seq, 'tRelax', 0);
if Seq.tRelax~=0, Seq.excitationFlipAngle = 90; end
Seq = set_EmptyField(Seq, 'excitationFlipAngle', 5);

warning('PD:sequence_InversionSnapshot:deprecated', ...
  ['This function is deprecated and might be removed in a future version ' ...
   'of OpenMatlab. Consider using "sequence_Recovery" directly.']);

[t1, T1, data, SeqOut] = sequence_Recovery(HW, Seq);

end
