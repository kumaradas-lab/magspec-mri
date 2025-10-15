function plotSeqTR(HW, Seq, AQ, TX, Grad)
%% Plot sequence with stacked "tRep"s
%
%   plotSeqTR(HW, Seq, AQ, TX, Grad)
%
% This function opens a figure showing the pulse program (rf pulses, acquisition
% windows, gradient pulses and digital IO signals) of the sequence. The "tRep"s
% starting at Seq.plotSeqStart and ending with Seq.plotSeqEnd are stacked upon
% each other such that each tRep begins at t=0s.
%
% INPUT:
%   See help for "plotSeq" with the following deviations:
%     Seq.plotSequence.wraps is always set to Inf
%     The default value for Seq.plotSequence.figureName is 'Pulse Program
%     (wrapped at each tRep)'.
%     The default value for Seq.plotSequence.Gradients is Seq.plotSeqTR (which
%     defaults to 1:3).
%
% OUTPUT:
%   none
%
% ------------------------------------------------------------------------
% (C) Copyright 2011-2019 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------

%% default input
hParent = 20;
Seq = set_EmptyField(Seq, 'plotSeqTR', 1:3);
% structure with configuration
if ~isfield(Seq, 'plotSequence'), Seq.plotSequence = []; end
Seq.plotSequence = set_EmptyField(Seq.plotSequence, 'hParent', hParent);

Seq.plotSequence.wraps = Inf; % always wrap after each tRep
Seq.plotSequence = set_EmptyField(Seq.plotSequence, 'Gradients', Seq.plotSeqTR);

Seq.plotSequence = set_EmptyField(Seq.plotSequence, 'figureName', 'Pulse Program (wrapped at each tRep)');

%% Call plotSeq
plotSeq(HW, Seq, AQ, TX, Grad);

end
