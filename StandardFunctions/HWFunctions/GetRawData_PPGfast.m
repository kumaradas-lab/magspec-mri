function [Seq, raw_data, num_averages] = GetRawData_PPGfast(Seq, iDevice)
%% Get raw data from MMRT
%
%       [Seq, raw_data, num_averages] = GetRawData_PPGfast(Seq)
%
%
% If Seq.IgnoreTimingError is set, no warning is issued on timing errors.
% If HW.RX.WarningNonValidData is set, a warning is issued for missing samples.
%
% OUTPUT:
%
%   raw_data
%       A 2d matrix with the raw data. The samples in each AQ window are
%       oriented along the 1st dimension. The second dimension corresponds to
%       the AQ windows in the pulse program.
%       Invalid data is replaced by HW.RX.NonValidData.
%       If the data is mixed to a center frequency (AQ.Frequency ~=
%       HW.RX.fSample), the function returns a cell array with one or two
%       elements. The first element corresponds to the data from the primary
%       mixer. If applicable, the second element corresponds to the data from
%       the secondary mixer (AQ.FrequencyX, AQ.PhaseX).
%
%   num_averages
%       Integer matrix with the number of actual measurements for each sample.
%       It has the same size as raw_data if there were any missing samples. It
%       is empty otherwise.
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function PD.Talker
%#function isemptyfield

