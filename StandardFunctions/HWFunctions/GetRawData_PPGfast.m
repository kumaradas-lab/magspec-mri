function [Seq, raw_data, num_averages, matrix_check_overflow, timestamps] ...
  = GetRawData_PPGfast(Seq, iDevice)
%% Get raw data from MMRT
%
%       [Seq, raw_data, num_averages, matrix_check_overflow, timestamps] ...
%          = GetRawData_PPGfast(Seq, iDevice)
%
%
% If HW.Function_After_Measurement is set, it is executed just before receiving
% the data. Its value has to be a function handle with the following signature:
%     HW = @(HW)
% If Seq.DataClassSingle is true, the return data will have single floating
% point precision. Otherwise, it will have double floating point precision.
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
%   matrix_check_overflow
%       Boolean matrix in which samples are marked for which the receiver was
%       saturated for at least one of the non-downsampled input. It has the same
%       size as raw_data if acquisitions are downsampled and at least one sample
%       is marked as overflowing. It is empty otherwise.
%
%   timestamps
%       Structure with the following fields containing the recorded timestamps
%       if at least one acquisition window was marked with AQ.GetDigitalInTime.
%       It is empty otherwise.
%     tRep
%         Timestamps of the start of the last three tReps in seconds.
%         Dimensions: [recorded timestamps, tRep]
%     AQ
%         Timestamps of the start of the last three AQ windows in seconds.
%         Dimensions: [recorded timestamps, AQ, tRep]
%     digitalInRising
%         Timestamps of the last three detected rising edges at the digital
%         inputs 1, 2, and 3 in seconds.
%         Dimensions: [recorded timestamps, AQ, tRep, input channel]
%     digitalInFalling
%         Timestamps of the last three detected falling edges at the digital
%         inputs 1, 2, and 3 in seconds.
%         Dimensions: [recorded timestamps, AQ, tRep, input channel]
%       All of these timestamps are recorded at the end of the respective
%       acquisition window.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2015-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function PD.Talker
%#function isemptyfield

