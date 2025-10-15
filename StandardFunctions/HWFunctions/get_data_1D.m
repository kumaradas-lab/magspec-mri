function [Seq, data_1D] = get_data_1D(Seq, data_S_AQs_TRs, Channel, iDevice, nucleusX)
%% Sort data into structure with all acquisition windows separated by NaN
%
%     [Seq, data_1D] = get_data_1D(Seq, data_S_AQs_TRs, Channel, nucleusX)
%
% INPUT:
%
%   Seq
%         Structure with sequence data. Among others, the following field is
%         used in this function:
%
%     AQ
%           Structure containing matrices describing the acquisition windows
%           in the pulse program. See "manual openMatlab" for further details.
%
%   data_S_AQs_TRs
%         Structure with measurement data as returned by "get_data".
%
%
%   Channel
%         For systems with multiple AQ channels, the number of the channel that
%         is used. (Default: 1)
%
%   iDevice
%         For systems with multiple MRT devices, the number of the device that
%         is used. (Default: 1)
%
%   nucleusX
%         For systems that support mixing acquisitions for two frequencies, this
%         Boolean value indicates if the data corresponds to the X-nucleus.
%         (Default: false)
%
%
% OUTPUT:
%
%   Seq
%         Exact same as input Seq.
%
%   data_1D
%         Structure with arrays containing the results with acquisition windows
%         separated by NaN. Amongst others the fields are:
%
%     time_of_tRep
%         Time of the sample acquisition starting in each tRep. The data points
%         are sorted in columns per tRep.
%
%     time_all
%         Same as time_of_tRep, but the cumulative time starting at the
%         beginning of the first tRep.
%
%     data
%         complex data (in T) acquired at each of the above samples.
%
%     AqFrequency
%         Acquisition frequency at each sample.
%
%     Amplitude2Norm
%         Amplification of the receiver chain in 1/T at each sample.
%
%     device
%         Device number that was used to acquire the data.
%
%     channel
%         Acquisition channel that was used to acquire the data.
%
%     nucleusX
%         Boolean value that indicates that the data was acquired for the
%         X-nucleus on devices that support mixing acquisitions for two
%         frequencies.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2022 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end
