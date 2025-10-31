function [Seq, data_S_AQs_TRs] = get_data(Seq, raw_data_AQs, averages, Channel, iDevice, nucleusX)
%% Sort raw data into structure with array of samples x AQ windows x tRep
%
%   [Seq, data_S_AQs_TRs] = get_data(Seq, raw_data_AQs, averages, Channel, iDevice, nucleusX)
%
% INPUT:
%
%   Seq
%         Structure with sequence data. Among others, the following fields are
%         used in this function:
%
%     AQ
%           Structure containing matrices describing the acquisition windows
%           in the pulse program. See "manual openMatlab" for further details.
%
%     HW
%           Structure holding all content of HW as used for the measurement.
%
%     fLarmor
%           Scalar with the Larmor frequency in Hz.
%
%     tRep
%           Vector with the duration of the tReps.
%
%   raw_data_AQs
%         2 dimensional complex numeric array containing the raw data as
%         acquired.
%
%   averages
%         Matrix with the same size as raw_data_AQs containing the number of
%         averages at each sample (or empty matrix).
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
%   data_S_AQs_TRs
%         Structure with arrays containing the sorted results. Amongst others
%         the fields are:
%
%     time_of_tRep
%         Time of the sample acquisition starting in each tRep. The data points
%         are ordered as:
%                 nSamples x nAQs x ntReps
%
%     time_all
%         Same as time_of_tRep, but the cumulative time starting at the
%         beginning of the first tRep.
%
%     data
%         complex data (in T) acquired at each of the above samples.
%
%     averages
%         Matrix with the same size as data_S_AQs_TRs.data containing the number
%         of averages at each sample if there were any missing samples. It is a
%         scalar with the number of averages, otherwise.
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
%     fft1_data
%         If Seq.AQ.CalculateFFTOfData is true, this field contains the
%         FFT of the data in each acquisition window corrected for the CIC
%         filter.
%
%     f_fft1_data
%         If Seq.AQ.CalculateFFTOfData is true, this field contains the
%         frequency corresponding to fft1_data.
%
%     cic_corr
%         If Seq.AQ.CalculateFFTOfData is true, this field contains the
%         correction values for the CIC filter. See also CIC_corr.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2016-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

end


%#function CIC_corr
%#function get_FFTGrid
%#function isemptyfield

