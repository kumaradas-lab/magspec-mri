function cic_corr_factor = CIC_corr(AQ, HW, iDevice)
%% Calculate the correction factors for the cascaded integrator comb filter
%
%   cic_corr_factor = CIC_corr(AQ, HW)
%
% BACKGROUND:
% The output frequency of the ADC (HW.RX.fSample) is mixed down to the input
% frequency of the CIC filter. The local oscillator frequency of this mixer is
% given by the acquisition frequency (AQ.Frequency).
% If the sample frequency of the ADC differs from the sample frequency of the
% acquisition window (AQ.fSample), a CIC filter is used to convert between these
% frequencies. The transfer function of the CIC filter is frequency dependent.
% Its absolute value is given by:
%
%     |H(f)| = (sin(pi*R*M*f))^N / (sin(pi*f))^N
%
% where f is the offset frequency to the acquisition frequency, R is the
% reduction factor, M is the differential delay in the comb stages, and N is the
% number of CIC stages.
%
%
% The function "CIC_corr" calculates the correction function that must be
% multiplied with the image (convolution theorem) to correct for this effect:
%
%     cic_corr_factor = (R*M)^N / |H(f)|
%
% where f is chosen at the centers of the voxels in read direction.
%
% EXAMPLE:
% If data.data contains one line of k-space (encoded in read direction!), the
% correction is applied on the corresponding profile with:
%
%     data.fft1_data = fftshift(ifft(ifftshift(data.data))) .* cic_corr_factor;
%
% -----------------------------------------------------------------------
% (C) Copyright 2011-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% -----------------------------------------------------------------------

end


%#function get_FFTGrid

