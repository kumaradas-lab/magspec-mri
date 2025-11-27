function [data, dataOut] = get_kSpaceAndImage(data, UseSample, UseAQWindow, UsetRep, MySizeOS, MySize)
%% Calculate oversampled and not oversampled k-Space and Image
%
%   [data, dataOut] = get_kSpaceAndImage(data, AQSlice)
% or (not recommended):
%   [data, dataOut] = get_kSpaceAndImage(data, UseSample, UseAQWindow, UsetRep, MySizeOS, MySize)
%
%
% INPUT:
%
%   data      Structure with measurement data as returned by set_sequence.
%
%   AQSlice   Array of structures with the settings for the slice orientation
%             and setup. If numel(AQSlice) > 1, the images corresponding to the
%             AQSlice structures are returned in an array of data of the same
%             size as AQSlice (un-tested).
%
%             Among others noticeable fields are:
%     UsetRep   Vector or array with the numbers of the tReps that contain the
%               information for the reconstruction of one image. If "UsetRep" is
%               an array, the second dimension corresponds to "independent"
%               images that are stored along the 5-th dimension of the returned
%               k-space and image arrays. The third dimension corresponds to
%               parts of one image (that are combined for corrections). See
%               "partsAverage" below.
%
%     UseAQWindow
%               Scalar, vector or array with the number of the acquisition
%               window in the corresponding tRep that contain the information
%               for the reconstruction of one image. The same rules for
%               multi-dimensional input applies as for "UsetRep". (Default: 1)
%
%     ReadDir
%               Scalar, vector or array. If its size matches AQSlice.UsetRep and
%               AQSlice.UseAQWindow, the corresponding (read encoded) samples
%               are used in reverse order if AQSlice.ReadDir is smaller or equal
%               to 0. (Optional setting with no default)
%
%     partsAverage
%               Scalar integer that indicates the number of images that are
%               created by averaging the k-spaces as defined by
%               partsAverageOrder. The k-space is averaged before reconstructing
%               the images. If 0 (or the size of the 6th dimension of the data),
%               the image parts are kept separately aligned along the 6-th
%               dimension of the returned k-space and image arrays. If 1 (the
%               default), all k-space data is averaged along its 6th dimension.
%
%     partsAverageOrder
%               2d array where each column contains the indices (in the 3rd
%               dimension of UsetRep and/or UseAQWindow, corresponding to the
%               6th dimension in the resulting k-space and image data) that are
%               used for each of the resulting images. By default, the data for
%               the images is equally spaced over the images defined along the
%               3rd dimension of UsetRep (and/or UseAQWindow).
%
%     ZeroFillFactor
%               Scalar or vector with the factors that are used for
%               zero-filling (also referred to as "zero-padding"). The resulting
%               size is the "outputSize" argument of the function "zeroFill".
%               The vector elements correspond to [read, phase(1:3)] directions.
%               (Default: 1, i.e. no zero-padding)
%
%     ZeroFillWindowSize
%               Scalar or vector with the factors that are used for damping high
%               k-space frequencies. The highest k-space frequency in each
%               dimension corresponds to 1. See also the input "winsize.size" of
%               the function "zeroFill". The vector elements correspond to
%               [read, phase(1:3)] directions. A cosine-square filter is used.
%               (Default: Inf, i.e. no k-space filtering)
%
%     ReadOSUsedForImage
%               For CSI or SPI images, number of (oversampled) samples that are
%               used for image reconstruction. If AQSlice.UseForImageStartIdx or
%               AQSlice.UseForImageStopIdx are also set, they take precedence.
%               (Default: AQSlice.ReadOS, i.e. the entire acquisition window.)
%
%     UseForImageStartIdx
%               For CSI or SPI images, index of first (oversampled) sample that
%               is used for image reconstruction. (Default: first sample around
%               the center according to ReadOSUsedForImage taking Spin Echo
%               encoding into account)
%
%     UseForImageStopIdx
%               For CSI or SPI images, index of last (oversampled) sample that
%               is used for image reconstruction. (Default: last sample around
%               the center according to ReadOSUsedForImage taking Spin Echo
%               encoding into account)
%
%     kLineStart, kLineEnd
%               For non-uniform distribution of k-space data (ZTE), these are
%               the k-space coordinates in 1/meter with the start and end of the
%               acquired rays. (Mandatory for non-uniformly distributed k-space
%               data.)
%
%     NUFFTReconstruction
%               If false, a basic interpolation algorithm is used to re-grid the
%               non-uniformly acquired k-space to a regular grid before
%               reconstructing the image with a Fast Fourier Transform (FFT). If
%               true, a non-uniform discrete Fourier transform is used instead.
%               This only applies to k-space data that is non-uniformly
%               distributed in k-space. For uniformly distributed k-space data,
%               a FFT is always used. (Default: true)
%
%
% OUTPUT:
%   data        Same as input "data" with (among others) the following
%               additional fields:
%               The k-space and image arrays follow this convention:
%               1st dimension: read direction
%               2nd dimension: phase 1 direction
%               3rd dimension: phase 2 direction
%               4th dimension: phase 3 direction
%               5th dimension: multiple (independent) images within AQSlice
%               6th dimension: parts of these images that are combined for
%                              corrections (optional, see
%                              "AQSlice.partsAverage")
%
%     kSpaceOsRaw Over-sampled k-space as acquired.
%     ImageOsRaw  Image(s) reconstructed from that k-space.
%     kSpaceOs    Over-sampled k-space(s) with zero-fill window size correction
%                 (also applied to all other k-spaces below).
%     ImageOs     Image(s) reconstructed from that k-space(s) including
%                 CIC-correction.
%     kSpaceOsZ   Over-sampled k-space(s) with zero filling.
%     ImageOsZ    Image(s) reconstructed from that k-space(s) including
%                 CIC-correction.
%     ImageZ      Center part of ImageOsZ corresponding to zero-filling.
%     kSpaceZ     k-space(s) corresponding to that image(s).
%     Image       Center part of ImageOs.
%     kSpace      k-space(s) corresponding to that image(s).
%
%   dataOut     Same as "data" but only the minimum number of fields to display
%               the k-space and image. This is useful to reduce computer memory
%               usage.
%
% EXAMPLES:
% function call e.g.
%  [data, dataOut] = get_kSpaceAndImage(data, SeqOut.AQSlice(1))
% or function call (not recommended) e.g.
%  [data, dataOut] = get_kSpaceAndImage(data, 1:SeqOut.AQ.nSamples(1), 1,  1:length(SeqOut.tRep),...
%                                          [SeqOut.AQSlice.nRead*SeqOut.AQSlice.ReadOS,  SeqOut.AQSlice.nPhase*SeqOut.AQSlice.PhaseOS, SeqOut.AQSlice.nPhase3D*SeqOut.AQSlice.PhaseOS3D],...
%                                          [SeqOut.AQSlice.nRead,                        SeqOut.AQSlice.nPhase,                        SeqOut.AQSlice.nPhase3D])
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2025 Pure Devices GmbH, Wuerzburg, Germany
%     www.pure-devices.com
% ------------------------------------------------------------------------------

if nargin == 2
  if isstruct(UseSample)
    AQSlice = UseSample;
  else
    error('second argument is not a struct');
  end
else
  AQSlice = struct();
end

if isemptyfield(AQSlice(1), 'partsAverage')
  % Calculate average over images defined along 3rd dimension of AQSlice.UsetRep
  AQSlice(1).partsAverage = 1;
end

for iAQSlice = 1:numel(AQSlice)
  if nargin == 2
    UsetRep = AQSlice(iAQSlice).UsetRep;
    UseAQWindow = AQSlice(iAQSlice).UseAQWindow;
  end
  % if isvector(UsetRep), UsetRep = UsetRep(:); end
  % if isvector(UseAQWindow), UseAQWindow = UseAQWindow(:); end

  [nKLines, nImages, nParts] = size(UsetRep);
  if isscalar(UseAQWindow)
    UseAQWindow(1:nKLines,1:nImages,1:nParts) = UseAQWindow;
  end
  if AQSlice(1).partsAverage>nParts || ~AQSlice(1).partsAverage
    AQSlice(1).partsAverage = nParts;
  end
  if isemptyfield(AQSlice(1), 'partsAverageOrder')
    nA = nParts/AQSlice(1).partsAverage;
    AQSlice(1).partsAverageOrder = zeros(nA, AQSlice(1).partsAverage);
    for iAverage = 1:AQSlice(1).partsAverage
      % Select images (along 6th dimension of data) that are used for each average
      AQSlice(1).partsAverageOrder(1:nA,iAverage) = iAverage + (0:AQSlice(1).partsAverage:nParts-1);
    end
  elseif size(AQSlice(1).partsAverageOrder, 2) ~= AQSlice(1).partsAverage
    error('PD:get_kSpaceAndImage:partsAverageOrderSize', ...
      'Size of AQSlice.partsAverageOrder does not match AQSlice.partsAverage.');
  end

  for iImage = 1:nImages
    for iPart = 1:nParts
      if nargin == 2
        if isemptyfield(AQSlice(iAQSlice), 'AQWindowDirection'), AQSlice(iAQSlice).AQWindowDirection = 1; end
        if isemptyfield(AQSlice(iAQSlice), 'ReadRadial'), AQSlice(iAQSlice).ReadRadial = 0; end

        if AQSlice(iAQSlice).ReadRadial
          MySizeOS = [AQSlice(iAQSlice).nRead*AQSlice(iAQSlice).ReadOS,  numel(AQSlice(iAQSlice).UsetRep)];
          MySize = MySizeOS;
        else
          MySizeOS = [ AQSlice(iAQSlice).nRead*AQSlice(iAQSlice).ReadOS,    AQSlice(iAQSlice).nPhase(1)*AQSlice(iAQSlice).PhaseOS(1), AQSlice(iAQSlice).nPhase(2)*AQSlice(iAQSlice).PhaseOS(2), AQSlice(iAQSlice).nPhase(3)*AQSlice(iAQSlice).PhaseOS(3)];
          MySize   = [ AQSlice(iAQSlice).nRead,                             AQSlice(iAQSlice).nPhase(1),                              AQSlice(iAQSlice).nPhase(2),                              AQSlice(iAQSlice).nPhase(3)                             ];
        end

        if isemptyfield(AQSlice(iAQSlice), 'kLineOrder')
          AQSlice(iAQSlice).kLineOrder = 1:prod(MySizeOS(2:end));
        end

        UseSample = 1:AQSlice(iAQSlice).nRead*AQSlice(iAQSlice).ReadOS;
        [szData(1),szData(2),szData(3),szData(4),szData(5)] = size(data(iAQSlice).data);

        % t1=(1:max(UseSample)).';
        % t2=(max(UseSample):-1:1).';
        % t3=zeros(max(UseSample),numel(AQSlice(iAQSlice).AQWindowDirection));
        % for t=1:numel(AQSlice(iAQSlice).AQWindowDirection)
        %   if AQSlice(iAQSlice).AQWindowDirection(t)==1
        %     t3(:,t)=t1;
        %   else
        %     t3(:,t)=t2;
        %   end
        % end
        t3 = repmat((1:max(UseSample)).', 1, numel(AQSlice(iAQSlice).AQWindowDirection));
        t3(:,AQSlice(iAQSlice).AQWindowDirection(:)~=1) = flipud(t3(:,AQSlice(iAQSlice).AQWindowDirection(:)~=1));

        AQSlice(iAQSlice).AQWindowSampleOrder = reshape(repmat(t3, [1,prod(MySizeOS(2:end))/numel(AQSlice(iAQSlice).AQWindowDirection)]), MySizeOS);
        AQSlice(iAQSlice).AQWindowSampleOrder(:) = AQSlice(iAQSlice).AQWindowSampleOrder(:) + reshape(max(UseSample)*ones(max(UseSample),1)*(0:prod(MySizeOS(2:end))-1), [], 1);
        if isemptyfield(AQSlice(iAQSlice), 'ReadOSUsedForImage'), AQSlice(iAQSlice).ReadOSUsedForImage = AQSlice(iAQSlice).ReadOS; end
        evenNumberOfSamples = ~mod(AQSlice(iAQSlice).nRead*AQSlice(iAQSlice).ReadOS, 2);
        if isemptyfield(AQSlice(iAQSlice), 'UseForImageStartIdx')
          AQSlice(iAQSlice).UseForImageStartIdx = ...
            round(numel(UseSample)/2) + evenNumberOfSamples ...
            + ceil(-AQSlice(iAQSlice).ReadOSUsedForImage/2);
        end
        if AQSlice(iAQSlice).UseForImageStartIdx < 1 ...
            || AQSlice(iAQSlice).UseForImageStartIdx > AQSlice(iAQSlice).nRead*AQSlice(iAQSlice).ReadOS
          error('PD:get_kSpaceAndImage:UseForImageStartIdxOutOfRange', ...
            'UseForImageStartIdx (%d) is out of range (1 - %d).', ...
            AQSlice(iAQSlice).UseForImageStartIdx, ...
            AQSlice(iAQSlice).nRead*AQSlice(iAQSlice).ReadOS);
        end
        if isemptyfield(AQSlice(iAQSlice), 'UseForImageStopIdx')
          AQSlice(iAQSlice).UseForImageStopIdx = ...
            round(numel(UseSample)/2) + evenNumberOfSamples ...
            + ceil(AQSlice(iAQSlice).ReadOSUsedForImage/2-1);
        end
        if AQSlice(iAQSlice).UseForImageStopIdx < 1 ...
            || AQSlice(iAQSlice).UseForImageStopIdx > AQSlice(iAQSlice).nRead*AQSlice(iAQSlice).ReadOS
          error('PD:get_kSpaceAndImage:UseForImageStopIdxOutOfRange', ...
            'UseForImageStopIdx (%d) is out of range (1 - %d).', ...
            AQSlice(iAQSlice).UseForImageStopIdx, ...
            AQSlice(iAQSlice).nRead*AQSlice(iAQSlice).ReadOS);
        end
        if isemptyfield(AQSlice(iAQSlice), 'UseAQWindow'), AQSlice(iAQSlice).UseAQWindow = 1; end
        if isemptyfield(AQSlice(iAQSlice), 'UsetRep'), AQSlice(iAQSlice).UsetRep = 1:size(data.data, 3); end
      % elseif and(AQSlice.nRead == 1, AQSlice.nPhase(1)>1)
      %   MySizeOS=[AQSlice.nPhase(1)*AQSlice.PhaseOS(1),  AQSlice.nPhase(2)*AQSlice.PhaseOS(2), AQSlice.nPhase(3)*AQSlice.PhaseOS(3)];
      %   MySize=[AQSlice.nPhase(1),                        AQSlice.nPhase(2),                        AQSlice.nPhase(3)];
      %   if AQSlice.ReadOS>16
      %     data.kSpaceCsiRaw=reshape(data.data(UseSample,UseAQWindow,UsetRep)...
      %                          ,[AQSlice.ReadOS,MySizeOS]);
      %     data.kSpaceCsiRaw=data.kSpaceCsiRaw(AQSlice.AQWindowSampleOrder);
      %     data.kSpaceCsiRaw=reshape(data.kSpaceCsiRaw,[AQSlice.ReadOS,MySizeOS]);
      %     data.ImageCsiRaw=fftshift(ifftn(ifftshift(data.kSpaceCsiRaw)));
      %     data.ImageCsi=data.ImageCsiRaw.*reshape(data.cic_corr(UseSample,UseAQWindow,UsetRep),[AQSlice.ReadOS,MySizeOS]);
      %     data.ImageCsiFrequency=reshape(data.f_fft1_data(UseSample,UseAQWindow,UsetRep),[AQSlice.ReadOS,MySizeOS]);
      %     if isemptyfield(AQSlice, 'ReadZeroFill')
      %       AQSlice.ReadZeroFill = 8;
      %     end
      %
      %       data.ImageCsiRawZero=fftshift(ifftn(ifftshift(zeroFill(data.kSpaceCsiRaw,[AQSlice.ReadOS*AQSlice.ReadZeroFill,MySizeOS]))));
      %       % data.ImageCsi=data.ImageCsiRawZero.*reshape(data.cic_corr(UseSample,UseAQWindow,UsetRep),[AQSlice.ReadOS,MySizeOS]);
      %       % data.ImageCsiFrequency=reshape(data.f_fft1_data(UseSample,UseAQWindow,UsetRep),[AQSlice.ReadOS,MySizeOS]);
      %       BW=data.f_fft1_data(UseSample(end),UseAQWindow,UsetRep(1))-data.f_fft1_data(1,UseAQWindow,UsetRep(1));
      %       if mod(AQSlice.nRead*AQSlice.ReadOS,2)
      %         data.ImageCsiFrequencyZero=data.f_fft1_data(ceil(numel(UseSample)/2),UseAQWindow,UsetRep(1))+BW/2*linspace(-1+1/(AQSlice.ReadZeroFill*AQSlice.ReadOS),1-1/(AQSlice.ReadZeroFill*AQSlice.ReadOS),AQSlice.ReadZeroFill*AQSlice.ReadOS).';
      %       else
      %         data.ImageCsiFrequencyZero=data.f_fft1_data(round(numel(UseSample)/2)+1,UseAQWindow,UsetRep(1))+BW/2*linspace(-1,1-2/(AQSlice.ReadZeroFill*AQSlice.ReadOS),AQSlice.ReadZeroFill*AQSlice.ReadOS).';
      %       end
      %       data.ImageCsiFrequencyZero=reshape(data.ImageCsiFrequencyZero*ones(1,prod(MySizeOS)),[AQSlice.ReadOS*AQSlice.ReadZeroFill,MySizeOS]);
      %     end
      %   end
      end
      if isemptyfield(AQSlice(iAQSlice), 'ZeroFillFactor')
        AQSlice(iAQSlice).ZeroFillFactor = 1;
      end
      if isemptyfield(AQSlice(iAQSlice), 'ZeroFillWindowSize')
        AQSlice(iAQSlice).ZeroFillWindowSize = Inf;
      end

      if isemptyfield(AQSlice(iAQSlice), 'kLineEnd')
        isNUFFT = false;
      else
        isNUFFT = true;
      end

      if isNUFFT
        if isscalar(AQSlice(iAQSlice).ZeroFillFactor)
          AQSlice(iAQSlice).ZeroFillFactor = ...
            [1, AQSlice(iAQSlice).ZeroFillFactor([1 1 1])];
        end
        if isscalar(AQSlice(iAQSlice).ZeroFillWindowSize)
          AQSlice(iAQSlice).ZeroFillWindowSize = ...
            [Inf, AQSlice(iAQSlice).ZeroFillWindowSize([1 1 1])];
        end
        AQSlice(iAQSlice).ZeroFillFactor(MySize<=1) = 1;
        AQSlice(iAQSlice).ZeroFillWindowSize(MySize<=1) = Inf;

        if isemptyfield(AQSlice(iAQSlice), 'NUFFTReconstruction')
          % image reconstruction with non-uniform Fourier transform
          % otherwise, basic re-gridding with scatteredInterp
          AQSlice(iAQSlice).NUFFTReconstruction = true;
        end

        usedkLine = UseAQWindow(:,iImage,iPart)~=0 & UsetRep(:,iImage,iPart)~=0;
        UseAQWindow_tRep = sub2ind(szData(2:3), UseAQWindow(usedkLine,iImage,iPart), UsetRep(usedkLine,iImage,iPart));

        % number of acquired rays in k-space
        numKRays = MySizeOS(1);

        % remove size in read direction for image reconstruction
        MySizeOSZ = ceil(MySizeOS.*(AQSlice(iAQSlice).ZeroFillFactor));
        MySizeZ = ceil(MySize.*(AQSlice(iAQSlice).ZeroFillFactor));
        MySizeOSZ(1) = 1;
        MySizeZ(1) = 1;
        MySizeOS(1) = 1;
        MySize(1) = 1;
        data(iAQSlice).ZeroFillFactorOS = MySizeOSZ./MySizeOS;
        data(iAQSlice).ZeroFillFactor = MySizeZ./MySize;

        if isemptyfield(AQSlice(iAQSlice), 'kLineStart')
          error('PD:get_kSpaceAndImage:NokLineStart', ...
            ['k-Space coordinates for the start and end points of the acquired rays ', ...
            '(AQSlice.kLineStart and AQSlice.kLineEnd) are needed for ', ...
            'non-uniform discrete Fourier transform.']);
        end
        % generate k-space coordinate for each sample in m^(-1) (or in s)
        k1 = cell2mat(arrayfun(@(s,e) linspace(s, e, numKRays).', ...
          AQSlice(iAQSlice).kLineStart(1,:), AQSlice(iAQSlice).kLineEnd(1,:), 'UniformOutput', false));
        k2 = cell2mat(arrayfun(@(s,e) linspace(s, e, numKRays).', ...
          AQSlice(iAQSlice).kLineStart(2,:), AQSlice(iAQSlice).kLineEnd(2,:), 'UniformOutput', false));
        k3 = cell2mat(arrayfun(@(s,e) linspace(s, e, numKRays).', ...
          AQSlice(iAQSlice).kLineStart(3,:), AQSlice(iAQSlice).kLineEnd(3,:), 'UniformOutput', false));

        usedDim = AQSlice(iAQSlice).nPhase>1;

        % resolution for scaling kx, ky, kz to (-pi, pi]
        resolution = AQSlice(iAQSlice).sizePhase ./ MySize(2:4);
        % or end times (in s) for non-encoded directions
        nonEncodedDirs = isinf(AQSlice(iAQSlice).sizePhase) & (AQSlice(iAQSlice).nPhase > 1);
        % FIXME: Is it save to assume that all acquisition windows have the
        % same duration? (It probably is...)
        resolution(nonEncodedDirs) = 1./(2*sqrt(sum(AQSlice(iAQSlice).kLineEnd(nonEncodedDirs,1).^2)));

        resolution(~isfinite(resolution)) = 1;

        if AQSlice(iAQSlice).NUFFTReconstruction
          % non-uniform Fourier transform

          if iImage == 1 && iPart == 1
            % FIXME: Do we need to reconstruct the other types of k-spaces?
            % data(iAQSlice).kSpaceOsRaw = zeros([MySizeOS nImages nParts], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageOsRaw = zeros([MySizeOS nImages nParts], 'like', data(iAQSlice).data);
            % data(iAQSlice).kSpaceOsZ = zeros([MySizeOSZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageOsZ = zeros([MySizeOSZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpaceOs = zeros([MySizeOS nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageOs = zeros([MySizeOS nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpaceZ = zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageZ = zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpace = zeros([MySize nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).Image = zeros([MySize nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
          end

          % normalize end points to ellipsoid in k-space
          if sum(usedDim) > 2
            % 3d
            k_abs_ellipsoid = sqrt((k1/max(abs(k1(:)))).^2 + ...
              (k2/max(abs(k2(:)))).^2 + ...
              (k3/max(abs(k3(:)))).^2);
          elseif sum(usedDim) > 1
            % 2d
            k_abs_ellipsoid = sqrt((k1/max(abs(k1(:)))).^2 + ...
              (k2/max(abs(k2(:)))).^2);
          else
            % FIXME: This hasn't been checked if it is working correctly.
            k_abs_ellipsoid = abs(k1/max(abs(k1(:))));
          end
          k_norm_ellipsoid = k_abs_ellipsoid / max(k_abs_ellipsoid(:));

          % correct for point density with distance to k-space center
          % Rays are closer to each other at the k-space center than in the
          % outer regions.
          weight_density = k_norm_ellipsoid .^ (sum(usedDim)-1);
          % increase weight to compensate for missing center of k-space
          weight_density(1,:) = 1/diff(k_norm_ellipsoid(1:2,1)) ./ (sum(usedDim)) .* ...
            (k_norm_ellipsoid(1)) .^ (sum(usedDim));

          if isemptyfield(AQSlice(iAQSlice), 'ReadDensityCorrection')
            AQSlice(iAQSlice).ReadDensityCorrection = ones(1, numKRays);
          end

          % FIXME: Double-check if the calculation of these weighting factors
          %        is correct!

          % FIXME: Do we need to normalize these weights to get meaningful image
          % amplitudes?
          weight = bsxfun(@mtimes, (weight_density/mean(weight_density(:))) .* ...
            cosd(k_norm_ellipsoid*90).^2, ...  % cos^2 filter; FIXME: Do we need this to be adjustable?
            AQSlice(iAQSlice).ReadDensityCorrection);  % non-uniform density of readout rays

          if 0
            % FIXME: Should there be a setting for enabling these plots?
            % graphical representation of the used weigths
            hf = figure(5501); clf(hf);
            hax = axes(hf);
            % pcolor(hax, k1, k2, weight);
            surf(hax, k1, k2, k3, weight, 'FaceAlpha', 0.5);
            shading(hax, 'interp');
            title(hax, 'total weight');
            hf = figure(5502); clf(hf);
            hax = axes(hf);
            plot(hax, AQSlice(iAQSlice).ReadDensityCorrection)
            shading(hax, 'interp');
            title(hax, 'sequence specific weights');
            hf = figure(5503); clf(hf);
            hax = axes(hf);
            plot(hax, weight_density)
            shading(hax, 'interp');
            title(hax, 'ray density (surface effect)');
            hf = figure(5504); clf(hf);
            hax = axes(hf);
            plot(hax, cosd(k_norm_ellipsoid*90).^2)
            shading(hax, 'interp');
            title(hax, 'cos^2 filter');
           end

          kSpaceRaw = data(iAQSlice).data(UseSample,UseAQWindow_tRep);
          if isa(kSpaceRaw, 'single')
            castFcn = @single;
            % FINUFFT requires double type
            kSpaceRaw = double(kSpaceRaw);
          else
            castFcn = @double;
          end

          numeric_stab_fac = max(abs(kSpaceRaw(:)));  % scale to 1 for numeric stability

          % FIXME: Can we always be sure that the dimensions are filled from the
          %        start?
          if sum(usedDim) > 2
            % 3d image
            data(iAQSlice).ImageOsRaw(:,:,:,:,iImage,iPart) = reshape(castFcn(finufft3d1(...
              k1(:) * resolution(1) * 2 * pi, ...
              k2(:) * resolution(2) * 2 * pi, ...
              k3(:) * resolution(3) * 2 * pi, ...
              kSpaceRaw(:) .* weight(:) / numeric_stab_fac, ...
              1, 1e-13, ...
              MySizeOS(2), MySizeOS(3), MySizeOS(4))), MySizeOS) * numeric_stab_fac / numel(k1);

            % FIXME: CIC correction?
            data(iAQSlice).ImageOs(:,:,:,:,iImage,iPart) = ...
              data(iAQSlice).ImageOsRaw(:,:,:,:,iImage,iPart) / AQSlice.VoxelVolume * AQSlice.AreaCoil;

            data(iAQSlice).ImageOsZ(:,:,:,:,iImage,iPart) = reshape(castFcn(finufft3d1(...
              k1(:) * resolution(1) * 2 * pi / data(iAQSlice).ZeroFillFactorOS(2), ...
              k2(:) * resolution(2) * 2 * pi / data(iAQSlice).ZeroFillFactorOS(3), ...
              k3(:) * resolution(3) * 2 * pi / data(iAQSlice).ZeroFillFactorOS(4), ...
              kSpaceRaw(:) .* weight(:) ./ numeric_stab_fac, ...
              1, 1e-13, ...
              MySizeOSZ(2), MySizeOSZ(3), MySizeOSZ(4))), MySizeOSZ) .* numeric_stab_fac / numel(k1) ...
              / AQSlice(iAQSlice).VoxelVolume * prod(data(iAQSlice).ZeroFillFactorOS) * AQSlice(iAQSlice).AreaCoil;
          elseif sum(usedDim) > 1
            % 2d image
            data(iAQSlice).ImageOsRaw(:,:,:,:,iImage,iPart) = reshape(castFcn(finufft2d1(...
              k1(:) * resolution(1) * 2 * pi, ...
              k2(:) * resolution(2) * 2 * pi, ...
              kSpaceRaw(:) .* weight(:) / numeric_stab_fac, ...
              1, 1e-13, ...
              MySizeOS(2), MySizeOS(3))), MySizeOS) * numeric_stab_fac / numel(k1);

            % FIXME: CIC correction?
            data(iAQSlice).ImageOs(:,:,:,:,iImage,iPart) = ...
              data(iAQSlice).ImageOsRaw(:,:,:,:,iImage,iPart) / AQSlice.VoxelVolume * AQSlice.AreaCoil;

            data(iAQSlice).ImageOsZ(:,:,:,:,iImage,iPart) = reshape(castFcn(finufft2d1(...
              k1(:) * resolution(1) * 2 * pi / data(iAQSlice).ZeroFillFactorOS(2), ...
              k2(:) * resolution(2) * 2 * pi / data(iAQSlice).ZeroFillFactorOS(3), ...
              kSpaceRaw(:) .* weight(:) ./ numeric_stab_fac, ...
              1, 1e-13, ...
              MySizeOSZ(2), MySizeOSZ(3))), MySizeOSZ) .* numeric_stab_fac / numel(k1) ...
              / AQSlice(iAQSlice).VoxelVolume * prod(data(iAQSlice).ZeroFillFactorOS) * AQSlice(iAQSlice).AreaCoil;
          else
            % 1d image
            % FIXME: This hasn't been checked if it is working correctly.
            data(iAQSlice).ImageOsRaw(:,:,:,:,iImage,iPart) = reshape(castFcn(finufft1d1(...
              k1(:) * resolution(1) * 2 * pi, ...
              kSpaceRaw(:) .* weight(:) / numeric_stab_fac, ...
              1, 1e-13, ...
              MySizeOS(2))), MySizeOS) * numeric_stab_fac;

            % FIXME: CIC correction?
            data(iAQSlice).ImageOs(:,:,:,:,iImage,iPart) = ...
              data(iAQSlice).ImageOsRaw(:,:,:,:,iImage,iPart) / AQSlice.VoxelVolume * AQSlice.AreaCoil;

            data(iAQSlice).ImageOsZ(:,:,:,:,iImage,iPart) = reshape(castFcn(finufft2d1(...
              k1(:) * resolution(1) * 2 * pi / data(iAQSlice).ZeroFillFactorOS(2), ...
              kSpaceRaw(:) .* weight(:) ./ numeric_stab_fac, ...
              1, 1e-13, ...
              MySizeOSZ(2))), MySizeOSZ) .* numeric_stab_fac ...
              / AQSlice(iAQSlice).VoxelVolume * prod(data(iAQSlice).ZeroFillFactorOS) * AQSlice(iAQSlice).AreaCoil;
          end
        else
          % Interpolate k-space data to cartesian grid.
          % This is pretty slow and consumes a lot of memory for 3d images.

          if AQSlice(iAQSlice).ZeroFillWindowSize > 1
            warning('PD:get_kSpaceAndImage:ZeroFillWindowSizeTooHigh', ...
              'A zero-fill window size larger than 1 might not make sense for ZTE.')
            % AQSlice(iAQSlice).ZeroFillWindowSize = 1;
          end

          % This messes with the order of fields in the "data" structure.
          % Assign to temporary structure instead of to original structure (at the
          % cost of duplicating memory).
          dataTmp = get_kSpaceAndImageTicks(data(iAQSlice), AQSlice(iAQSlice));

          % Normalize data. The interpolation seems to be unstable otherwise.
          numeric_factor = max(abs(reshape(data(iAQSlice).data(UseSample,UseAQWindow_tRep), 1, [])));

          if iImage == 1 && iPart == 1
            data(iAQSlice).kSpaceOsRaw = zeros([MySizeOS nImages nParts], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageOsRaw = zeros([MySizeOS nImages nParts], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpaceOsZ = zeros([MySizeOSZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageOsZ = zeros([MySizeOSZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpaceOs = zeros([MySizeOS nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageOs = zeros([MySizeOS nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpaceZ = zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageZ = zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpace = zeros([MySize nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).Image = zeros([MySize nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
          end

          % For non-encoded directions, k1,k2,k3 are in seconds.
          % Convert to "phase samples".
          if (nonEncodedDirs(1))
            k1 = k1 * resolution(1) * MySizeOS(2);
          end
          if (nonEncodedDirs(2))
            k2 = k2 * resolution(2) * MySizeOS(3);
          end
          if (nonEncodedDirs(2))
            k3 = k3 * resolution(3) * MySizeOS(4);
          end
          if sum(usedDim) > 2
            % 3d
            kSpaceInterp = scatteredInterpolant(k1(:), k2(:), k3(:), ...
              reshape(data(iAQSlice).data(UseSample,UseAQWindow_tRep) / numeric_factor, [], 1), 'natural', 'none');
            [k1_grid, k2_grid, k3_grid] = ndgrid(dataTmp.kTicks(1).PhaseOs, ...
              dataTmp.kTicks(2).PhaseOs, ...
              dataTmp.kTicks(3).PhaseOs);
            clear dataTmp
            data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart) = ...
              reshape(kSpaceInterp(k1_grid, k2_grid, k3_grid) * numeric_factor, MySizeOS);
            clear k1_grid k2_grid k3_grid kSpaceInterp
          elseif sum(usedDim) > 1
            % 2d
            kSpaceInterp = scatteredInterpolant(k1(:), k2(:), ...
              reshape(data(iAQSlice).data(UseSample,UseAQWindow_tRep) / numeric_factor, [], 1), 'natural', 'none');
            [k1_grid, k2_grid] = ndgrid(dataTmp.kTicks(1).PhaseOs, ...
              dataTmp.kTicks(2).PhaseOs);
            clear dataTmp
            data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart) = ...
              reshape(kSpaceInterp(k1_grid, k2_grid) * numeric_factor, MySizeOS);
            clear k1_grid k2_grid kSpaceInterp
          else
            % FIXME: Does 1d even make sense?
            data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart) = ...
              reshape(interp1(k1(:), ...
                              reshape(data(iAQSlice).data(UseSample,UseAQWindow_tRep) / numeric_factor, [], 1), ...
                              dataTmp.kTicks(1).PhaseOs, 'natural', 'none') * numeric_factor, ...
                      MySizeOS);
            clear dataTmp
          end
          % set extrapolated values to flat 0
          data(iAQSlice).kSpaceOsRaw(isnan(data(iAQSlice).kSpaceOsRaw)) = 0;

          data(iAQSlice).kSpaceOs(:,:,:,:,iImage,iPart) = ...
            zeroFill(data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart), MySizeOS, AQSlice(iAQSlice).ZeroFillWindowSize) ...
            / AQSlice(iAQSlice).VoxelVolume * AQSlice(iAQSlice).AreaCoil;

          % Fourier transform the re-gridded k-space
          data(iAQSlice).ImageOs(:,:,:,:,iImage,iPart) = ...
            fftshift(ifftn(ifftshift(data(iAQSlice).kSpaceOs(:,:,:,:,iImage,iPart))));
          data(iAQSlice).kSpaceOsZ(:,:,:,:,iImage,iPart) = ...
            zeroFill(data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart), MySizeOSZ, AQSlice(iAQSlice).ZeroFillWindowSize) ...
            / AQSlice.VoxelVolume * prod(data(iAQSlice).ZeroFillFactorOS) * AQSlice(iAQSlice).AreaCoil;
          data(iAQSlice).ImageOsZ(:,:,:,:,iImage,iPart) = ...
            fftshift(ifftn(ifftshift(data(iAQSlice).kSpaceOsZ(:,:,:,:,iImage,iPart))));
        end

        data = cutImageOs(data(iAQSlice), AQSlice(iAQSlice), MySizeOSZ, MySizeZ, MySizeOS, MySize, iImage, iPart);

        if AQSlice(iAQSlice).NUFFTReconstruction
          % FIXME: Do we want to (optionally) return the k-spaces that
          %        correspond to the reconstructed images also with basic
          %        regridding?
          data.kSpaceOs(:,:,:,:,iImage,iPart) = fftshift(fftn(ifftshift(data.ImageOs(:,:,:,:,iImage,iPart))));
        end

      else
        % RAW data of the oversampled k-Space
        if AQSlice(iAQSlice).nRead == 1 && prod(AQSlice(iAQSlice).nPhase) > 1
          %% only phase encoding (SPI/CSI)
          if isscalar(AQSlice(iAQSlice).ZeroFillFactor)
            AQSlice(iAQSlice).ZeroFillFactor = ...
              [1, AQSlice(iAQSlice).ZeroFillFactor([1 1 1])];
          end
          if isscalar(AQSlice(iAQSlice).ZeroFillWindowSize)
            AQSlice(iAQSlice).ZeroFillWindowSize = ...
              [Inf, AQSlice(iAQSlice).ZeroFillWindowSize([1 1 1])];
          end
          AQSlice(iAQSlice).ZeroFillFactor(end+1:4) = 1;
          AQSlice(iAQSlice).ZeroFillWindowSize(end+1:4) = Inf;
          MySizeCsiOS = MySizeOS;
          MySizeOS(1) = 1;                    % for mean Image
          MySizeCsiOSZ = ceil(MySizeCsiOS.*(AQSlice(iAQSlice).ZeroFillFactor));
          MySizeCsiOSZ(MySizeCsiOSZ<=1) = 1;
          MySizeOSZ = [1 MySizeCsiOSZ(2:4)];  % for mean Image
          MySizeZ = ceil(MySize.*[1 AQSlice(iAQSlice).ZeroFillFactor(2:4)]);
          MySizeZ(MySize<=1) = 1;
          data(iAQSlice).ZeroFillFactorCsiOS = MySizeCsiOSZ./MySizeCsiOS;
          data(iAQSlice).ZeroFillFactorOS = MySizeOSZ./MySizeOS;
          data(iAQSlice).ZeroFillFactor = MySizeZ./MySize;
          UseAQWindow_tRep = sub2ind(szData(2:3), UseAQWindow(:,iImage,iPart), UsetRep(:,iImage,iPart));

          if AQSlice(iAQSlice).ReadOS > ...
              (AQSlice(iAQSlice).UseForImageStopIdx - AQSlice(iAQSlice).UseForImageStartIdx + 1)
            if iImage == 1 && iPart == 1
              data(iAQSlice).kSpaceCsiRaw = zeros([MySizeCsiOS nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).ImageCsiRaw = zeros([MySizeCsiOS nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).ImageCsi = zeros([MySizeCsiOS nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).ImageCsiFrequency = zeros([MySizeCsiOS nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).ImageCsiRawZero = zeros([MySizeCsiOSZ nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).ImageCsiFrequencyZero = zeros([MySizeCsiOSZ(1) 1 1 1 nImages nParts], 'like', data(iAQSlice).data);
            end
            kSpaceCsiRaw = reshape(data(iAQSlice).data(UseSample,UseAQWindow_tRep), MySizeCsiOS);
            kSpaceCsiRaw = kSpaceCsiRaw(AQSlice(iAQSlice).AQWindowSampleOrder);
            data(iAQSlice).kSpaceCsiRaw(:,:,:,:,iImage,iPart) = reshape(kSpaceCsiRaw, MySizeCsiOS);
            data(iAQSlice).ImageCsiRaw(:,:,:,:,iImage,iPart) = fftshift(ifftn(ifftshift(data(iAQSlice).kSpaceCsiRaw(:,:,:,:,iImage,iPart))));
            data(iAQSlice).ImageCsi(:,:,:,:,iImage,iPart) = data(iAQSlice).ImageCsiRaw(:,:,:,:,iImage,iPart) .* reshape(data(iAQSlice).cic_corr(UseSample,UseAQWindow_tRep), MySizeCsiOS);
            data(iAQSlice).ImageCsiFrequency(:,:,:,:,iImage,iPart) = reshape(data(iAQSlice).f_fft1_data(UseSample,UseAQWindow_tRep), MySizeCsiOS);

            data(iAQSlice).ImageCsiRawZero(:,:,:,:,iImage,iPart) = fftshift(ifftn(ifftshift(zeroFill(data(iAQSlice).kSpaceCsiRaw(:,:,:,:,iImage,iPart), MySizeCsiOSZ))));
            % data.ImageCsi=data.ImageCsiRawZero.*reshape(data(iAQSlice).cic_corr(UseSample,UseAQWindow,UsetRep),[AQSlice.ReadOS,MySizeOS]);
            % data.ImageCsiFrequency=reshape(data(iAQSlice).f_fft1_data(UseSample,UseAQWindow,UsetRep),[AQSlice.ReadOS,MySizeOS]);
            BW = (data(iAQSlice).f_fft1_data(UseSample(2),UseAQWindow_tRep(1)) - data(iAQSlice).f_fft1_data(UseSample(1),UseAQWindow_tRep(1))) .* length(UseSample);
            if mod(AQSlice(iAQSlice).nRead*AQSlice(iAQSlice).ReadOS, 2)
              ImageCsiFrequencyZero = data(iAQSlice).f_fft1_data(ceil(numel(UseSample)/2),UseAQWindow_tRep(1)) + ...
                BW/2*linspace(-1+1/MySizeCsiOSZ(1), 1-1/MySizeCsiOSZ(1), MySizeCsiOSZ(1)).';
            else
              ImageCsiFrequencyZero = data(iAQSlice).f_fft1_data(round(numel(UseSample)/2)+1,UseAQWindow_tRep(1)) + ...
                BW/2*linspace(-1, 1-2/MySizeCsiOSZ(1), MySizeCsiOSZ(1)).';
            end
            data(iAQSlice).ImageCsiFrequencyZero(:,1,1,1,iImage,iPart) = ImageCsiFrequencyZero;

            if iImage == 1 && iPart == 1
              data(iAQSlice).kSpaceOsRawFrequency = zeros([MySizeCsiOS nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).kSpaceOsRaw =          zeros([MySizeOS nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).ImageZ =               zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
              data(iAQSlice).kSpaceZ =              zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            end

            data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart) = ...
              reshape(mean(data(iAQSlice).data(AQSlice(iAQSlice).UseForImageStartIdx:AQSlice(iAQSlice).UseForImageStopIdx,UseAQWindow_tRep), 1), MySizeOS);
            data(iAQSlice).kSpaceOsRawFrequency(:,:,:,:,iImage,iPart) = reshape(data(iAQSlice).data(UseSample,UseAQWindow_tRep), MySizeCsiOS);

          else
            if iImage == 1 && iPart == 1
              data(iAQSlice).kSpaceOsRawFrequency = zeros([MySizeCsiOS nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).kSpaceOsRaw =          zeros([MySizeOS nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).ImageZ =               zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
              data(iAQSlice).kSpaceZ =              zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            end
            data(iAQSlice).kSpaceOsRawFrequency(:,:,:,:,iImage,iPart) = reshape(data(iAQSlice).data(UseSample,UseAQWindow_tRep), MySizeCsiOS);
            data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart) = reshape(mean(data(iAQSlice).data(UseSample,UseAQWindow_tRep), 1), MySizeOS);

          end

          % Calculate mean frequency per Voxel
          if length(UseSample) > 1
            if iImage == 1 && iPart == 1
              data(iAQSlice).ImageOsRawFftPhase = zeros([MySizeCsiOS nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).kSpaceOsZPhase = zeros([MySizeCsiOSZ nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).ImageOsZFftPhase = zeros([MySizeCsiOSZ nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).ImageOsFrequency = zeros([MySizeOS nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).ImageOsZFrequency = zeros([MySizeOSZ nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).ImageFrequency = zeros([MySize nImages nParts], 'like', data(iAQSlice).data);
              data(iAQSlice).ImageZFrequency = zeros([MySizeZ nImages nParts], 'like', data(iAQSlice).data);
            end
            for t = 1:size(data(iAQSlice).kSpaceOsRawFrequency,1)
              data(iAQSlice).ImageOsRawFftPhase(t,:,:,:,iImage,iPart) = fftshift(ifftn(ifftshift(data(iAQSlice).kSpaceOsRawFrequency(t,:,:,:,iImage,iPart))));
              data(iAQSlice).kSpaceOsZPhase(t,:,:,:,iImage,iPart) = zeroFill(data(iAQSlice).kSpaceOsRawFrequency(t,:,:,:,iImage,iPart), [1,MySizeOSZ(2:end)], [1,AQSlice(iAQSlice).ZeroFillWindowSize(2:end)]);
              data(iAQSlice).ImageOsZFftPhase(t,:,:,:,iImage,iPart) = fftshift(ifftn(ifftshift(data(iAQSlice).kSpaceOsZPhase(t,:,:,:,iImage,iPart))));
            end
            data(iAQSlice).ImageOsFrequency(:,:,:,:,iImage,iPart) = -get_MeanPhaseDiffWeighted(data(iAQSlice).ImageOsRawFftPhase(:,:,:,:,iImage,iPart), 1)/2/pi ...
              / squeeze(diff(data(iAQSlice).time_of_tRep(UseSample(1:2),UseAQWindow_tRep(1))));
            data(iAQSlice).ImageOsZFrequency(:,:,:,:,iImage,iPart) = -get_MeanPhaseDiffWeighted(data(iAQSlice).ImageOsZFftPhase(:,:,:,:,iImage,iPart), 1)/2/pi ...
              / squeeze(diff(data(iAQSlice).time_of_tRep(UseSample(1:2),UseAQWindow_tRep(1))));

            data(iAQSlice).ImageZFrequency(:,:,:,:,iImage,iPart) = data(iAQSlice).ImageOsZFrequency(1,...
                                                          floor(MySizeOSZ(2)/2) + (1-floor(MySizeZ(2)/2):ceil(MySizeZ(2)/2)),...
                                                          floor(MySizeOSZ(3)/2) + (1-floor(MySizeZ(3)/2):ceil(MySizeZ(3)/2)),...
                                                          floor(MySizeOSZ(4)/2) + (1-floor(MySizeZ(4)/2):ceil(MySizeZ(4)/2)),...
                                                          iImage,iPart);
            data(iAQSlice).ImageFrequency(:,:,:,:,iImage,iPart) = data(iAQSlice).ImageOsFrequency(1,...
                                                        floor(MySizeOS(2)/2) + (1-floor(MySize(2)/2):ceil(MySize(2)/2)),...
                                                        floor(MySizeOS(3)/2) + (1-floor(MySize(3)/2):ceil(MySize(3)/2)),...
                                                        floor(MySizeOS(4)/2) + (1-floor(MySize(4)/2):ceil(MySize(4)/2)),...
                                                        iImage,iPart);
          end

          if iImage == 1 && iPart == 1
            data(iAQSlice).kSpaceOsZ = zeros([MySizeOSZ nImages nParts], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpaceOs = zeros([MySizeOS nImages nParts], 'like', data(iAQSlice).data);
          end
          % zero-fill and smooth and convert to signal density
          data(iAQSlice).kSpaceOsZ(:,:,:,:,iImage,iPart) = ...
            zeroFill(data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart), MySizeOSZ, [Inf, AQSlice(iAQSlice).ZeroFillWindowSize(2:end)]) ...
            / AQSlice(iAQSlice).VoxelVolume * prod(data(iAQSlice).ZeroFillFactorOS) * AQSlice(iAQSlice).AreaCoil;
          % apply k-space filter (ZeroFillWindowSize)
          data(iAQSlice).kSpaceOs(:,:,:,:,iImage,iPart) = ...
            zeroFill(data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart), MySizeOS, [Inf, AQSlice(iAQSlice).ZeroFillWindowSize(2:end)]) ...
            / AQSlice(iAQSlice).VoxelVolume * AQSlice(iAQSlice).AreaCoil;

          if iImage == 1 && iPart == 1
            data(iAQSlice).ImageOsRaw = zeros([MySizeOS nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageOsZ = zeros([MySizeOSZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageOs = zeros([MySizeOS nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpaceZ = zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageZ = zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpace = zeros([MySize nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).Image = zeros([MySize nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
          end

          % Calculate oversampled Image
          data(iAQSlice).ImageOsRaw(:,:,:,:,iImage,iPart) = ...
            fftshift(ifftn(ifftshift(data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart))));
          if AQSlice(iAQSlice).partsAverage == nParts
            data(iAQSlice).ImageOsZ(:,:,:,:,iImage,iPart) = ...
              fftshift(ifftn(ifftshift(data(iAQSlice).kSpaceOsZ(:,:,:,:,iImage,iPart))));
            data(iAQSlice).ImageOs(:,:,:,:,iImage,iPart) = ...
              data(iAQSlice).ImageOsRaw(:,:,:,:,iImage,iPart);
          end

          if AQSlice(iAQSlice).partsAverage == nParts
            if isemptyfield(AQSlice(iAQSlice), 'ReadRadial'), AQSlice(iAQSlice).ReadRadial = 0; end
            if isemptyfield(data(iAQSlice), 'kSpaceZ'), data(iAQSlice).kSpaceZ = zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data); end
            if isemptyfield(data(iAQSlice), 'Image'),   data(iAQSlice).Image = zeros([MySize nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data); end
            if isemptyfield(data(iAQSlice), 'kSpace'),  data(iAQSlice).kSpace = zeros([MySize nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data); end
            data(iAQSlice) = cutImageOs(data(iAQSlice), AQSlice(iAQSlice), MySizeOSZ, MySizeZ, MySizeOS, MySize, iImage, iPart);
          end

        else
          %% read and phase encoding
          if isscalar(AQSlice(iAQSlice).ZeroFillFactor)
            AQSlice(iAQSlice).ZeroFillFactor = ...
              AQSlice(iAQSlice).ZeroFillFactor([1 1 1 1]);
          end
          if isscalar(AQSlice(iAQSlice).ZeroFillWindowSize)
            AQSlice(iAQSlice).ZeroFillWindowSize = ...
              AQSlice(iAQSlice).ZeroFillWindowSize([1 1 1 1]);
          end
          AQSlice(iAQSlice).ZeroFillFactor(end+1:4) = 1;
          AQSlice(iAQSlice).ZeroFillWindowSize(end+1:4) = Inf;
          AQSlice(iAQSlice).ZeroFillFactor(MySize<=1) = 1;
          AQSlice(iAQSlice).ZeroFillWindowSize(MySize<=1) = Inf;

          MySizeOSZ = ceil(MySizeOS.*(AQSlice(iAQSlice).ZeroFillFactor));
          MySizeZ = ceil(MySize.*(AQSlice(iAQSlice).ZeroFillFactor));
          data(iAQSlice).ZeroFillFactorOS = MySizeOSZ./MySizeOS;
          data(iAQSlice).ZeroFillFactor = MySizeZ./MySize;

          if iImage == 1 && iPart == 1
            data(iAQSlice).kSpaceOsRaw = zeros([MySizeOS nImages nParts], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageOsRaw = zeros([MySizeOS nImages nParts], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpaceOsZ = zeros([MySizeOSZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageOsZ = zeros([MySizeOSZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpaceOs = zeros([MySizeOS nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageOs = zeros([MySizeOS nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpaceZ = zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).ImageZ = zeros([MySizeZ nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).kSpace = zeros([MySize nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
            data(iAQSlice).Image = zeros([MySize nImages AQSlice(iAQSlice).partsAverage], 'like', data(iAQSlice).data);
          end
          usedkLine = UseAQWindow(:,iImage,iPart)~=0 & UsetRep(:,iImage,iPart)~=0;
          UseAQWindow_tRep = sub2ind(szData(2:3), UseAQWindow(usedkLine,iImage,iPart), UsetRep(usedkLine,iImage,iPart));
          kSpaceOsRaw = NaN(MySizeOS);
          kSpaceOsRaw(:,sort(AQSlice(iAQSlice).kLineOrder)) = ...
            reshape(data(iAQSlice).data(UseSample,UseAQWindow_tRep), [MySizeOS(1), numel(AQSlice(iAQSlice).kLineOrder)]);
          if isfield(AQSlice(iAQSlice), 'ReadDir') && size(AQSlice(iAQSlice).ReadDir, 1) == numel(UseAQWindow_tRep)
            kSpaceOsRaw(:,AQSlice(iAQSlice).ReadDir(:,iImage,iPart)<=0) = flip(kSpaceOsRaw(:,AQSlice(iAQSlice).ReadDir(:,iImage,iPart)<=0), 1);
          end
          kSpaceOsRaw = kSpaceOsRaw(AQSlice(iAQSlice).AQWindowSampleOrder);
          data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart) = reshape(kSpaceOsRaw, MySizeOS);
          if ~isemptyfield(data(iAQSlice), 'fft1_data') && ~isemptyfield(data(iAQSlice), 'f_fft1_data')
            fft1_data = NaN(MySizeOS, 'like', data(iAQSlice).data);
            fft1_data(:,sort(AQSlice(1).kLineOrder(:))) = ...
              data(iAQSlice).fft1_data(UseSample,UseAQWindow_tRep);
            data(iAQSlice).fft1_dataCut(:,:,:,:,iImage,iPart) = fft1_data;
            f_fft1_data = NaN(MySizeOS);
            f_fft1_data(:,sort(AQSlice(1).kLineOrder(:))) = ...
              data(iAQSlice).f_fft1_data(UseSample,UseAQWindow_tRep);
            data(iAQSlice).f_fft1_dataCut(:,:,:,:,iImage,iPart) = f_fft1_data;
          else
            data(iAQSlice).fft1_dataCut = [];
            data(iAQSlice).f_fft1_dataCut = [];
          end

          % zero-fill and smooth and convert to signal density
          % FIXME: Ist es sinnvoll den k-Raum zu skalieren? Oder besser nur
          % Bild-Raum?
          data(iAQSlice).kSpaceOsZ(:,:,:,:,iImage,iPart) = ...
            zeroFill(data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart), MySizeOSZ, AQSlice(iAQSlice).ZeroFillWindowSize) ...
            / AQSlice(iAQSlice).VoxelVolume * prod(data(iAQSlice).ZeroFillFactorOS) * AQSlice(iAQSlice).AreaCoil;
          data(iAQSlice).kSpaceOs(:,:,:,:,iImage,iPart) = ...
            zeroFill(data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart), MySizeOS, AQSlice(iAQSlice).ZeroFillWindowSize) ...
            / AQSlice(iAQSlice).VoxelVolume * AQSlice(iAQSlice).AreaCoil;

          % CIC correction
          N = data(iAQSlice).cic_N;
          M = data(iAQSlice).cic_M;
          n = permute(size(data(iAQSlice).kSpaceOsZ(:,:,:,:,iImage,iPart), 1), [4,1,2,3]);
          R = permute(1000, [4,1,2,3]);
          f = get_FFTGrid(1./R, n);
          R = repmat(R, [size(f,1),1,1,1]);
          Hf = (sin(pi.*R.*M.*f)./(sin( pi.*f))).^N;
          Hf(isnan(Hf)) = (R(isnan(Hf)).*M).^N;
          data(iAQSlice).cic_corrOsZ = ((R.*M).^N)./Hf;  % FIXME: How large is this going to be?
          firstkLineIdx = find(UseAQWindow(:,iImage,iPart)~=0 & UsetRep(:,iImage,iPart)~=0, 1, 'first');
          data(iAQSlice).cic_corrOs = data(iAQSlice).cic_corr(UseSample,UseAQWindow(firstkLineIdx,iImage,iPart),UsetRep(firstkLineIdx,iImage,iPart));

          % Calculate oversampled Image
          data(iAQSlice).ImageOsRaw(:,:,:,:,iImage,iPart) = ...
            fftshift(ifftn(ifftshift(data(iAQSlice).kSpaceOsRaw(:,:,:,:,iImage,iPart))));
          if AQSlice(iAQSlice).partsAverage == nParts
            data(iAQSlice).ImageOsZ(:,:,:,:,iImage,iPart) = ...
              bsxfun(@times, fftshift(ifftn(ifftshift(data(iAQSlice).kSpaceOsZ(:,:,:,:,iImage,iPart)))), data(iAQSlice).cic_corrOsZ);
            cic_corr = NaN(MySizeOS);
            cic_corr(:,sort(AQSlice(1).kLineOrder(:))) = ...
              data(iAQSlice).cic_corr(UseSample,UseAQWindow_tRep);
            data(iAQSlice).ImageOs(:,:,:,:,iImage,iPart) = ...
              fftshift(ifftn(ifftshift(data(iAQSlice).kSpaceOs(:,:,:,:,iImage,iPart)))) .* cic_corr;
          end

          if AQSlice(iAQSlice).partsAverage == nParts
            if isemptyfield(AQSlice(iAQSlice), 'ReadRadial'), AQSlice(iAQSlice).ReadRadial = 0; end
            data(iAQSlice) = cutImageOs(data(iAQSlice), AQSlice(iAQSlice), MySizeOSZ, MySizeZ, MySizeOS, MySize, iImage, iPart);
          end
        end
      end

    end
  end

  if AQSlice(iAQSlice).partsAverage ~= nParts
    %% Calculate average over several image parts

    for iPart = 1:AQSlice(iAQSlice).partsAverage
      PartIndex = AQSlice(1).partsAverageOrder(:,iPart);
      PartIndex = PartIndex(PartIndex>0);
      data(iAQSlice).kSpaceOsZ(:,:,:,:,:,iPart) = mean(data(iAQSlice).kSpaceOsZ(:,:,:,:,:,PartIndex), 6);
      data(iAQSlice).kSpaceOs(:,:,:,:,:,iPart) = mean(data(iAQSlice).kSpaceOs(:,:,:,:,:,PartIndex), 6);
    end
    if ~isempty(iPart)
      % remove data that was used for averaging
      data(iAQSlice).kSpaceOsZ(:,:,:,:,:,iPart+1:end) = [];
      data(iAQSlice).kSpaceOs(:,:,:,:,:,iPart+1:end) = [];
    end
    for iImage = 1:nImages
      for iPart = 1:AQSlice(iAQSlice).partsAverage
        % calculate over-sampled image from averaged k-space
        data(iAQSlice).ImageOsZ(1:MySizeOSZ(1),1:MySizeOSZ(2),1:MySizeOSZ(3),1:MySizeOSZ(4),iImage,iPart) = ...
          fftshift(ifftn(ifftshift(data(iAQSlice).kSpaceOsZ(:,:,:,:,iImage,iPart))));
        data(iAQSlice).ImageOs(1:MySizeOS(1),1:MySizeOS(2),1:MySizeOS(3),1:MySizeOS(4),iImage,iPart) = ...
          fftshift(ifftn(ifftshift(data(iAQSlice).kSpaceOs(:,:,:,:,iImage,iPart))));
        % FIXME: Complete implementation of averaging for CSI. There are more
        % fields in data that need special treatment.
        if ~(AQSlice(iAQSlice).nRead == 1 && prod(AQSlice(iAQSlice).nPhase) > 1) % read and phase encoding
          % CIC filter correction along read
          data(iAQSlice).ImageOsZ(1:MySizeOSZ(1),1:MySizeOSZ(2),1:MySizeOSZ(3),1:MySizeOSZ(4),iImage,iPart) = ...
            bsxfun(@times, data(iAQSlice).ImageOsZ(:,:,:,:,iImage,iPart), data(iAQSlice).cic_corrOsZ);
          data(iAQSlice).ImageOs(1:MySizeOS(1),1:MySizeOS(2),1:MySizeOS(3),1:MySizeOS(4),iImage,iPart) = ...
            bsxfun(@times, data(iAQSlice).ImageOs(:,:,:,:,iImage,iPart), data(iAQSlice).cic_corrOs);
        end
        % Cut image from over-sampled image and re-construct matching k-space
        data(iAQSlice) = cutImageOs(data(iAQSlice), AQSlice(iAQSlice), MySizeOSZ, MySizeZ, MySizeOS, MySize, iImage, iPart);
      end
    end
  end
end

if nargout == 2
  dataOut.Image = data.Image;
  dataOut.ImageZ = data.ImageZ;
  % dataOut.ImageOs = data.ImageOs;
  % dataOut.ImageOsRaw = data.ImageOsRaw;
  dataOut.kSpace = data.kSpace;
  dataOut.kSpaceOs = data.kSpaceOs;
  % dataOut.kSpaceOsRaw = data.kSpaceOsRaw;
  dataOut.ZeroFillFactorOS = data.ZeroFillFactorOS;
  dataOut.ZeroFillFactor = data.ZeroFillFactor;
  dataOut.Amplitude2Uin = data.Amplitude2Uin;
  if AQSlice.plotFft1_data ~= 0
    dataOut.fft1_dataCut = data.fft1_dataCut;
    dataOut.f_fft1_dataCut = data.f_fft1_dataCut;
  end
end

end


function data = cutImageOs(data, AQSlice, MySizeOSZ, MySizeZ, MySizeOS, MySize, iImage, iPart)
%% Cut image from over-sampled image and re-construct matching k-space

if ~AQSlice.ReadRadial
  % cut oversampled data (ReadOS, PhaseOS and PhaseOS )
  if isfield(data, 'ImageOsZ')
    data.ImageZ(:,:,:,:,iImage,iPart) = data.ImageOsZ(floor(MySizeOSZ(1)/2) + (1-floor(MySizeZ(1)/2):ceil(MySizeZ(1)/2)),...
                                                      floor(MySizeOSZ(2)/2) + (1-floor(MySizeZ(2)/2):ceil(MySizeZ(2)/2)),...
                                                      floor(MySizeOSZ(3)/2) + (1-floor(MySizeZ(3)/2):ceil(MySizeZ(3)/2)),...
                                                      floor(MySizeOSZ(4)/2) + (1-floor(MySizeZ(4)/2):ceil(MySizeZ(4)/2)),...
                                                      iImage,iPart);
    data.kSpaceZ(:,:,:,:,iImage,iPart) = fftshift(fftn(ifftshift(data.ImageZ(:,:,:,:,iImage,iPart))));
  end
  data.Image(:,:,:,:,iImage,iPart)    = data.ImageOs( floor(MySizeOS(1)/2) + (1-floor(MySize(1)/2):ceil(MySize(1)/2)),...
                                                  floor(MySizeOS(2)/2) + (1-floor(MySize(2)/2):ceil(MySize(2)/2)), ...
                                                  floor(MySizeOS(3)/2) + (1-floor(MySize(3)/2):ceil(MySize(3)/2)), ...
                                                  floor(MySizeOS(4)/2) + (1-floor(MySize(4)/2):ceil(MySize(4)/2)),...
                                                  iImage,iPart);
  data.kSpace(:,:,:,:,iImage,iPart) = fftshift(fftn(ifftshift(data.Image(:,:,:,:,iImage,iPart))));

  % if AQSlice.nRead == 1
  %   data.kSpaceOs=data.kSpaceOsRaw;
  %
  % elseif numel(UseSample)==MySizeOS(1)
  %   data.kSpaceOs=fftshift(fft(ifftshift(reshape(data.fft1_data(UseSample,UseAQWindow,UsetRep), ...
  %                                                MySizeOS), 1), [], 1), 1);
  % elseif isscalar(UseSample)
  %   data.kSpaceOs=data.kSpaceOsRaw;
  % else
  %   warning('no CIC correction available')
  %   data.kSpaceOs=data.kSpaceOsRaw;
  %   % data_S_AQs_TRs.fft1_data(1:AQ.nSamples(t,TR),t,TR)=(fftshift(ifft(ifftshift(raw_data_AQs(1:AQ.nSamples(t,TR),nAQ).*AQ.Norm2Amplitude(TR).*AQ.rawData2Norm)))).*data_S_AQs_TRs.cic_corr(1:AQ.nSamples(t,TR),t,TR);
  %   % fftshift(ifft(ifftshift(data.data(UseSample,UseAQWindow,UsetRep)).*data.cic_corr(UseSample,UseAQWindow,UsetRep
  % end
else
  if AQSlice.ZeroFillFactor == 1
    if ~isemptyfield(data, 'ImageOsZ')
      data.ImageZ(:,:,:,:,iImage,iPart) = data.ImageOsZ(:,:,:,:,iImage,iPart);
      data.kSpaceZ(:,:,:,:,iImage,iPart) = fftshift(fftn(ifftshift(data.ImageZ(:,:,:,:,iImage,iPart))));
    end
    data.Image(:,:,:,:,iImage,iPart) = data.ImageOs(:,:,:,:,iImage,iPart);
    data.kSpace(:,:,:,:,iImage,iPart) = fftshift(fftn(ifftshift(data.Image(:,:,:,:,iImage,iPart))));
  end
end

end
