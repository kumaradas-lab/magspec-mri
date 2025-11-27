classdef TXMaxDef < handle
  %% Class managing the Max and Def properties of the class PD.TXClass
  %
  %     PD.TXMaxDef(TX, type, freq)
  %
  % ----------------------------------------------------------------------------
  % (C) Copyright 2020-2025 Pure Devices GmbH, Wuerzburg, Germany
  % www.pure-devices.com
  % ----------------------------------------------------------------------------

  properties (Dependent)
    %% input values
    Norm
    Mmrt
    Uout
    PaUout
    Amplitude
  end


  properties (SetAccess = private, Dependent)
    LimitMode
  end


  properties (SetAccess = private)
    %% output values
    NormCalibrated
    MmrtUoutCalibrated
    UoutCalibrated
    PaUoutCalibrated
    AmplitudeCalibrated
  end


  properties (GetAccess = private, SetAccess = immutable)
    %% immutable internal settings
    TX
    typeMax  % 0: Def, 1: Max

    Frequency

    fLarmor
    fLarmorX
  end


  properties (GetAccess = private, SetAccess = ?PD.TXClass)
    %% private property that is allowed to be updated on fSystem change (by HW.TX)
    FrequencyGrid
  end


  properties (GetAccess = private, SetAccess = private)
    %% calculate set of output values for each input value
    % The matrices in this section contain 5 rows. They correspond to the values
    % derived from the following values:
    % 1: Norm
    % 2: Mmrt
    % 3: Uout
    % 4: PaUout
    % 5: Amplitude
    PaUoutCalibratedAll
    UoutCalibratedAll
    MmrtUoutCalibratedAll
    NormCalibratedAll
    AmplitudeCalibratedAll
  end


  methods

    function this = TXMaxDef(tx, typ, freq)
      %% constructor for class PD.TXMaxDef

      this.TX = tx;
      this.typeMax = strcmp(typ, 'Max');
      typeCal = strcmp(typ, 'Cal');
      if nargin > 2
        % If this is set, the object doesn't change on Larmor frequency changes.
        this.Frequency = freq;
      end

      if isa(tx, 'PD.TXClass')
        % This only works for HW object!
        this.fLarmor = tx.HWInstance.fLarmor;
        this.fLarmorX = tx.HWInstance.fLarmorX;
      end

      if (isa(tx, 'PD.TXClass') && ~isempty(tx.fSample) && ~isempty(tx.DdsPicBits)) ...
          || (isstruct(tx) && ~isemptyfield(tx, 'fSample') && ~isemptyfield(tx, 'DdsPicBits'))
        this.FrequencyGrid = tx.fSample / 2^tx.DdsPicBits;
      else
        % FIXME: This is the default value. The actual value might be different
        %        for some configurations.
        this.FrequencyGrid = 125e6 / 2^32;
      end

      % default input values
      if typeCal
        n = this.TX.n;
        this.Norm = Inf(1, n);
        this.Mmrt = Inf(1, n);
        this.Uout = Inf(1, n);
        this.PaUout = Inf(1, n);
        this.Amplitude = Inf(1, n);
      end
    end


    function stru = struct(this)
      %% overload struct function

      stru.Norm = this.Norm;
      stru.Mmrt = this.Mmrt;
      stru.Uout = this.Uout;
      stru.PaUout = this.PaUout;
      stru.Amplitude = this.Amplitude;
      stru.PaUoutCalibrated = this.PaUoutCalibrated;
      stru.UoutCalibrated = this.UoutCalibrated;
      stru.MmrtUoutCalibrated = this.MmrtUoutCalibrated;
      stru.NormCalibrated = this.NormCalibrated;
    end


    function set.Amplitude(this, val)
      this.AmplitudeCalibratedAll(5,:) = val;

      this.UpdatePaUoutCalibratedByAmplitude(5);
      if ~this.typeMax
        this.TX.UpdateEstPaUoutCalibrated();
      end
      this.UpdateUoutCalibratedByPaUout(5);
      if ~this.typeMax
        this.TX.UpdateEstUoutCalibrated();
      end
      this.UpdateMmrtUoutCalibratedByUout(5);
      if ~this.typeMax
        this.TX.UpdateEstMmrtUoutCalibrated();
      end
      this.UpdateNormCalibratedByMmrtUout(5);
    end

    function val = get.Amplitude(this)
      if size(this.AmplitudeCalibratedAll, 1) < 5
        val = [];
      else
        val = this.AmplitudeCalibratedAll(5,:);
      end
    end


    function set.PaUout(this, val)
      this.PaUoutCalibratedAll(4,:) = val;

      this.UpdateAmplitudeCalibratedByPaUout(4);

      if ~this.typeMax
        this.TX.UpdateEstPaUoutCalibrated();
      end
      this.UpdateUoutCalibratedByPaUout(4);
      if ~this.typeMax
        this.TX.UpdateEstUoutCalibrated();
      end
      this.UpdateMmrtUoutCalibratedByUout(4);
      if ~this.typeMax
        this.TX.UpdateEstMmrtUoutCalibrated();
      end
      this.UpdateNormCalibratedByMmrtUout(4);


      this.TX.UpdateUout2PaUoutCalibrationGainMaxDef();  % FIXME: Needed?
    end

    function val = get.PaUout(this)
      if size(this.PaUoutCalibratedAll, 1) < 4
        val = [];
      else
        val = this.PaUoutCalibratedAll(4,:);
      end
    end


    function set.Uout(this, val)
      this.UoutCalibratedAll(3,:) = val;

      this.UpdatePaUoutCalibratedByUout(3);
      if ~this.typeMax
        this.TX.UpdateEstPaUoutCalibrated();
      end
      this.UpdateAmplitudeCalibratedByPaUout(3);

      if ~this.typeMax
        this.TX.UpdateEstUoutCalibrated();
      end
      this.UpdateMmrtUoutCalibratedByUout(3);
      if ~this.typeMax
        this.TX.UpdateEstMmrtUoutCalibrated();
      end
      this.UpdateNormCalibratedByMmrtUout(3);

      this.TX.UpdateUout2PaUoutCalibrationGainMaxDef();  % FIXME: Needed?
    end

    function val = get.Uout(this)
      if size(this.UoutCalibratedAll, 1) < 3
        val = [];
      else
        val = this.UoutCalibratedAll(3,:);
      end
    end


    function set.Mmrt(this, val)
      this.MmrtUoutCalibratedAll(2,:) = val;
      if ~this.typeMax
        this.TX.UpdateEstMmrtUoutCalibrated();
      end

      this.UpdateUoutCalibratedByMmrtUout(2);
      this.UpdatePaUoutCalibratedByUout(2);
      if ~this.typeMax
        this.TX.UpdateEstPaUoutCalibrated();
      end
      this.UpdateAmplitudeCalibratedByPaUout(2);

      this.UpdateNormCalibratedByMmrtUout(2);
    end

    function val = get.Mmrt(this)
      if size(this.MmrtUoutCalibratedAll, 1) < 2
        val = [];
      else
        val = this.MmrtUoutCalibratedAll(2,:);
      end
    end


    function set.Norm(this, val)
      this.NormCalibratedAll(1,:) = val;

      this.UpdateMmrtUoutCalibratedByNorm(1);
      if ~this.typeMax
        this.TX.UpdateEstMmrtUoutCalibrated();
      end
      this.UpdateUoutCalibratedByMmrtUout(1);
      this.UpdatePaUoutCalibratedByUout(1);
      if ~this.typeMax
        this.TX.UpdateEstPaUoutCalibrated();
      end
      this.UpdateAmplitudeCalibratedByPaUout(1);
    end

    function val = get.Norm(this)
      if size(this.NormCalibratedAll, 1) < 1
        val = [];
      else
      val = this.NormCalibratedAll(1,:);
      end
    end


    function val = get.NormCalibrated(this)
      this.NormCalibratedAll(1,:) = this.Norm;
      if any(this.NormCalibratedAll == 0)
        fieldNames = {'Norm', 'Mmrt', 'Uout', 'PaUout', 'Amplitude'};
        fieldOrder = [1, 2, 3, 4, 5];
        for iRow = fieldOrder
          if any(this.NormCalibratedAll(iRow,:)==0)
            % trigger the update functions (twice :-( )
            oldVal = this.(fieldNames{iRow});
            this.(fieldNames{iRow}) = oldVal+eps;
            this.(fieldNames{iRow}) = oldVal;
          end
        end
      end
      val = min(this.NormCalibratedAll, [], 1);
    end

    function val = get.MmrtUoutCalibrated(this)
      this.MmrtUoutCalibratedAll(2,:) = this.Mmrt;
      if any(this.MmrtUoutCalibratedAll == 0)
        fieldNames = {'Norm', 'Mmrt', 'Uout', 'PaUout', 'Amplitude'};
        fieldOrder = [2, 3, 4, 5, 1];
        for iRow = fieldOrder
          if any(this.MmrtUoutCalibratedAll(iRow,:)==0)
            oldVal = this.(fieldNames{iRow});
            this.(fieldNames{iRow}) = oldVal+eps;
            this.(fieldNames{iRow}) = oldVal;
          end
        end
      end
      val = min(this.MmrtUoutCalibratedAll, [], 1);
    end

    function val = get.UoutCalibrated(this)
      this.UoutCalibratedAll(3,:) = this.Uout;
      if any(this.UoutCalibratedAll == 0)
        fieldNames = {'Norm', 'Mmrt', 'Uout', 'PaUout', 'Amplitude'};
        fieldOrder = [3, 4, 5, 2, 1];
        for iRow = fieldOrder
          if any(this.UoutCalibratedAll(iRow,:)==0)
            oldVal = this.(fieldNames{iRow});
            this.(fieldNames{iRow}) = oldVal+eps;
            this.(fieldNames{iRow}) = oldVal;
          end
        end
      end
      val = min(this.UoutCalibratedAll, [], 1);
    end

    function val = get.PaUoutCalibrated(this)
      this.PaUoutCalibratedAll(4,:) = this.PaUout;
      if any(this.PaUoutCalibratedAll == 0)
        fieldNames = {'Norm', 'Mmrt', 'Uout', 'PaUout', 'Amplitude'};
        fieldOrder = [4, 5, 3, 2, 1];
        for iRow = fieldOrder
          if any(this.PaUoutCalibratedAll(iRow,:)==0)
            oldVal = this.(fieldNames{iRow});
            this.(fieldNames{iRow}) = oldVal+eps;
            this.(fieldNames{iRow}) = oldVal;
          end
        end
      end
      val = min(this.PaUoutCalibratedAll, [], 1);
    end

    function val = get.AmplitudeCalibrated(this)
      this.AmplitudeCalibratedAll(5,:) = this.Amplitude;
      if any(this.AmplitudeCalibratedAll == 0)
        fieldNames = {'Norm', 'Mmrt', 'Uout', 'PaUout', 'Amplitude'};
        fieldOrder = [5, 4, 3, 2, 1];
        for iRow = fieldOrder
          if any(this.AmplitudeCalibratedAll(iRow,:)==0)
            oldVal = this.(fieldNames{iRow});
            this.(fieldNames{iRow}) = oldVal+eps;
            this.(fieldNames{iRow}) = oldVal;
          end
        end
      end
      val = min(this.AmplitudeCalibratedAll, [], 1);
    end


    function val = get.LimitMode(this)
      %% Indicate which of the input values limit the output values
      fieldNames = {'Norm', 'Mmrt', 'Uout', 'PaUout', 'Amplitude'};
      val = cell(1, size(this.NormCalibratedAll, 2));
      for iChannel = 1:size(this.NormCalibratedAll, 2)
        val(iChannel) ...
          = fieldNames(find(this.NormCalibratedAll(:,iChannel) ...
                       == min(this.NormCalibratedAll(:,iChannel)), 1, 'first'));
      end
    end
  end


  methods (Access = private)

    function UpdatePaUoutCalibratedByAmplitude(this, row)
      %% step from Amplitude to PaUout
      if size(this.AmplitudeCalibratedAll, 1) < row || ...
          isempty(this.TX.PaUout2AmplitudeCalibrationGain) || ...
          isempty(this.TX.PaUout2Amplitude)
        return;
      end
      if isempty(this.TX.PaUout2AmplitudeX) || isequal(this.fLarmor, this.fLarmorX)
        useXNucleus = false;
      else
        useXNucleus = abs(this.Frequency - this.fLarmorX) < abs(this.Frequency - this.fLarmor);
      end
      if useXNucleus
        paUoutCalibrated = this.AmplitudeCalibratedAll(row,:) ./ ...
          this.TX.PaUout2AmplitudeCalibrationGain ./ this.TX.PaUout2AmplitudeX;
      else
        paUoutCalibrated = this.AmplitudeCalibratedAll(row,:) ./ ...
          this.TX.PaUout2AmplitudeCalibrationGain ./ this.TX.PaUout2Amplitude;
      end
      this.PaUoutCalibratedAll(row,:) = paUoutCalibrated;
      if ~isempty(this.PaUout)
        validIdx = isfinite(this.PaUout) & (this.PaUoutCalibratedAll(row,:)==0);
        this.PaUoutCalibratedAll(row,validIdx) = this.PaUout(validIdx);
      end
    end


    function UpdateAmplitudeCalibratedByPaUout(this, row)
      %% step from PaUout to Amplitude
      if size(this.PaUoutCalibratedAll, 1) < row || ...
          isempty(this.TX.PaUout2AmplitudeCalibrationGain) || ...
          isempty(this.TX.PaUout2Amplitude)
        return;
      end
      if isempty(this.TX.PaUout2AmplitudeX) || isequal(this.fLarmor, this.fLarmorX)
        useXNucleus = false;
      else
        useXNucleus = abs(this.Frequency - this.fLarmorX) < abs(this.Frequency - this.fLarmor);
      end
      if useXNucleus
        amplitudeCalibrated = this.PaUoutCalibratedAll(row,:) .* ...
          this.TX.PaUout2AmplitudeCalibrationGain .* this.TX.PaUout2AmplitudeX;
      else
        amplitudeCalibrated = this.PaUoutCalibratedAll(row,:) .* ...
          this.TX.PaUout2AmplitudeCalibrationGain .* this.TX.PaUout2Amplitude;
      end
      this.AmplitudeCalibratedAll(row,:) = amplitudeCalibrated;
      if ~isempty(this.Amplitude)
        validIdx = isfinite(this.Amplitude) & (this.AmplitudeCalibratedAll(row,:)==0);
        this.AmplitudeCalibratedAll(row,validIdx) = this.Amplitude(validIdx);
      end
    end


    function UpdateUoutCalibratedByPaUout(this, row)
      %% step from PaUout to Uout
      if size(this.PaUoutCalibratedAll, 1) < row || ...
          isempty(this.TX.Uout2PaUout)
        return;
      end

      % get calibration value for amplification of RF amplifier
      nChannels = size(this.UoutCalibratedAll, 2);
      uout2PaUoutCalibrationGain = ones(1, nChannels);
      paUoutCalibrated = this.PaUoutCalibratedAll(row,:);
      for iChannel = 1:nChannels
        if isfield(this.TX.CalibrationRfAmp, 'TriScatteredAmp') && ...
            ~isempty(this.TX.CalibrationRfAmp(iChannel).TriScatteredAmp)
          % FIXME: This might not limit the PaUout correctly at the current
          % frequency.
          validMap = ~isnan(this.TX.CalibrationRfAmp(iChannel).InputAmplitude);
          paUoutCalibrated(iChannel) = min(paUoutCalibrated(iChannel), ...
            max(abs(this.TX.CalibrationRfAmp(iChannel).OutputAmplitude(validMap))));
          if isempty(this.Frequency)
            frequency = this.TX.HWInstance.fLarmor;
          else
            frequency = this.Frequency;
          end
          frequency = round(frequency/this.FrequencyGrid) * this.FrequencyGrid;
          uoutCalibrated = ...
            this.TX.CalibrationRfAmp(iChannel).TriScatteredAmp(frequency, ...
                                                               paUoutCalibrated(iChannel));
          uout2PaUoutCalibrationGain(iChannel) = ...
            (paUoutCalibrated(iChannel) / uoutCalibrated) / this.TX.Uout2PaUout(iChannel);
        elseif isfield(this.TX.CalibrationRfAmp, 'Gain') && ...
            ~isempty(this.TX.CalibrationRfAmp(t).Gain) && ...
            isfield(this.TX.CalibrationRfAmp, 'Frequency') && ...
            ~isempty(this.TX.CalibrationRfAmp(t).Frequency)
          if isempty(this.Frequency)
            frequency = this.TX.HWInstance.fLarmor;
          else
            frequency = this.Frequency;
          end
          frequency = round(frequency/this.FrequencyGrid) * this.FrequencyGrid;
          uout2PaUoutCalibrationGain(iChannel) = ...
            interp1(this.TX.CalibrationRfAmp(iChannel).Frequency, ...
                    this.TX.CalibrationRfAmp(iChannel).Gain, ...
                    frequency, 'pchip', 1);
        end
      end

      uoutCalibrated = paUoutCalibrated ./ ...
        uout2PaUoutCalibrationGain ./ this.TX.Uout2PaUout;
      uoutCalibrated(isnan(uoutCalibrated)) = Inf;
      this.UoutCalibratedAll(row,:) = uoutCalibrated;
      if ~isempty(this.Uout)
        validIdx = isfinite(this.Uout) & (this.UoutCalibratedAll(row,:)==0);
        this.UoutCalibratedAll(row,validIdx) = this.Uout(validIdx);
      end
    end


    function UpdatePaUoutCalibratedByUout(this, row)
      %% step from Uout to PaUout
      if size(this.UoutCalibratedAll, 1) < row || ...
          isempty(this.TX.Uout2PaUout)
        return;
      end

      % get (inverse) calibration value for amplification of RF amplifier
      nChannels = size(this.UoutCalibratedAll, 2);
      uout2PaUoutCalibrationGain = ones(1, nChannels);
      uOutCalibrated = this.UoutCalibratedAll(row,:);
      for iChannel = 1:nChannels
        if isfield(this.TX.CalibrationRfAmp, 'TriScatteredAmp') && ...
            ~isempty(this.TX.CalibrationRfAmp(iChannel).TriScatteredAmp)
          validMap = ~isnan(this.TX.CalibrationRfAmp(iChannel).InputAmplitude);
          paUoutVec = linspace(0, max(abs(this.TX.CalibrationRfAmp(iChannel).OutputAmplitude(validMap))), 1000);
          if isempty(this.Frequency)
            frequency = this.TX.HWInstance.fLarmor;
          else
            frequency = this.Frequency;
          end
          frequency = round(frequency/this.FrequencyGrid) * this.FrequencyGrid;
          uoutCalibratedVec = ...
            this.TX.CalibrationRfAmp(iChannel).TriScatteredAmp(frequency*ones(size(paUoutVec)), ...
                                                               paUoutVec);
          if this.UoutCalibratedAll(row,iChannel) > max(uoutCalibratedVec)
            [uOutCalibrated(iChannel), imax] = max(uoutCalibratedVec);
            paUoutCalibrated = paUoutVec(imax);
            % That means that we'll "fake" a calibration factor to have the
            % maximum output amplitude of the amplifier.
          elseif this.UoutCalibratedAll(row,iChannel) < min(uoutCalibratedVec)
            [uOutCalibrated(iChannel), imin] = min(uoutCalibratedVec);
            paUoutCalibrated = paUoutVec(imin);
          else
            isValid = ~isnan(uoutCalibratedVec);
            if any(isValid)
              paUoutCalibrated = ...
                interp1(uoutCalibratedVec(isValid), ...
                        paUoutVec(isValid), ...
                        min(this.UoutCalibratedAll(row,iChannel), 1e5));
            else
              paUoutCalibrated = NaN;
            end
          end
          uout2PaUoutCalibrationGain(iChannel) = ...
            (paUoutCalibrated / uOutCalibrated(iChannel)) / this.TX.Uout2PaUout(iChannel);

          uout2PaUoutCalibrationGain(uOutCalibrated==0) = 0;

        elseif isfield(this.TX.CalibrationRfAmp, 'Gain') && ...
            ~isempty(this.TX.CalibrationRfAmp(t).Gain) && ...
            isfield(this.TX.CalibrationRfAmp, 'Frequency') && ...
            ~isempty(this.TX.CalibrationRfAmp(t).Frequency)
          if isempty(this.Frequency)
            frequency = this.TX.HWInstance.fLarmor;
          else
            frequency = this.Frequency;
          end
          frequency = round(frequency/this.FrequencyGrid) * this.FrequencyGrid;
          uout2PaUoutCalibrationGain(iChannel) = ...
            interp1(this.TX.CalibrationRfAmp(iChannel).Frequency, ...
            this.TX.CalibrationRfAmp(iChannel).Gain, frequency, 'pchip', 1);
        end
      end

      paUoutCalibrated = uOutCalibrated .* ...
        uout2PaUoutCalibrationGain .* this.TX.Uout2PaUout;
      paUoutCalibrated(isnan(paUoutCalibrated)) = Inf;
      this.PaUoutCalibratedAll(row,:) = paUoutCalibrated;
      if ~isempty(this.PaUout)
        validIdx = isfinite(this.PaUout) & this.PaUoutCalibratedAll(row,:)==0;
        this.PaUoutCalibratedAll(row,validIdx) = this.PaUout(validIdx);
      end
    end


    function UpdateMmrtUoutCalibratedByUout(this, row)
      %% step from Uout to MmrtUout
      if size(this.UoutCalibratedAll, 1) < row || ...
          isempty(this.TX.MmrtUout2UoutCalibrationGain) || ...
          isempty(this.TX.MmrtUout2Uout)
        return;
      end
      nChannels = size(this.MmrtUoutCalibratedAll, 2);
      mmrtUout2UoutCalibrationGain = ones(1, nChannels);
      % mmrtUout2UoutCalibrationPhaseOffset = zeros(1, nChannels);
      uoutCalibrated = this.UoutCalibratedAll(row,:);
      for iChannel = 1:nChannels
        if isfield(this.TX.CalibrationUout, 'TriScatteredAmp') && ...
            ~isempty(this.TX.CalibrationUout(iChannel).TriScatteredAmp)
          % gain (and phase) dependent on frequency and amplitude
          % FIXME: This is not yet tested.
          % FIXME: This might not limit the PaUout correctly at the current
          % frequency.
          validMap = ~isnan(this.TX.CalibrationUout(iChannel).InputAmplitude);
          uoutCalibrated(iChannel) = min(uoutCalibrated(iChannel), ...
            max(abs(this.TX.CalibrationUout(iChannel).OutputAmplitude(validMap))));
          if isempty(this.Frequency)
            frequency = this.TX.HWInstance.fLarmor;
          else
            frequency = this.Frequency;
          end
          frequency = round(frequency/this.FrequencyGrid) * this.FrequencyGrid;
          mmrtUoutCalibrated = ...
            this.TX.CalibrationUout(iChannel).TriScatteredAmp(frequency, ...
                                                              uoutCalibrated(iChannel));
          mmrtUout2UoutCalibrationGain(iChannel) = ...
            (uoutCalibrated(iChannel) / mmrtUoutCalibrated) / this.TX.MmrtUout2Uout(iChannel);

          % mmrtUout2UoutCalibrationPhaseOffset = ...
          %   this.TX.CalibrationUout(iChannel).TriScatteredPhase(frequency, ...
          %                                                       uoutCalibrated(iChannel));
        elseif isfield(this.TX.CalibrationUout, 'Gain') && ...
            ~isempty(this.TX.CalibrationUout(iChannel).Gain) && ...
            isfield(this.TX.CalibrationUout, 'Frequency') && ...
            ~isempty(this.TX.CalibrationUout(iChannel).Frequency)
          if isempty(this.Frequency)
            frequency = this.TX.HWInstance.fLarmor;
          else
            frequency = this.Frequency;
          end
          frequency = round(frequency/this.FrequencyGrid) * this.FrequencyGrid;
          mmrtUout2UoutCalibrationGain(iChannel) = ...
            interp1(this.TX.CalibrationUout(iChannel).Frequency, ...
                    this.TX.CalibrationUout(iChannel).Gain, frequency, 'pchip', 1);
          % mmrtUout2UoutCalibrationPhaseOffset(iChannel) = ...
          %   interp1(this.TX.CalibrationUout(iChannel).Frequency, ...
          %           unwrap(this.TX.CalibrationUout(iChannel).Phase./180.*pi).*180./pi, ...
          %           frequency, 'pchip', 0);
        end
      end

      mmrtUoutCalibrated = this.UoutCalibratedAll(row,:) ./ ...
        mmrtUout2UoutCalibrationGain ./ this.TX.MmrtUout2Uout;
      this.MmrtUoutCalibratedAll(row,:) = mmrtUoutCalibrated;
      if ~isempty(this.Mmrt)
        validIdx = isfinite(this.Mmrt) & this.MmrtUoutCalibratedAll(row,:)==0;
        this.MmrtUoutCalibratedAll(row,validIdx) = this.Mmrt(validIdx);
      end
    end


    function UpdateUoutCalibratedByMmrtUout(this, row)
      %% step from MmrtUout to Uout
      if size(this.MmrtUoutCalibratedAll, 1) < row || ...
          isempty(this.TX.MmrtUout2UoutCalibrationGain) || ...
          isempty(this.TX.MmrtUout2Uout)
        return;
      end
      nChannels = size(this.MmrtUoutCalibratedAll, 2);
      mmrtUout2UoutCalibrationGain = ones(1, nChannels);
      % mmrtUout2UoutCalibrationPhaseOffset = zeros(1, nChannels);
      mmrtUoutCalibrated = this.MmrtUoutCalibratedAll(row,:);
      for iChannel = 1:nChannels
        if isfield(this.TX.CalibrationUout, 'TriScatteredAmp') && ...
            ~isempty(this.TX.CalibrationUout(iChannel).TriScatteredAmp)
          % gain (and phase) dependent on frequency and amplitude
          % FIXME: This is not yet tested.

          validMap = ~isnan(this.TX.CalibrationUout(iChannel).InputAmplitude);
          uoutVec = linspace(0, max(abs(this.TX.CalibrationUout(iChannel).OutputAmplitude(validMap))), 1000);
          if isempty(this.Frequency)
            frequency = this.TX.HWInstance.fLarmor;
          else
            frequency = this.Frequency;
          end
          frequency = round(frequency/this.FrequencyGrid) * this.FrequencyGrid;
          mmrtUoutCalibratedVec = ...
            this.TX.CalibrationUout(iChannel).TriScatteredAmp(frequency*ones(size(uoutVec)), ...
                                                              uoutVec);
          if this.MmrtUoutCalibratedAll(row,iChannel) > max(mmrtUoutCalibratedVec)
            [mmrtUoutCalibrated(iChannel), imax] = max(mmrtUoutCalibratedVec);
            uoutCalibrated = uoutVec(imax);
            % That means that we'll "fake" a calibration factor to have the
            % maximum output amplitude of the amplifier.
          elseif this.MmrtUoutCalibratedAll(row,iChannel) < min(mmrtUoutCalibratedVec)
            [mmrtUoutCalibrated(iChannel), imin] = min(mmrtUoutCalibratedVec);
            uoutCalibrated = uoutVec(imin);
          else
            isValid = ~isnan(mmrtUoutCalibratedVec);
            if any(isValid)
              uoutCalibrated = ...
                interp1(mmrtUoutCalibratedVec(isValid), ...
                        uoutVec(isValid), ...
                        min(this.MmrtUoutCalibratedAll(row,iChannel), 1e5));
            else
              uoutCalibrated = NaN;
            end
          end
          mmrtUout2UoutCalibrationGain(iChannel) = ...
            (uoutCalibrated / mmrtUoutCalibrated(iChannel)) / this.TX.MmrtUout2Uout(iChannel);

          mmrtUout2UoutCalibrationGain(mmrtUoutCalibrated==0) = 0;

        elseif isfield(this.TX.CalibrationUout, 'Gain') && ...
            ~isempty(this.TX.CalibrationUout(iChannel).Gain) && ...
            isfield(this.TX.CalibrationUout, 'Frequency') && ...
            ~isempty(this.TX.CalibrationUout(iChannel).Frequency)
          if isempty(this.Frequency)
            frequency = this.TX.HWInstance.fLarmor;
          else
            frequency = this.Frequency;
          end
          frequency = round(frequency/this.FrequencyGrid) * this.FrequencyGrid;
          mmrtUout2UoutCalibrationGain(iChannel) = ...
            interp1(this.TX.CalibrationUout(iChannel).Frequency, ...
                    this.TX.CalibrationUout(iChannel).Gain, frequency, 'pchip', 1);
          % mmrtUout2UoutCalibrationPhaseOffset(iChannel) = ...
          %   interp1(this.TX.CalibrationUout(iChannel).Frequency, ...
          %           unwrap(this.TX.CalibrationUout(iChannel).Phase./180.*pi).*180./pi, ...
          %           frequency, 'pchip', 0);
        end
      end

      uoutCalibrated = mmrtUoutCalibrated .* ...
        mmrtUout2UoutCalibrationGain .* this.TX.MmrtUout2Uout;
      this.UoutCalibratedAll(row,:) = uoutCalibrated;
      if ~isempty(this.Uout)
        validIdx = isfinite(this.Uout) & this.UoutCalibratedAll(row,:)==0;
        this.UoutCalibratedAll(row,validIdx) = this.Uout(validIdx);
      end
    end


    function UpdateNormCalibratedByMmrtUout(this, row)
      %% step from MmrtUout to normalized DAC
      if size(this.MmrtUoutCalibratedAll, 1) < row || ...
          isempty(this.TX.Norm2MmrtUoutCalibrationGain) || ...
          isempty(this.TX.Norm2MmrtUout)
        return;
      end
      normCalibrated = this.MmrtUoutCalibratedAll(row,:) ./ ...
        this.TX.Norm2MmrtUoutCalibrationGain ./ this.TX.Norm2MmrtUout;
      this.NormCalibratedAll(row,:) = normCalibrated;
      if ~isempty(this.Norm)
        validIdx = isfinite(this.Norm) & this.NormCalibratedAll(row,:)==0;
        this.NormCalibratedAll(row,validIdx) = this.Norm(validIdx);
      end
    end


    function UpdateMmrtUoutCalibratedByNorm(this, row)
      %% step from normalized DAC to MmrtUout
      if size(this.NormCalibratedAll, 1) < row || ...
          isempty(this.TX.Norm2MmrtUoutCalibrationGain) || ...
          isempty(this.TX.Norm2MmrtUout)
        return;
      end
      mmrtUoutCalibrated = this.NormCalibratedAll(row,:) .* ...
        this.TX.Norm2MmrtUoutCalibrationGain .* this.TX.Norm2MmrtUout;
      this.MmrtUoutCalibratedAll(row,:) = mmrtUoutCalibrated;
      if ~isempty(this.Mmrt)
        validIdx = isfinite(this.Mmrt) & (this.MmrtUoutCalibratedAll(row,:)==0);
        this.MmrtUoutCalibratedAll(row,validIdx) = this.Mmrt(validIdx);
      end
    end

  end


  methods (Access = {?PD.TXClass, ?PD.HWClass})

    function UpdatePaUout2Amplitude(this)
      %% TX.PaUout2Amplitude has changed. Trigger recalculation
      for row = 1:4
        this.UpdateAmplitudeCalibratedByPaUout(row);
      end

      row = 5;
      this.UpdatePaUoutCalibratedByAmplitude(row);
      this.UpdateUoutCalibratedByPaUout(row);
      this.UpdateMmrtUoutCalibratedByUout(row);
      this.UpdateNormCalibratedByMmrtUout(row);
    end


    function UpdateUout2PaUout(this)
      %% TX.Uout2PaUout has changed. Trigger recalculation
      for row = 1:3
        this.UpdatePaUoutCalibratedByUout(row);
        this.UpdateAmplitudeCalibratedByPaUout(row);
      end

      for row = 4:5
        this.UpdateUoutCalibratedByPaUout(row);
        this.UpdateMmrtUoutCalibratedByUout(row);
        this.UpdateNormCalibratedByMmrtUout(row);
      end
    end


    function UpdateMmrtUout2Uout(this)
      %% TX.MmrtUout2Uout has changed. Trigger recalculation
      for row = 1:2
        this.UpdateUoutCalibratedByMmrtUout(row);
        this.UpdatePaUoutCalibratedByUout(row);
        this.UpdateAmplitudeCalibratedByPaUout(row);
      end

      for row = 3:5
        this.UpdateMmrtUoutCalibratedByUout(row);
        this.UpdateNormCalibratedByMmrtUout(row);
      end
    end


    function UpdateNorm2MmrtUout(this)
      %% TX.Norm2MmrtUout has changed. Trigger recalculation
      row = 1;
      this.UpdateMmrtUoutCalibratedByNorm(row);
      this.UpdateUoutCalibratedByMmrtUout(row);
      this.UpdatePaUoutCalibratedByUout(row);
      this.UpdateAmplitudeCalibratedByPaUout(row);

      for row = 2:5
        this.UpdateNormCalibratedByMmrtUout(row);
      end
    end

  end

end
