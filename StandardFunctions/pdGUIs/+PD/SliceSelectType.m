classdef SliceSelectType < uint32
  %% Enumeration class for identification of sequence in Slice Select GUI
  %
  % ----------------------------------------------------------------------------
  % (C) Copyright 2018-2024 Pure Devices GmbH, Wuerzburg, Germany
  % www.pure-devices.com
  % ----------------------------------------------------------------------------

  enumeration
    % The values in brackets correspond to the position in the popupmenu.
    % They must be consecutive starting at 1.
    SliceSelect (3)
    SpinEcho (1)
    GradEcho (2)
  end

  methods
    function str = char(this)
      %% Convert enumeration to string
      % These are the strings that are displayed in the popupmenu.
      str = cell(size(this));

      for iElem = 1:numel(this)
        switch this(iElem)
          case PD.SliceSelectType.SliceSelect
            str{iElem} = 'Slice Selection';
          case PD.SliceSelectType.SpinEcho
            str{iElem} = 'Spin Echo';
          case PD.SliceSelectType.GradEcho
            str{iElem} = 'Grad Echo';
          otherwise
            str{iElem} = 'undefined name';
        end
      end

      if isscalar(this), str = str{1}; end
    end

  end

end
