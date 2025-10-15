classdef MySave < PD.DynamicHandle
  %% Singleton object that stores persistent information about the device
  %
  % ----------------------------------------------------------------------------
  % (C) Copyright 2022 Pure Devices GmbH, Wuerzburg, Germany
  % www.pure-devices.com
  % ----------------------------------------------------------------------------


  properties
    % additional properties are added "on-the-fly"

    RootPath  % path to the root of the OpenMatlab installation
    UserPath  % path to the User folder (absolute or relative to RootPath)

    FPGA_Libs  % used driver library version
    libsLoaded  % .NET assemblies are loaded

    DummySerial  % if positive, "emulate" selected device serial
    reg  % quick access to some information stored in Windows registry

  end


  methods (Access = private)

    function this = MySave(varargin)
      %% (protected) constructor
      %
      %   this = MySave(InStruct)
      %
      % The input argument is optional. If it is a structure, it is converted to
      % a PD.MySave object. If it already is a PD.DynamicHandle, it is taken as
      % the (new) singleton reference.

      this = this@PD.DynamicHandle();
      if isempty(varargin)
        return;
      end
      if isstruct(varargin{1})
        this.SetStruct(varargin{1});
      elseif isa(varargin{1}, 'PD.DynamicHandle')
        this = varargin{1};
      elseif ~isempty(varargin{1})
        warning('PD.MySave: Input of class "%s" not supported.', class(varargin{1}))
      end

    end

  end


  methods (Static)

    function this = GetInstance(varargin)
      %% Get singleton instance of class

      % Use a persistent variable with "mlock" to store the singleton object to
      % avoid loosing it if the file is re-parsed.
      persistent pdMySave
      mlock();

      if isempty(pdMySave) || ~isvalid(pdMySave)
        pdMySave = PD.MySave(varargin{:});
      elseif nargin > 0 && isstruct(varargin{1})
        % override fields with values from input structure
        mySaveIn = varargin{1};
        pdMySave.SetStruct(mySaveIn);
      end

      this = pdMySave;
    end

  end


  methods

    function SetStruct(this, inStruct)
      % Update properties with fields in struct
      fnames = fieldnames(inStruct);

      s.type = '.';
      for iFields = 1:numel(fnames)
        s.subs = fnames{iFields};
        this = this.subsasgn(s, inStruct.(fnames{iFields}));
      end
    end

  end

end
