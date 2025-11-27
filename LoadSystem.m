function LoadSystem(varargin)
%% LOADSYSTEM - connect to MRT device and load variables
%
%   LoadSystem(doForce, useOop, useStruct, useDummy, userName)
%
% When called the first time, the connection to the MRT device is established
% and the necessary variables are loaded into the caller's (and base) workspace
% (imitating behavior of a script).
% On later calls, the structures "Seq", "TX", "AQ", and "Grad" are reset and the
% "HW" (and if applicable the "ShimMatrix" and the "DDS") objects are
% re-initialized.
% See optional input below to change standard behavior.
% On each call, the search paths are checked and set if necessary (see addpath),
% the connection to the device is (re-)established, and the files inside the
% "User" folder are loaded.
%
%
% INPUT:
%
%   doForce
%       optional, logical or string. If set to true, the connection to the MRT
%       device is re-opened and the HW (and ShimMatix) objects are re-created.
%       The strings 'true', 'force', 'reload' or 'init' are treated as true.
%       Everything else is treated as false.
%       Default: false.
%
%   useOop
%       optional, string. If set to 'oop', an object of class PD.HWClass is
%       created.
%       Default: not set.
%
%   useStruct
%       optional, string. If set to 'struct', a structure is created.
%       Default: not set.
%
%   useDummy
%       optional, string. If set to 'dummy', a dummy MRI device is loaded. If
%       the string has a numeric suffix (e.g., 'dummy123') the dummy MRI device
%       behaves as if a MRI device with that serial number were connected.
%       Default: not set.
%
%   userName
%       optional, string. Any string that doesn't match any of the above is used
%       as a user name.
%       Default: 'default' or the previously used user name.
%
%   All input parameters can be passed in arbitrary order but logical true (for
%   doForce) can only be passed as first argument.
%   The input arguments 'struct' and 'oop' are mutually exclusive. If none
%   are set, the type of "HW" is maintained as before the call to LoadSystem.
%   If "doForce" is set or "HW" is not present when calling LoadSystem, "HW" is
%   created as a PD.HWClass object.
%   If "userName" is 'default', %OpenMatlabRoot%\User (where %OpenMatlabRoot% is
%   the folder where OpenMatlab is installed) is used for settings and
%   calibration values. Otherwise, the folder %OpenMatlabRoot%\User\%UserName%
%   is used (where %UserName% is the value of userName).
%   The default values for the above input arguments can be overridden with a
%   script "LoadSystemDefaultInput.m" located at %OpenMatlabRoot%\User. In that
%   file, the default values for "doForce", "useOop", "useStruct", and
%   "useDummy" can be set as Boolean values; "userName" can be set to a valid
%   string. Input arguments of LoadSystem take presedence over any values that
%   are set in that file.
%
%
% OUTPUT:
%
%   none
%       However, the following variables are created (or potentially
%       overwritten) in the scope of the caller's context:
%           * HW
%           * mySave
%           * Seq
%           * TX
%           * AQ
%           * Grad
%           * (ShimMatrix)
%           * (DDS)
%       If the caller's context is different from the 'base' context, the
%       following variables are also assigned (and potentially overwritten) in
%       the 'base' context:
%           * HW
%           * mySave
%           * (ShimMatrix)
%           * (DDS)
%
%
% EXAMPLES:
%
%   LoadSystem
%     * Maintains type of HW. If HW is not present when calling LoadSystem,
%     HW is created as a PD.HWClass object.
%
%   LoadSystem oop
%     * Forces HW to be a PD.HWClass object.
%
%   LoadSystem struct
%     * Forces HW to be a structure.
%
%   LoadSystem true
%   LoadSystem force
%   LoadSystem reload
%   LoadSystem init
%     * Independent of whether HW already exists or of which type it is, HW
%       and all other respective variables are re-initialized or created in
%       the scope of the caller's workspace. HW is of type PD.HWClass.
%
%   LoadSystem oop true
%   LoadSystem oop force
%   LoadSystem oop reload
%   LoadSystem oop init
%   LoadSystem true oop
%   LoadSystem force oop
%   LoadSystem reload oop
%   LoadSystem init oop
%     * Independent of whether HW already exists or of which type it is, HW
%       and all other respective variables are re-initialized or created in
%       the scope of the caller's workspace. HW is of type PD.HWClass.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2025 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% un-documented input arguments
% NOTE: These arguments might change at any point in the future.
%
% "un-paired" arguments:
%   'createFullLoadMySystem'
%       If one of the input arguments is this string and the LoadMySystem file
%       for the selected user doesn't exist yet, create it with its full content
%       even when running in compiled code. This argument has no effect if not
%       running in compiled code.
%
% Property-value pairs:
%   regKeySuffix
%       The string that is following this input argument is used as a suffix to
%       the Windows registry keys "LibPath" and "InstallPath".
%       (Default: '')


%% Note:
% Instead of using input and output arguments for the following objects and
% structures, they are passed differently:
% evalin('caller', ...), evalin('base', ...) if the former is unsuccessful,
% assignin('caller', ...), and assignin('base', ...) are used for "HW",
% "mySave", "ShimMatrix", and "DDS" to make sure that these variables always
% have their correct names.
% assignin('caller', ...) is used for the structures "Seq", "TX", "AQ", and
% "Grad" for the same reason.


%% override default values for input with user settings
if isdeployed()
  % FIXME: What should we do on non-Windows?
  userRootPath = fullfile(getenv('localappdata'), 'Pure Devices', 'pdGUI', 'User');
  defInputScript = fullfile(userRootPath, 'LoadSystemDefaultInput.m');
  if exist(defInputScript, 'file')
    % Deployed executables cannot run plain-text .m files.
    % Read and eval content instead.
    defInputString = fileread(defInputScript);
    if ~isempty(defInputString)
      % eval content of script
      eval(defInputString);
    end
  end
else
  rootPath = fileparts(mfilename('fullpath'));
  defInputScript = fullfile(rootPath, 'User', 'LoadSystemDefaultInput.m');
  if exist(defInputScript, 'file')
    run(defInputScript);
  end
end


%% default input
if ~exist('doForce', 'var')
  doForce = false;
end
if ~exist('useOop', 'var')
  useOop = false;
end
if ~exist('useStruct', 'var')
  useStruct = false;
end
if ~exist('useDummy', 'var')
  useDummy = false;
end
if ~exist('userName', 'var')
  userName = '';
end

defaultUserName = true;
createFullLoadMySystem = false;
regKeySuffix = '';


%% get HW and mySave from context of caller (or base) if possible
HW = [];
mySave = [];
% prefer the objects in the caller's context
if evalin('caller', 'exist(''HW'', ''var'')')
  HW = evalin('caller', 'HW');
end
if evalin('caller', 'exist(''mySave'', ''var'')')
  mySave = evalin('caller', 'mySave');
end
% also try base context
if isempty(HW) && evalin('base', 'exist(''HW'', ''var'')')
  HW = evalin('base', 'HW');
end
if isempty(mySave) && evalin('base', 'exist(''mySave'', ''var'')')
  mySave = evalin('base', 'mySave');
end
% check if the types are usable
if ~isstruct(HW) && ~isa(HW, 'PD.HWClass')
  HW = [];
end
if ~isstruct(mySave) && ~isa(mySave, 'PD.MySave')
  mySave = [];
end


%% parse input
if ~isempty(varargin) && (islogical(varargin{1}) || isnumeric(varargin{1}))
  doForce = varargin{1};
  varargin(1) = [];
end

if ~isempty(varargin) && ischar(varargin{1})
  % check for property-value pairs first
  idxArg = find(~cellfun(@isempty, regexpi(varargin, '(regKeySuffix)', 'start')), 1, 'last');
  if ~isempty(idxArg) && idxArg < numel(varargin) && ischar(varargin{idxArg+1})
    regKeySuffix = varargin{idxArg+1};
    varargin([idxArg, idxArg+1]) = [];
  end

  % check for un-paired arguments
  isMatching = ~cellfun(@isempty, regexpi(varargin, '(true|force|reload|init)', 'start'));
  if any(isMatching)
    doForce = true;
    varargin(isMatching) = [];
  end

  isMatching = ~cellfun(@isempty, regexpi(varargin, '(oop)', 'start'));
  if any(isMatching)
    useOop = true;
    varargin(isMatching) = [];
  end

  isMatching = ~cellfun(@isempty, regexpi(varargin, '(struct)', 'start'));
  if any(isMatching)
    useStruct = true;
    varargin(isMatching) = [];
  end

  isMatching = ~cellfun(@isempty, regexpi(varargin, '(dummy)', 'start'));
  if any(isMatching)
    useDummy = true;
    serialString = regexpi(varargin{find(isMatching, 1, 'last')}, 'dummy([0-9]*)', 'tokens');
    if ~isempty(serialString) && ~isempty(serialString{1})
      dummySerial = str2double(serialString{1}{1});
      if isfinite(dummySerial)
        mySave.DummySerial = dummySerial;
      end
    end
    varargin(isMatching) = [];
  end

  isMatching = ~cellfun(@isempty, regexpi(varargin, '(createFullLoadMySystem)', 'start'));
  if any(isMatching)
    createFullLoadMySystem = true;
    varargin(isMatching) = [];
  end
end

if useStruct && useOop
  error('OpenMatlab:LoadSystem:type', ...
    '"struct" and "oop" are exclusive arguments.');
end

if numel(varargin) > 1
  error('OpenMatlab:LoadSystem:oneUserName', ...
    'Only one user name can be specified at a time.');
elseif numel(varargin) > 0
  if ischar(varargin{1})
    userName = varargin{1};
    defaultUserName = false;
  else
    error('OpenMatlab:LoadSystem:wrongUserName', ...
      '"userName" must be a string');
  end
end

if useDummy && (~isfield(mySave, 'DummySerial') || isempty(mySave.DummySerial))
  mySave.DummySerial = 1;
end


%% get optionally loaded objects from context of caller (or base) if possible
keepShimMatrix = false;
keepDDS = false;
if doForce
  evalin('caller', 'clear ShimMatrix');
  evalin('base', 'clear ShimMatrix');
  if evalin('caller', 'exist(''DDS'', ''var'')') && ...
      evalin('caller', 'isvalid(''DDS'')')
    evalin('caller', 'delete(DDS)');
  end
  evalin('caller', 'clear DDS');
  if evalin('base', 'exist(''DDS'', ''var'')') && ...
      evalin('base', 'isvalid(''DDS'')')
    evalin('base', 'delete(DDS)');
  end
  evalin('base', 'clear DDS');
else
  if evalin('caller', 'exist(''ShimMatrix'', ''var'')')
    ShimMatrix = evalin('caller', 'ShimMatrix');
    keepShimMatrix = true;
  end
  if ~exist('ShimMatrix', 'var') && ...
      evalin('base', 'exist(''ShimMatrix'', ''var'')')
    ShimMatrix = evalin('base', 'ShimMatrix');
    keepShimMatrix = true;
  end
  if evalin('caller', 'exist(''DDS'', ''var'')')
    DDS = evalin('caller', 'DDS');
    keepDDS = true;
  end
  if ~exist('DDS', 'var') && ...
      evalin('base', 'exist(''DDS'', ''var'')')
    DDS = evalin('base', 'DDS');
    keepDDS = true;
  end
end


%% add path to LoadHW
if isdeployed()
  mySave.RootPath = getOpenMatlabRootPath();
else
  if exist(fullfile(fileparts(mfilename('fullpath')), 'StandardFunctions', 'HWFunctions', 'getOpenMatlabRootPath.m'), 'file')
    % function was found in folder relative to current folder
    movePathToFront(fullfile(fileparts(mfilename('fullpath')), 'StandardFunctions', 'HWFunctions'));
    mySave.RootPath = getOpenMatlabRootPath();
  elseif ~isempty(which('getOpenMatlabRootPath'))
    mySave.RootPath = getOpenMatlabRootPath();
    movePathToFront(fullfile(mySave.RootPath, 'StandardFunctions', 'HWFunctions'));
  else
    warning('OpenMatlab:LoadSystem:Path', ...
      ['Change your current directory to the root folder of the ', ...
      'OpenMatlab toolbox before calling function "%s" the first time.\n', ...
      'Trying to load anyway. Might fail.'], mfilename)
    addpath('StandardFunctions/HWFunctions');
    mySave.RootPath = fileparts([mfilename('fullpath'),'.m']);
  end
end


%% check user
if defaultUserName
  if isempty(userName)
    userName = 'default';
  end
  if ((isa(HW, 'PD.HWClass') && isvalid(HW)) || isfield(HW, 'UserName')) && ~isempty(HW.UserName)
    userName = HW.UserName;
  end
end


%% copy UserPath from HW if available
if ((isa(HW, 'PD.HWClass') && isvalid(HW)) || isfield(HW, 'UserPath')) && ~isempty(HW.UserPath)
  mySave.UserPath = HW.UserPath;
end


%% decide which type
if ~doForce && ~useStruct && ~useOop && isa(HW, 'struct')
  % keep as struct
  useOop = false;
elseif ~useStruct
  useOop = true;
end


%% create HW and mySave (and ShimMatrix if applicable)
% delete HW before creating new one
if (doForce || useStruct) && isa(HW, 'PD.HWClass')
  delete(HW);
end

if keepShimMatrix
  ShMatPV = {'ShimMatrix', ShimMatrix};
else
  ShMatPV = {};
end
if keepDDS
  DDSPV = {'DDS', DDS};
else
  DDSPV = {};
end
[HW, mySave, ShimMatrix, DDS] = LoadHW('HW', HW, 'mySave', mySave, ...
          ShMatPV{:}, DDSPV{:}, ...
          'doForce', doForce, 'oop', useOop, 'user', userName, ...
          'usingDefaultUserName', defaultUserName, ...
          'createFullLoadMySystem', createFullLoadMySystem, ...
          'regKeySuffix', regKeySuffix);
evalin('caller', 'clear AQSlice SliceSelect');  % moved here from LoadHW

% return early in case LoadHW exitted without loading the device
if isempty(HW) || (isfield(mySave, 'isUpdating') && mySave.isUpdating)
  return;
end

% assign variables in context of caller (and base)
% FIXME: 'base' might be the same context as 'caller'. Should we try to avoid a
% possible duplicate assigment here?
assignin('caller', 'HW', HW);
assignin('base', 'HW', HW);
assignin('caller', 'mySave', mySave);
assignin('base', 'mySave', mySave);
if ~keepShimMatrix && ~isempty(ShimMatrix)
  assignin('caller', 'ShimMatrix', ShimMatrix);
  assignin('base', 'ShimMatrix', ShimMatrix);
end
if ~keepDDS && ~isempty(DDS)
  assignin('caller', 'DDS', DDS);
  assignin('base', 'DDS', DDS);
end


%% reset structures with parameters for demos and standard sequences
ResetStructs;


%% assign variables in context of caller
assignin('caller', 'Seq', Seq(:));  %#ok<NODEF>
assignin('caller', 'TX', TX(:));  %#ok<NODEF>
assignin('caller', 'AQ', AQ(:));  %#ok<NODEF>
assignin('caller', 'Grad', Grad(:));  %#ok<NODEF>

end


function movePathToFront(path1)
%% add path if it is not on front of search paths
mp = matlabpath;
if ~strncmpi(mp, [path1 pathsep], length(path1)+1)
  addpath(path1);
end

end
