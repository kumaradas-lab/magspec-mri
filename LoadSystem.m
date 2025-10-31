function LoadSystem(varargin)
%% LOADSYSTEM - connect to MRT device and load variables
%
%   LoadSystem(doForce, useOop, useStruct, useDummy, userName)
%
% When called the first time, the connection to the MRT device is
% established and the necessary variables are loaded into the caller's
% workspace (imitating behavior of a script).
% On later calls, the structures Seq, TX, AQ and Grad are reset and the HW
% (and ShimMatrix) objects are re-initialized. See optional input below to
% change standard behavior.
% On each call, the search paths are checked and set if necessary, the
% connection to the device is (re-)established and the files inside the
% "User" folder are loaded.
%
%
% INPUT:
%
%   doForce
%         optional, logical or string. If set to true, the connection to the MRT
%         device is re-opened and the HW (and ShimMatix) objects are re-created.
%         The strings 'true', 'force', 'reload' or 'init' are treated as true.
%         Everything else is treated as false.
%         Default: false.
%
%   useOop
%         optional, string. If set to 'oop', an object of class PD.HW is
%         created.
%         Default: not set.
%
%   useStruct
%         optional, string. If set to 'struct', a structure is created.
%         Default: not set.
%
%   useDummy
%         optional, string. If set to 'dummy', a dummy MRI device is loaded.
%         Default: not set.
%
%   userName
%         optional, string. Any string that doesn't match any of the above is
%         used as a user name.
%         Default: 'default' or the previously used user name.
%
%   All input parameters can be passed in arbitrary order but logical true (for
%   doForce) can only be passed as first argument.
%   The input arguments 'struct' and 'oop' are mutually exclusive. If none
%   are set, the type of HW maintained as before the call to LoadSystem.
%   If "doForce" is set or HW is not present when calling LoadSystem, HW is
%   created as a PD.HW object.
%   If "userName" is 'default', %OpenMatlabRoot%\User (where %OpenMatlabRoot% is
%   the folder where OpenMatlab is installed) is used for settings and
%   calibration values. Otherwise, the folder %OpenMatlabRoot%\User\%UserName%
%   is used (where %UserName% is the value of userName).
%   The default values for the above input arguments can be overridden with a
%   script "LoadSystemDefaultInput.m" located at %OpenMatlabRoot%\User. In that
%   file, the default values for "doForce", "useOop", "useStruct", and
%   "useDummy" can be set as Boolean values; "userName" can be set to a valid
%   string. This only applies if LoadSystem is called without input arguments.
%
%
% OUTPUT:
%
%   none
%         However, the following variables are created (or potentially
%         overwritten) in the scope of the caller's workspace:
%           * HW
%           * mySave
%           * Seq
%           * TX
%           * AQ
%           * Grad
%           * (ShimMatrix)
%
%
% EXAMPLES:
%
%   LoadSystem
%     * Maintains type of HW. If HW is not present when calling LoadSystem,
%     HW is created as a PD.HW object.
%
%   LoadSystem oop
%     * Forces HW to be a PD.HW object.
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
%       the scope of the caller's workspace. HW is of type PD.HW.
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
%       the scope of the caller's workspace. HW is of type PD.HW.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2011-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% Note:
% evalin('caller', ...) and assignin('caller', ...) are used for "HW", "mySave"
% and "ShimMatrix" to make sure that these variables always have their correct
% names.

%% default input
doForce = false;
useOop = false;
useStruct = false;
useDummy = false;
userName = '';
defaultUserName = true;

%% override default values for input with user settings
if nargin == 0
  rootPath = fileparts(mfilename('fullpath'));
  defInputScript = fullfile(rootPath, 'User', 'LoadSystemDefaultInput.m');
  if exist(defInputScript, 'file')
    run(defInputScript);
  end
end

%% parse input
if ~isempty(varargin) && (islogical(varargin{1}) || isnumeric(varargin{1}))
  doForce = varargin{1};
  varargin(1) = [];
end

if ~isempty(varargin) && ischar(varargin{1})
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

%% get HW and mySave from context of caller if possible
HW = [];
mySave = [];
if evalin('caller', 'exist(''HW'', ''var'')')
  HW = evalin('caller', 'HW');
end
if evalin('caller', 'exist(''mySave'', ''var'')')
  mySave = evalin('caller', 'mySave');
end
if useDummy
  mySave.DummySerial = 1; % set it to 1 or [] depending on if real MRI is connected or not
end
keepShimMatrix = false;
keepDDS = false;
if doForce
  evalin('caller', 'clear ShimMatrix');
  if evalin('caller', 'exist(''DDS'', ''var'')')
    evalin('caller', 'delete(DDS)');
  end
  evalin('caller', 'clear DDS');
else
  if evalin('caller', 'exist(''ShimMatrix'', ''var'')')
    ShimMatrix = evalin('caller', 'ShimMatrix');
    keepShimMatrix = true;
  end
  if evalin('caller', 'exist(''DDS'', ''var'')')
    DDS = evalin('caller', 'DDS');
    keepDDS = true;
  end
end

%% add path to LoadHW
if isdeployed()
  mySave.RootPath = getOpenMatlabRootPath();
else
  if exist(fullfile(pwd, 'StandardFunctions', 'HWFunctions', 'getOpenMatlabRootPath.m'), 'file')
    % function was found in folder relative to current folder
    movePathToFront(fullfile(pwd, 'StandardFunctions', 'HWFunctions'));
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
  if ((isa(HW, 'PD.HW') && isvalid(HW)) || isfield(HW, 'UserName')) && ~isempty(HW.UserName)
    userName = HW.UserName;
  end
end

%% copy UserPath from HW if available
if ((isa(HW, 'PD.HW') && isvalid(HW)) || isfield(HW, 'UserPath')) && ~isempty(HW.UserPath)
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
if (doForce || useStruct) && isa(HW, 'PD.HW')
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
          'doForce', doForce, 'oop', useOop, 'user', userName);
evalin('caller', 'clear AQSlice SliceSelect');  % moved here from LoadHW

% return early in case LoadHW exitted without loading the device
if isempty(HW) || (isfield(mySave, 'isUpdating') && mySave.isUpdating)
  return;
end

% assign variables in context of caller
assignin('caller', 'HW', HW);
assignin('caller', 'mySave', mySave);
if ~keepShimMatrix && ~isempty(ShimMatrix)
  assignin('caller', 'ShimMatrix', ShimMatrix);
end
if ~keepDDS && ~isempty(DDS)
  assignin('caller', 'DDS', DDS);
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
