function [varargout] = Find_Shim_Readjust(varargin)
%% Re-adjust an already good magnet shim
%
%   [HW, mySave, SliceSelectOut] = Find_Shim_Readjust(HW, mySave, minTime, doPlot, Seq, SliceSelect)
%
% This function calls "Find_Shim" with settings that allow for a quick
% re-shimming if the magnet shim is already close to the optimum values.
%
% This function is typically called with:
%
%   [HW, mySave] = Find_Shim_Readjust(HW, mySave, 0);
%
%
% INPUT:
%
% The input is the same as for the function "Find_Shim" but some default
% properties differ:
%
%   minTime
%         (default: 0)
%
%   Seq
%         The following fields in the structure Seq have different default
%         values from the ones in "Find_Shim":
%     iterations
%           (default: 50)
%     ShimStart
%           (default: HW.MagnetShim)
%     ShimStep
%           (default: [0.5e-3, 0.5e-3, 0.5e-3, 0.1e-6])
%
%
% OUTPUT:
%
% See documentation for "Find_Shim".
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2020-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------

%% default input parameters
if nargin < 1  % HW
  error('PD:Find_Shim_Readjust:noHW', 'The first input argument "HW" is mandatory.');
end
if nargin < 6  % SliceSelect
  [varargin{(nargin+1):6}] = deal([]);
end

HW = varargin{1};
minTime = varargin{3};
Seq = varargin{5};
SliceSelect = varargin{6};

if isempty(minTime), minTime = 0; end
if isemptyfield(Seq, 'iterations'), Seq.iterations = 50; end
% FIXME: Optionally use all available gradient channels for shimming
if isemptyfield(SliceSelect, 'iDevice'), SliceSelect.iDevice = 1; end
if isemptyfield(Seq, 'ShimStart'), Seq.ShimStart = HW.Grad(SliceSelect.iDevice).AmpOffset; end
if isemptyfield(Seq, 'ShimStep')
  Seq.ShimStep = [0.5e-3, 0.5e-3, 0.5e-3, repmat(0.1e-6, 1, numel(HW.Grad(SliceSelect.iDevice).B))];
end

varargin{3} = minTime;
varargin{5} = Seq;

%% pass operation on to "Find_Shim"
[varargout{1:nargout}] = Find_Shim(varargin{:});

end
