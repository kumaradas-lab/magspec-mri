function hAxes = plot_data_1D_TR(HW, data_1D, hAxes, raiseWindow, config)
%% Plot RF signal (abs, real, imag), its phase and its frequency offset to AQ frequency
%
%       hAxes = plot_data_1D_TR(HW, data_1D, hAxes, raiseWindow, config)
%
% The "tRep"s are stacked upon each other such that each "tRep" starts at t=0s.
%
% INPUT:
%   HW:           HW structure or object
%   data_1D:      structure with measurement data
%   hAxes:        array with graphics handles for the three axes containing the
%                 RF signal (abs, real, imag), its phase and its frequency
%                 offset to AQ frequency (see corresponding output argument).
%                 Or: graphics handle to a valid parent for these three axes
%                 (uipanel or figure handle).
%                 If omitted, empty or 1, figure 99 is used.
%   raiseWindow:  boolean. If false and hAxes is empty, the figure window does
%                 not steal the focus if it already exists. (Default: true)
%
% OUTPUT:
%   hAxes:        array with graphics handles for the three axes containing the
%                 RF signal (abs, real, imag), its phase and its frequency
%                 offset to AQ frequency (see corresponding input argument).
%
% ------------------------------------------------------------------------
% (C) Copyright 2011-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------

%% Input check
if nargin < 3, hAxes = []; end
if nargin < 4, raiseWindow = true; end
if nargin < 5, config = []; end

config = set_EmptyField(config, 'defaultFigure', 99);
config = set_EmptyField(config, 'figureTitle', 'Timeline wrapped at TR');
config = set_EmptyField(config, 'timeFieldname', 'time_of_tRep');

hAxes = plot_data_1D(HW, data_1D, hAxes, raiseWindow, config);

end
