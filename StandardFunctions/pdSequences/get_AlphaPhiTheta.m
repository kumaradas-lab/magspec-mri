function AQSlice = get_AlphaPhiTheta(AQSlice, orientation_SPR)
%% Get Euler angles for cardinal directions of image coordinate system
%
%   AQSlice = get_AlphaPhiTheta(AQSlice, orientation_SPR)
%
% This function helps to get the Euler angles to orient an image along the
% cardinal axes of the magnet coordinate system.
%
% INPUT:
%
%   AQSlice
%         AQSlice structure as used by the imaging sequences
%         "sequence_Spin_Echo" and "sequence_Flash".
%
%   orientation_SPR
%         Character vector describing the order of imaging encoding directions
%         (slice/phase(1), phase(2), read/phase(3)) in the magnet coordinate
%         system.
%         For 3d images/volumes, the following (right-handed coordinate system)
%         orientations are supported:
%              'xyz',    'yzx',    'zxy'
%             '-xzy',   'zy-x',   'y-xz'
%             '-yxz',   'xz-y',   'z-yx'
%             '-zyx',   'yx-z',   'x-zy'
%            'x-y-z',  '-y-zx',  '-zx-y'
%            'y-z-x',  '-z-xy',  '-xy-z'
%            'z-x-y',  '-x-yz',  '-yz-x'
%           '-z-y-x', '-y-x-z', '-x-z-y'
%         For 2d images with (phase(2), phase(3)/read) encoding, the following
%         orientations are supported:
%             'xy',   'xz',   'yx',   'yz',   'zx',   'zy'
%           '-x-y', '-x-z', '-y-x', '-y-z', '-z-x', '-z-y'
%            'x-y',  'x-z',  'y-x',  'y-z',  'z-x',  'z-y'
%            '-xy',  '-xz',  '-yx',  '-yz',  '-zx',  '-zy'
%         For 1d images/profiles with (phase(3)/read) encoding, the following
%         orientations are supported:
%           'x', 'y', 'z', '-x', '-y', '-z'
%
%
% OUTPUT:
%
%   AQSlice
%         Same as the input "AQSlice" but with the following additional field
%         with values that result in an image orientation as defined by
%         "orientation_SPR":
%
%     alfa
%           First Euler angle "alfa" (alpha) with the rotation around the x axis
%           in radians.
%
%     phi
%           Second Euler angle "phi" with the rotation around the rotated y axis
%           in radians.
%
%     theta
%           Third Euler angle "theta" with the rotation around the rotated z
%           axis in radians.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2024 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
%-------------------------------------------------------------------------------


%% default input
if nargin < 2
    orientation_SPR = 'xyz';
end


%% map input to Euler angles
switch lower(orientation_SPR)
  case {'xyz', 'yz', 'z'}
    AQSlice(1).alfa = 0.0*pi;
    AQSlice(1).phi  = 0.0*pi;
    AQSlice(1).theta= 0.0*pi;

  case {'yzx', 'zx', 'x'}
    AQSlice(1).alfa = 0.5*pi;
    AQSlice(1).phi  = 0.0*pi;
    AQSlice(1).theta= 0.5*pi;

  case {'zxy', 'xy', 'y'}
    AQSlice(1).alfa =-0.5*pi;
    AQSlice(1).phi  =-0.5*pi;
    AQSlice(1).theta= 0.0*pi;

  case {'-xzy', 'zy'}
    AQSlice(1).alfa =-0.5*pi;
    AQSlice(1).phi  = 1.0*pi;
    AQSlice(1).theta= 0.0*pi;

  case {'zy-x', 'y-x'}
    AQSlice(1).alfa = 1.0*pi;
    AQSlice(1).phi  =-0.5*pi;
    AQSlice(1).theta= 1.0*pi;

  case {'y-xz', '-xz'}
    AQSlice(1).alfa = 1.0*pi;
    AQSlice(1).phi  = 1.0*pi;
    AQSlice(1).theta=-0.5*pi;

  case {'-yxz', 'xz'}
    AQSlice(1).alfa = 0.0*pi;
    AQSlice(1).phi  = 0.0*pi;
    AQSlice(1).theta=-0.5*pi;

  case {'xz-y', 'z-y'}
    AQSlice(1).alfa =-0.5*pi;
    AQSlice(1).phi  = 1.0*pi;
    AQSlice(1).theta= 1.0*pi;

  case {'z-yx', '-yx'}
    AQSlice(1).alfa = 0.0*pi;
    AQSlice(1).phi  =-0.5*pi;
    AQSlice(1).theta= 1.0*pi;

  case {'-zyx', 'yx'}
    AQSlice(1).alfa = 0.0*pi;
    AQSlice(1).phi  = 0.5*pi;
    AQSlice(1).theta= 0.0*pi;

  case {'yx-z', 'x-z'}
    AQSlice(1).alfa = 1.0*pi;
    AQSlice(1).phi  = 0.0*pi;
    AQSlice(1).theta= 0.5*pi;

  case {'x-zy', '-zy'}
    AQSlice(1).alfa = 0.5*pi;
    AQSlice(1).phi  = 1.0*pi;
    AQSlice(1).theta= 1.0*pi;

  case {'x-y-z','-y-z','-z'}
    AQSlice(1).alfa = 1.0*pi;
    AQSlice(1).phi  = 0.0*pi;
    AQSlice(1).theta= 0.0*pi;

  case {'-y-zx', '-zx'}
    AQSlice(1).alfa =-0.5*pi;
    AQSlice(1).phi  = 0.0*pi;
    AQSlice(1).theta=-0.5*pi;

  case {'-zx-y', 'x-y'}
    AQSlice(1).alfa = 0.5*pi;
    AQSlice(1).phi  = 0.5*pi;
    AQSlice(1).theta= 0.0*pi;

  case {'y-z-x','-z-x','-x'}
    AQSlice(1).alfa =-0.5*pi;
    AQSlice(1).phi  = 0.0*pi;
    AQSlice(1).theta= 0.5*pi;

  case {'-z-xy', '-xy'}
    AQSlice(1).alfa = 0.0*pi;
    AQSlice(1).phi  = 0.5*pi;
    AQSlice(1).theta= 0.5*pi;

  case {'-xy-z', 'y-z'}
    AQSlice(1).alfa = 0.0*pi;
    AQSlice(1).phi  = 1.0*pi;
    AQSlice(1).theta= 0.0*pi;

  case {'z-x-y', '-x-y', '-y'}
    AQSlice(1).alfa = 0.0*pi;
    AQSlice(1).phi  =-0.5*pi;
    AQSlice(1).theta= 0.5*pi;

  case {'-x-yz', '-yz'}
    AQSlice(1).alfa = 0.0*pi;
    AQSlice(1).phi  = 0.0*pi;
    AQSlice(1).theta= 1.0*pi;

  case {'-yz-x', 'z-x'}
    AQSlice(1).alfa = 0.5*pi;
    AQSlice(1).phi  = 0.0*pi;
    AQSlice(1).theta=-0.5*pi;

  case {'-z-y-x', '-y-x'}
    AQSlice(1).alfa = 1.0*pi;
    AQSlice(1).phi  = 0.5*pi;
    AQSlice(1).theta= 0.0*pi;

  case {'-y-x-z', '-x-z'}
    AQSlice(1).alfa = 1.0*pi;
    AQSlice(1).phi  = 0.0*pi;
    AQSlice(1).theta=-0.5*pi;

  case {'-x-z-y', '-z-y'}
    AQSlice(1).alfa = 0.5*pi;
    AQSlice(1).phi  = 1.0*pi;
    AQSlice(1).theta= 0.0*pi;

  otherwise
    if ~ischar(orientation_SPR)
      orientation_SPR = num2str(orientation_SPR);
    end
    error('PD:get_AlphaPhiTheta:UnsupportedOrientation', ...
      ['Unsupported orientation "%s". ', ...
      'See "<a href="matlab:doc get_AlphaPhiTheta">doc get_AlphaPhiTheta</a>" ', ...
      'for a list of supported coordinate systems.'], ...
      orientation_SPR);

end
