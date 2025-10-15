function res = sphericalConnection(startPoint, endPoint, props)
%% Find points on "spherical" trajectory that connects two points
%
%   res = sphericalConnection(startPoint, endPoint, props)
%
% The trajectory for the connection between the two points is the shortest
% connection on a sphere if both start and end points have the same amplitude.
% If the amplitude of the two points differ, the trajectory is chosen such the
% amplitude changes uniformly.
%
%
% INPUT:
%
%   startPoint
%       3 element vector with starting point of the trajectory.
%
%   endPoint
%       3 element vector with ending point of the trajectory.
%
%   props
%       structure with properties of the trajectory. The following fields can be
%       used. If they are omitted or empty, default values are used:
%
%     nPoints
%         Number of points on the resulting trajectory (default: 21).
%
%
% OUTPUT:
%
%   res
%       Nx3 matrix with the points on the trajectory.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
if nargin < 3 || isempty(props), props = struct(); end

if isemptyfield(props, 'nPoints'),  props.nPoints = 21;  end

normStart = norm(startPoint);
normEnd = norm(endPoint);


%% special cases if one or both of the end points are 0
if normStart < eps && normEnd < eps
  res = zeros(props.nPoints, 3);
  return;
end

if normStart/normEnd < 1e-6
  amps = linspace(normStart, normEnd, props.nPoints).';
  direction = endPoint(:).' - startPoint(:).';
  res = bsxfun(@times, direction/norm(direction), amps);
  return;
end

if normEnd/normStart < 1e-6
  amps = linspace(normStart, normEnd, props.nPoints).';
  direction = startPoint(:).' - endPoint(:).';
  res = bsxfun(@times, direction/norm(direction), amps);
  return;
end

%% rotate such that startPoint is in x-y plane
rot1_axis = cross(startPoint, [0, 0, 1]);

if norm(rot1_axis)/norm(startPoint) < 1/1000
  % startPoint is approximately parallel to z-axis
  rot1_axis = [1, 0, 0];  % arbitrarily chose x-axis for first rotation
else
  rot1_axis = rot1_axis / norm(rot1_axis);
end

rot1_angle = asin(startPoint(3) / normStart);

R1 = rot_mat(rot1_axis, rot1_angle);

startPoint1 = startPoint(:).' * R1;
endPoint1 = endPoint(:).' * R1;


%% rotate around startPoint1 such that endPoint1 is in x-y plane
rot2_axis = startPoint1 / norm(startPoint1);
rot2_axis(3) = 0;

n2 = cross(rot2_axis, endPoint1);
n2 = n2 / norm(n2);
rot2_angle = acos(([0, 0, 1] * n2.'));

R2 = rot_mat(rot2_axis, rot2_angle);

endPoint2 = endPoint1(:).' * R2;

if abs(endPoint2(3)) > 10*eps
  % FIXME: What is the "real" condition for the other cosine branch?
  rot2_angle = -acos(([0, 0, 1] * n2.'));
  R2 = rot_mat(rot2_axis, rot2_angle);

  endPoint2 = endPoint1(:).' * R2;
end
% startPoint2 = startPoint1(:).' * R2;  % shouldn't change anything


%% interpolate trajectory in x-y plane
% angle
startAngle = atan2(startPoint1(2), startPoint1(1));
endAngle = atan2(endPoint2(2), endPoint2(1));

angles = unwrap([startAngle, endAngle]);

angles = linspace(angles(1), angles(end), props.nPoints).';

% amp
amps = linspace(normStart, normEnd, props.nPoints).';

res2 = bsxfun(@times, [cos(angles), sin(angles), zeros(size(angles))], amps);


%% rotate back to initial coordinate system
res = res2 * (R2.' * R1.');


end


function res = cross(a, b)
%% cross product of vectors a and b
res = [a(2)*b(3)-a(3)*b(2), a(3)*b(1)-a(1)*b(3), a(1)*b(2)-a(2)*b(1)];

end


function R = rot_mat(rot_axis, rot_angle)
%% rotation matrix for rotation by angle around axis

cosA = cos(rot_angle);
sinA = sin(rot_angle);

R = cosA * eye(3) + ...
  sinA * [0, -rot_axis(3), rot_axis(2); rot_axis(3), 0, -rot_axis(1); -rot_axis(2), rot_axis(1), 0] + ...
  (1-cosA) * rot_axis(:)*rot_axis(:).';

end
