function [Rx, Ry, Rz] = get_aptDegRotationMatrix(alfa, phi, theta)
% get Rx Ry Rz form alfa phi theta in DEG
% The angles might be vectors (e.g. sequence_radial)

  n=max([numel(alfa) numel(phi) numel(theta)]);
  
  % Rx=[ 1         ,  0           ,  0; ...
  %      0         ,  cosd(alpha)  , -sind(alpha); ...
  %      0         ,  sind(alpha)  ,  cosd(alpha)];
  Rx = zeros(3,3,n);
  Rx(1,1,:) = 1;
  Rx(2,2,:) = cosd(alfa);
  Rx(3,3,:) = Rx(2,2,:);
  Rx(3,2,:) = sind(alfa);
  Rx(2,3,:) = -Rx(3,2,:);

  % Ry=[ cosd(phi)  ,  0          ,  sind(phi); ...
  %      0         ,  1          ,  0; ...
  %     -sind(phi)  ,  0          ,  cosd(phi)];
  Ry = zeros(3,3,n);
  Ry(1,1,:) = cosd(phi);
  Ry(2,2,:) = 1;
  Ry(3,3,:) = Ry(1,1,:);
  Ry(1,3,:) = sind(phi);
  Ry(3,1,:) = -Ry(1,3,:);

  % Rz=[ cosd(theta), -sind(theta), 0; ...
  %      sind(theta),  cosd(theta), 0; ...
  %      0         ,  0          , 1];
  Rz = zeros(3,3,n);
  Rz(1,1,:) = cosd(theta);
  Rz(2,2,:) = Rz(1,1,:);
  Rz(3,3,:) = 1;
  Rz(2,1,:) = sind(theta);
  Rz(1,2,:) = -Rz(2,1,:);
    
end