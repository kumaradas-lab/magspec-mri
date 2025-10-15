function [XX_magnet, YY_magnet, ZZ_magnet] = RotateImage2Magnet(AQSlice)
%% rotate image coordinate system to magnet coordinate system

if sum([AQSlice.nPhase(:)>1; AQSlice.nRead>1]) < 3
    error('only for 3d images')
end

% grid in image coordinate system
if AQSlice.nRead>1
    x = AQSlice.CenterRot(1)+linspace(-AQSlice.sizeRead/2,     AQSlice.sizeRead/2,     AQSlice.nRead);
    y = AQSlice.CenterRot(2)+linspace(-AQSlice.sizePhase(1)/2, AQSlice.sizePhase(1)/2, AQSlice.nPhase(1));
    z = AQSlice.CenterRot(3)+linspace(-AQSlice.sizePhase(2)/2, AQSlice.sizePhase(2)/2, AQSlice.nPhase(2));
else
    x = AQSlice.CenterRot(1)+linspace(-AQSlice.sizePhase(2)/2, AQSlice.sizePhase(2)/2, AQSlice.nPhase(2));
    y = AQSlice.CenterRot(2)+linspace(-AQSlice.sizePhase(1)/2, AQSlice.sizePhase(1)/2, AQSlice.nPhase(1));
    z = AQSlice.CenterRot(3)+linspace(-AQSlice.sizePhase(3)/2, AQSlice.sizePhase(3)/2, AQSlice.nPhase(3));
end

[YY_image, XX_image, ZZ_image] = meshgrid (y, x, z);

Rx = [ 1                   , 0                     , 0; ...
       0                   , cos(AQSlice.alfa)     ,-sin(AQSlice.alfa); ...
       0                   , sin(AQSlice.alfa)     , cos(AQSlice.alfa)];

Ry = [ cos(AQSlice.phi)    , 0                     , sin(AQSlice.phi); ...
       0                   , 1                     , 0; ...
      -sin(AQSlice.phi)    , 0                     , cos(AQSlice.phi)];

Rz = [ cos(AQSlice.theta)  ,-sin(AQSlice.theta)    , 0; ...
       sin(AQSlice.theta)  , cos(AQSlice.theta)    , 0; ...
       0                   , 0                     , 1];


coor_magnet = (Rz*Ry*Rx*[YY_image(:),ZZ_image(:),XX_image(:)]')';

XX_magnet = reshape(coor_magnet(:,1), size(XX_image));
YY_magnet = reshape(coor_magnet(:,2), size(XX_image));
ZZ_magnet = reshape(coor_magnet(:,3), size(XX_image));

end
