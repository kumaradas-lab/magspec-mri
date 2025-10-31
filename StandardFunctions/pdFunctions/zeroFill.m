function filled = zeroFill(toFill, outputSize, winsize)
%% Apply zero-filling and filter to k-space
%
%   filled = zeroFill(toFill, outputSize, winsize)
%
% Zero filling is used to essentially enhance the resolution of an acquired
% image. The complex data allows to locate borders of a measured structure on a
% sub-pixel level. For this, the k-space is expanded to a higher frequency range
% and filled with zeros before calculating the image via FFT.
% Additionally, this function allows to smooth out images by reducing the
% amplitude of high frequency portions (which is likely due to noise).
%
%
% INPUT:
%
%   toFill
%         complex double multi-dimensional array containing the k-space data
%         (frequency 0 at the center).
%   outputSize
%         vector containing the dimensions of the zero filled output k-space.
%   winsize
%         scalar or vector containing the relative size of the k-space filter or
%         structure with the following fields:
%     size
%           scalar or vector containing the relative size of the k-space filter.
%           The highest k-space frequency in each dimension corresponds to 1.
%           (no default, mandatory)
%     winexp
%           scalar with the exponent used in the filter function (default: 1)
%     used
%           Boolean value indicating whether a filter should be applied
%           (default: true if nargin >=3; false otherwise)
%     winType
%           string indicating the type of the k-space filter function (default:
%           'RaisedCos'). All filter functions are symmetric with respect to the
%           center of the k-space and only depend on the kartesian distance (R)
%           from the k-space center. The multiplicative filter is applied before
%           zero-filling.
%           The following strings are supported:
%       'RaisedCos':
%             Raised cosine function (cos(R*pi).^(2.*winsize.winexp) where the
%             relative kartesian distance to the k-space center is smaller
%             than 0.5 (see winsize.size) and 0 otherwise.
%       'Gaussian':
%             Gaussian bell curve:
%                 exp(- ((R*exp(1/2)*2).^2 .* winsize.winexp)
%
%
% OUTPUT:
%
%   filled
%         The zero-filled and filtered k-space with the same number of
%         dimensions as input toFill and size according to outputSize.
%
%
% ------------------------------------------------------------------------------
% (C) Copyright 2014-2021 Pure Devices GmbH, Wuerzburg, Germany
% www.pure-devices.com
% ------------------------------------------------------------------------------


%% default input
if nargin < 3
  winsize.used = 0;
end
if ~isa(winsize, 'struct')
  t = winsize;
  clear winsize
  winsize.size = t;
else
  if isfield(winsize, 'i1'), winsize.size(1) = winsize.i1; end
  if isfield(winsize, 'i2'), winsize.size(2) = winsize.i2; end
  if isfield(winsize, 'i3'), winsize.size(3) = winsize.i3; end
  if isfield(winsize, 'i4'), winsize.size(4) = winsize.i4; end
  if isfield(winsize, 'i5'), winsize.size(5) = winsize.i5; end
end

if isfield(winsize, 'size') && isscalar(winsize.size)
  winsize.size(1:5) = winsize.size;
end

if ~isfield(winsize, 'winexp'), winsize.winexp = 1; end
if ~isfield(winsize, 'used'), winsize.used = 1; end
if ~isfield(winsize, 'winType'), winsize.winType = 'RaisedCos'; end


%% apply filter
if winsize.used
  r = zeros(max(size(toFill)), 5);
  for t = 1:numel(size(toFill))
    if mod(size(toFill, t), 2)
      r(1:size(toFill,t),t) = linspace((-0.5)/winsize.size(t), (0.5)/winsize.size(t), size(toFill,t)).';
      if size(toFill, t) == 1, r(1,t) = 0; end
    else
      r(1:size(toFill,t),t) = linspace((-0.5)/winsize.size(t), (0.5-1/size(toFill,t))/winsize.size(t), size(toFill,t)).';
    end
  end

  [myr1, myr2, myr3, myr4, myr5] = ...
    ndgrid(r(1:size(toFill,1),1), r(1:size(toFill,2),2), r(1:size(toFill,3),3), r(1:size(toFill,4),4), r(1:size(toFill,5),5));
  myr = (myr1.^2 + myr2.^2 + myr3.^2 + myr4.^2 + myr5.^2).^0.5;
  switch winsize.winType
    case 'RaisedCos'
      myr(myr>0.5) = 0.5;
      mywin = cos(myr*pi).^(2.*winsize.winexp);
    case 'Gaussian'
      mywin = exp(- ((myr*exp(0.5)*2).^2).*winsize.winexp);
    otherwise
      error('PD:zeroFill:UnknownFilter', 'Use one of the filter window functions: "RaisedCos" or "Gaussian".')
  end

  mywin(mywin<eps('double')/1e6) = 0;
  toFill = toFill.*mywin;
end


%% pad with zeros
filled = zeros(outputSize);
filled(...
  (floor(end/2)-floor(size(toFill,1)/2))+(1:size(toFill,1)),...
  (floor(end/2)-floor(size(toFill,2)/2))+(1:size(toFill,2)),...
  (floor(end/2)-floor(size(toFill,3)/2))+(1:size(toFill,3)),...
  (floor(end/2)-floor(size(toFill,4)/2))+(1:size(toFill,4)),...
  (floor(end/2)-floor(size(toFill,5)/2))+(1:size(toFill,5))) = toFill;


%%
if 0

%%
i1=11;
i2=12;
i3=17;
i1o=5;
i2o=5;
i3o=5;
i1z=11*5;
i2z=12*4;
i3z=17*2;
winsize.i1=1;
winsize.i2=1;
winsize.i3=1;

myImage=zeros(i1,i2,i3);
myImage(floor(end/2)+(1:i1o)-floor(i1o/2),floor(end/2)+(1:i2o)-floor(i2o/2),floor(end/2)+(1:i3o)-floor(i3o/2))=ones(i1o,i2o,i3o)+2i;
[Filled_image,kSpaceZ,kSpace]=zeroFill_image(myImage,[i1z,i2z,i3z],winsize);
%             sliceomatic(angle(myImage))
%             sliceomatic(abs(kSpace))
%             sliceomatic(abs(kSpaceZ))
            sliceomatic(abs(Filled_image))
%             sliceomatic(angle(Filled_image))
%             sliceomatic(angle(kSpace))
% sliceomatic(angle(kSpaceZ))
%%




dataOut(1).image_meanZ=zeroFill_image(dataOut(1).image_mean,size(dataOut(1).image_mean)*2,2);
dataOut(1).image_mean=(dataOut(1).image+dataOut(2).image+dataOut(2).image)/3;
dataOut(1).kSpace_mean=fftshift(ifftn(fftshift(dataOut(1).image_mean)));
dataOut(1).kSpaceZ_mean=zeroFill(dataOut(1).kSpace_mean,size(dataOut(1).kSpace_mean),0.5);
dataOut(1).kSpaceZ_mean=zeroFill(ones([9,7,8].*[1,2,1]),[10,8,9]*2+[1,0,1]);
dataOut(1).image_meanZ=(fftshift(fftn(fftshift(dataOut(1).kSpaceZ_mean))));
dataOut(1).kSpaceZLF_mean=zeroFill(dataOut(1).kSpace_mean,size(dataOut(1).kSpace_mean)*2,0.2);
dataOut(1).imageLF_meanZ=(fftshift(fftn(fftshift(dataOut(1).kSpaceZLF_mean))));
% dataOut(1).meanZ=(fftshift(fftn(fftshift(ifftn(fftshift(dataOut(1).mean))),size(dataOut(1).mean)*2)));
%             for t=1:size(dataOut(1).meanZ,3)
%                dataOut(1).meanZ_YXZ(:,:,t)=dataOut(1).meanZ(:,:,t).';
%             end
            sliceomatic(angle(dataOut(1).image_mean))
            sliceomatic(abs(dataOut(1).image_mean))

            sliceomatic(abs(dataOut(1).kSpace))
            sliceomatic(angle(dataOut(1).kSpace))
            sliceomatic(abs(dataOut(1).kSpaceZ_mean))
            sliceomatic(angle(dataOut(1).kSpaceZ_mean))
            sliceomatic(abs(dataOut(1).image_meanZ))
            figure

            imagesc(squeeze(mean(abs(dataOut(1).image_meanZ(:,:,270:300)),3)))
            set(gca,'YDir','normal')
            sliceomatic(abs(dataOut(1).imageLF_meanZ))
            sliceomatic(angle(dataOut(1).image_meanZ))
            tEcho=15e-3;
            Frequency=23e6;
%%
dataOut(1).ppmDiff=(unwrap3Dmiddle(angle(dataOut(1).image_meanZ))-unwrap3Dmiddle(angle(dataOut(1).imageLF_meanZ)))/tEcho/2/pi/(Frequency/1e6)

limit=2
dataOut(1).ppmDiff(dataOut(1).ppmDiff>limit)=limit;
dataOut(1).ppmDiff(dataOut(1).ppmDiff<-limit)=-limit;
            sliceomatic(dataOut(1).ppmDiff)


%%
%%
dataOut(1).image_mean=(dataOut(1).image+dataOut(2).image+dataOut(3).image+dataOut(4).image)/4;
dataOut(1).image_meanabs=abs(dataOut(1).image_mean);

dataOut(1).image_meanZ=zeroFill_image(dataOut(1).image_meanabs,size(dataOut(1).image_mean)*4,1.4);

tEcho=15e-3;
winsize.i1=2;
winsize.i2=2;
winsize.i3=2;
dataOut(1).image_meanZHF=unwrap3Dmiddle(angle(zeroFill_image(dataOut(1).image_mean,size(dataOut(1).image_mean)*1,winsize)))/tEcho/2/pi/(Frequency/1e6);
winsize.i1=0.1;
winsize.i2=0.1;
winsize.i3=0.1;
dataOut(1).image_meanZLF=real(zeroFill_image(dataOut(1).image_meanZHF,size(dataOut(1).image_mean)*1,winsize));
%             sliceomatic((dataOut(1).image_meanZLF))
%             sliceomatic((dataOut(1).image_meanZHF))
%             sliceomatic(dataOut(1).image_meanZHF-dataOut(1).image_meanZLF)
dataOut(1).ppmDiff=(dataOut(1).image_meanZHF-dataOut(1).image_meanZLF);
limit=0.5;
dataOut(1).ppmDiff(dataOut(1).ppmDiff>limit)=limit;
dataOut(1).ppmDiff(dataOut(1).ppmDiff<-limit)=-limit;
            sliceomatic(dataOut(1).ppmDiff)


end

end
