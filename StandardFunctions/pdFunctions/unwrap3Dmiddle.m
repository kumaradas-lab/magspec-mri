function unwraped = unwrap3Dmiddle(tounwrap)
% sliceomatic(tounwrap)
[s(1),s(2),s(3),s(4),s(5),s(6)]=size(tounwrap);
tounwrap=squeeze(tounwrap);
unwraped12=zeros(size(tounwrap));
unwraped=zeros(size(tounwrap));
for n=1:size(tounwrap,3)
    unwraped12(:,:,n) = unwrap2Dmiddle(tounwrap(:,:,n));
end
unwraped3m=unwrap(squeeze(unwraped12(round(end/2),round(end/2),:)))-(squeeze(unwraped12(round(end/2),round(end/2),:)));
for n=1:size(tounwrap,3)
    unwraped(:,:,n) = unwraped12(:,:,n)+unwraped3m(n);
end
unwraped=unwraped+(tounwrap(round(end/2),round(end/2),round(end/2))-unwraped(round(end/2),round(end/2),round(end/2))); % Set phase of middle voxel
unwraped=reshape(unwraped,s(1),s(2),s(3),s(4),s(5),s(6));



% diff1=diff(tounwrap,1,1);
% diff1(abs(diff1)<pi)=0;
% diff1s=diff1;
% diff1s((diff1)<-pi)=2*pi;
% diff1s((diff1)>pi)=-2*pi;
% sliceomatic(diff1s)
% 
% cumstep1=cumsumNAN(diff1s,1);
% cumstep1(isnan(diff1s))=0;
% sliceomatic(cumstep1)
% 
% u1a=(tounwrap);
% u1a(2:end,:,:)=u1a(2:end,:,:)+cumstep1;
% sliceomatic(u1a,image3D.xV*1000000,image3D.yV*1000000,image3D.zV*1000000)
% 
% sliceomatic(u1a)
% diff2=diff(u1a,1,2);
% diff2(abs(diff2)<pi)=0;
% diff2s=diff2;
% diff2s((diff2)<-pi)=2*pi;
% diff2s((diff2)>pi)=-2*pi;
% sliceomatic(diff2s)
% 
% cumstep2=cumsumNAN(diff2s,2);
% cumstep2(isnan(diff2s))=0;
% sliceomatic(cumstep2)
% 
% u2a=(u1a);
% u2a(:,2:end,:)=u2a(:,2:end,:)+cumstep1;
% sliceomatic(u2a,image3D.xV*1000000,image3D.yV*1000000,image3D.zV*1000000)
% 
% 
% 
% 
% sliceomatic(tounwrap,image3D.xV*1000000,image3D.yV*1000000,image3D.zV*1000000)
% 
% u1=unwrap((tounwrap),[],1);
% u1m=unwrap(u1(end/2,:,:),[],2)-u1(end/2,:,:);
% u1a=u1+reshape(ones(size(u1,1),1)*u1m(:).',size(u1));
% 
% sliceomatic(u1a,image3D.xV*1000000,image3D.yV*1000000,image3D.zV*1000000)
% 
% u2=unwrap((u1a),[],2);
% u2m=unwrap(u2(end/2,:,:),[],3)-u2(end/2,:,:);
% u2a=u2+reshape(ones(size(u2,1),1)*u2m(:).',size(u2));
% 
% sliceomatic(u2a,image3D.xV*1000000,image3D.yV*1000000,image3D.zV*1000000)
% 
% u3=unwrap((u2a),[],3);
% u3m=unwrap(u3(end/2,:,:),[],1)-u3(end/2,:,:);
% u3a=u3+reshape(ones(size(u3,1),1)*u3m(:).',size(u3));
% 
% sliceomatic(u3a,image3D.xV*1000000,image3D.yV*1000000,image3D.zV*1000000)
% 
% unwraped=u3a-(u3a(end/2,end/2,end/2)-(tounwrap(end/2,end/2,end/2)));
