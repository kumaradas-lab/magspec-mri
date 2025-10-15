function [fh,interp_image,interp_image_corr]=plot_interp_image(image,AQSlice,winexp,winsize,PlotPhase,fh)
%% function [fh,interp_image,interp_image_corr]=plot_interp_image(image,AQSlice,winexp,winsize,PlotPhase,fh)

    if nargin<=2
        winexp=1;
    end
    if nargin<=3
        winsize=1.4;
    end
    if nargin<=4
        PlotPhase=0;
    end
    if nargin<=5
        fh=3330;
    end
    temp=winsize;
    clear winsize
    winsize.i1=temp;
    winsize.i2=temp;
    winsize.i3=temp;
    winsize.used=1;
    winsize.winexp=winexp;
    clear temp
    
    if ~exist('AQSlice','var'); AQSlice=[]; end
    if ~isfield(AQSlice,'sizePhase'); AQSlice.sizePhase=1; end
    if ~isfield(AQSlice,'sizeRead'); AQSlice.sizeRead=1; end

    
    fh=figure(fh);
    clf;
%     ringingread=linspace(-1/winsize,1/winsize-1/winsize/size(image,1),size(image,1)).';
%     ringingphase=linspace(-1/winsize,1/winsize-1/winsize/size(image,2),size(image,2));
%     myr=(ringingread.^2*ones(1,size(image,2))+ones(size(image,1),1)*ringingphase.^2).^0.5;
%     myr(myr>1)=1;
%     mywin=cos(myr*pi/2).^2.^winexp;
% %     imagesc(abs(mywin))
% %     pause
%     interp_image=((fft2((fftshift(ifft2(image)).*mywin),2^10,2^10)));
%     

    [interp_image] = zeroFill_image(squeeze(image),[2^10,2^10],winsize);

    if PlotPhase;subplot(2,1,1); else subplot(1,1,1);end
    imagesc(abs(interp_image));
%     colorbar
    colormap(gray);
%      xlabel('Phase');
%     ylabel('READ')    
    if nargin>=2
        set(gca,'DataAspectRatio',[1 AQSlice.sizePhase/AQSlice.sizeRead 1]);
    else
        axis equal
        axis tight
    end

    if PlotPhase 
        subplot(2,1,2)
        imagesc(unwrap2Dmiddle(angle(interp_image)));
    %     colorbar
        colormap(gray);
    %      xlabel('Phase');
    %     ylabel('READ')    
        if nargin>=2
            set(gca,'DataAspectRatio',[1 AQSlice.sizePhase/AQSlice.sizeRead 1]);
        else
            axis equal
            axis tight
        end
    end
    drawnow


    if nargout > 2
		 factor=[1 1];
		if AQSlice.nRead > AQSlice.nPhase
		   factor=[1 AQSlice.nPhase/AQSlice.nRead];
		elseif AQSlice.nRead < AQSlice.nPhase
		   factor=[AQSlice.nRead/AQSlice.nPhase 1];
		end
% 		interp_image_corr=((fft2((fftshift(ifft2(image)).*mywin),round(2^10 * factor(1)),round(2^10  * factor(2)))));
        [interp_image_corr] = zeroFill_image(image,[round(2^10 * factor(1)),round(2^10  * factor(2))],winsize);
    end
end

