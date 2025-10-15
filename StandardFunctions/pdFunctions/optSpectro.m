function output=optSpectro(mydatafft,fOffsetFft,pOffsetFft,x)

    spectrum=mydatafft.*exp(-1i.*((fOffsetFft+x(1)).*(1:numel(mydatafft)).'+(pOffsetFft+x(2))));
    output=mean(abs(spectrum))-mean(real(spectrum));
%     disp([num2str(output ) '      ' num2str(x )]) 
    
