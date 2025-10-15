function mywin=RaisedCosine2(n1,n2)
    r1=linspace(-1+1/n1,1-1/n1,n1).';
    r2=linspace(-1+1/n1,1-1/n1,n2);
    myr=(r1.^2*ones(1,n2)+ones(n1,1)*r2.^2).^0.5;
    myr(myr>1)=1;
    mywin=cosd(myr*90).^2;

