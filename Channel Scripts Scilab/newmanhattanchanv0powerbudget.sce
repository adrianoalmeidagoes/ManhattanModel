buildingwidth=200
Streetwidth=50
buildingheight=30;
heightcar=2;
heightBS=20;
dman=Streetwidth+buildingwidth;
speedlight=3e8;
er=1.2;
e0=8.854e-12;
e=er*e0;
nground=sqrt(e/e0);
lambda=0.255;
TxpowdB=30;
noisedensitydBHz=-140;
BW=1e5;
noisefiguredB=8;
noisepowerdB=BW*noisedensitydBHz;
noisepowerlin=10^(noisepowerdB+noisefiguredB);
TxgaindB=15;
RxgaindB=8;
totalpowlinear=10^((TxpowdB+TxgaindB+RxgaindB)/10);
xBs=520; yBs=275;
xcar=525; ycar=400;


Tsymbol=1e-3;  //Subcarrier spacing.

nreflec=0:1:10;
Mpath=zeros(1,30);
KdBprofileLOS=[30 10 0 -3 -6 zeros(1,25)];
KdBprofileNLOS=[-3 -10 -15 -20 -30 zeros(1,25)];

if abs(xcar-xBs)<Streetwidth || abs(ycar-yBs)<Streetwidth then
   dground=sqrt((xcar-xBs)^2+(ycar-yBs)^2+(heightcar+heightBS)^2);  
   dlos=sqrt((xcar-xBs)^2+(ycar-yBs)^2+(heightcar-heightBS)^2);
   if abs(xcar-xBs)<Streetwidth then
        distancexlos=abs(ycar-yBs);
        a=Streetwidth/2;
        b=abs(xcar-floor(xcar/dman)*dman);     
   else if abs(ycar-yBs)<Streetwidth then
        distancexlos=abs(xcar-xBs);
        a=Streetwidth/2;
        b=abs(ycar-floor(ycar/dman)*dman);
        end
    end
 
    y=a+nreflec*Streetwidth+b;
    distances=sqrt(y.^2+distancexlos.^2+(heightBS-heightcar)^2);
    proywall=sqrt(distancexlos.^2+(heightBS-heightcar)^2);
    cosbeta=proywall./distances;
    sinbeta=y./distances;
    refcoef=(sinbeta-(sqrt(real(nground)^2-cosbeta.^2)+%i*imag(nground)))./(sinbeta+(sqrt(real(nground)^2-cosbeta.^2)+%i*imag(nground)));
    y2=y-2*a-2*b+2*Streetwidth;
    distances2=sqrt(y2.^2+distancexlos.^2+(heightBS-heightcar)^2);
    cosbeta=proywall./distances2;
    sinbeta=y2./distances2;
    refcoef2=(sinbeta-(sqrt(real(nground)^2-cosbeta.^2)+%i*imag(nground)))./(sinbeta+(sqrt(real(nground)^2-cosbeta.^2)+%i*imag(nground)));
    RC1=refcoef.^(nreflec);
    RC2=refcoef2.^(nreflec);
    
    delaylos=distancexlos/speedlight;
    taplos=floor(delaylos/Tsymbol)+1;
    delayground=dground/speedlight;
    tapground=floor(delayground/Tsymbol)+1;
    delayreflecs1=distances/speedlight;
    tapref1=floor(delayreflecs1/Tsymbol)+1;
    delayreflecs2=distances2/speedlight;
    tapref2=floor(delayreflecs2/Tsymbol)+1;
    Mpath(taplos)=Mpath(taplos)+exp(-2*%pi*%i*distancexlos/lambda)/(4*%pi*distancexlos/lambda);
    Mpath(tapground)=Mpath(tapground)+exp(-2*%pi*%i*dground/lambda)/(4*%pi*dground/lambda);
    Mpath(tapref1)=Mpath(tapref1)+RC1.*exp(-2*%pi*%i*distances/lambda)./(4*%pi*distances/lambda);
    Mpath(tapref2)=Mpath(tapref2)+RC2.*exp(-2*%pi*%i*distances2/lambda)./(4*%pi*distances2/lambda);
    KdB=KdBprofileLOS;
else
   xcar0=xcar-Streetwidth/2; ycar0=ycar-Streetwidth/2;
   uy=abs(ycar-round(ycar0/dman)*dman-Streetwidth/2); ux=abs(xcar-round(xcar0/dman)*dman-Streetwidth/2);
   distancex=abs(round(xcar/dman)*dman-xBs)+Streetwidth/2; distancey=abs(round(ycar/dman)*dman-yBs)+Streetwidth/2;
   if uy<=Streetwidth && ux>=Streetwidth
       distance1=distancex; distance2=abs(xcar-xBs); b=abs(ycar-round(ycar/dman)*dman);
   elseif uy>=Streetwidth && ux<=Streetwidth
           distance1=distancey; distance2=abs(ycar-yBs); b=abs(xcar-round(xcar/dman)*dman);
   end
 
   a=Streetwidth/2;
   y1=a+nreflec*Streetwidth;  w1=b+nreflec*Streetwidth;
   miny2=((distance1-Streetwidth)./(y1+Streetwidth)).*(y1+distance2);
   minw2=((distance2-Streetwidth)./(w1+Streetwidth)).*(w1+distance1);
   
        expectednreflecy2=min(20, ceil((miny2-b)./Streetwidth)); 
        expectednreflecw2=min(20, ceil((minw2-a)./Streetwidth));   
        y2=b+expectednreflecy2*Streetwidth;
        w2=a+expectednreflecw2*Streetwidth
        distancea=sqrt((y1+distance2).^2+(y2+distance2)^2+(heightBS-heightcar)^2); 
        proywall=sqrt((y1+distance2).^2+(heightBS-heightcar)^2);
        cosbeta=proywall./distancea;
        sinbeta=y./distancea;
        refcoef=(sinbeta-(sqrt(real(nground)^2-cosbeta.^2)+%i*imag(nground)))./(sinbeta+(sqrt(real(nground)^2-cosbeta.^2)+%i*imag(nground)));
        distanceb=sqrt((w1+distance2).^2+(w2+distance2)^2); 
        proywall=sqrt((y2+distance2).^2+(heightBS-heightcar)^2);
        cosbeta=proywall./distanceb;
        sinbeta=y./distanceb;
        refcoef2=(sinbeta-(sqrt(real(nground)^2-cosbeta.^2)+%i*imag(nground)))./(sinbeta+(sqrt(real(nground)^2-cosbeta.^2)+%i*imag(nground)));
        RC1=refcoef.^(nreflec);
        RC2=refcoef2.^(nreflec);
        
        delayreflecs1=distancea/speedlight;
    tapref1=floor(delayreflecs1/Tsymbol)+1;
    delayreflecs2=distanceb/speedlight;
    tapref2=floor(delayreflecs2/Tsymbol)+1;
    Mpath(tapref1)=Mpath(tapref1)+RC1.*exp(-2*%pi*%i*distancea)./(4*%pi*distancea);
    Mpath(tapref2)=Mpath(tapref2)+RC2.*exp(-2*%pi*%i*distanceb)./(4*%pi*distanceb);
    KdB=KdBprofileNLOS;
   end
   K=10.^(KdB/10);
   scat=grand(1,length(Mpath),'nor',0,sqrt(1/2))+%i*grand(1,length(Mpath),'nor',0,sqrt(1/2));
   
   CHAN=Mpath+scat./sqrt(K);
   
   SNR=sum(abs(CHAN(1:5)).^2)*txpowlinear/noisepowerlinear;
   
   



