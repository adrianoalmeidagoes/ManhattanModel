    clear
    clc();
Streetwidth=50
heightcar=2;
heightBS=20;

speedlight=3e8;
erground=1.2;
e0=8.854e-12;
eground=erground*e0;
nground=sqrt(eground/e0);


erwall=1.5;
ewall=erwall*e0;
nwall=sqrt(ewall/e0);

lambda=0.255;

TxpowdB=30;
noisedensitydBHz=-140;
BW=1e9;
noisefiguredB=10;
noisepowerlin=BW*10^((noisedensitydBHz+noisefiguredB)/10);
TxgaindB=7;
RxgaindB=7;
totalpowlinear=10^((TxpowdB+TxgaindB+RxgaindB)/10);


Tsymbol=1e-3;

nreflec=2:1:10;
Mpath=zeros(1,30);

KdBprofileNLOS=[-3 -10 -15 -20 -30 zeros(1,25)];

distance1=250;
dwally1=20;
dwally2=30;
distance2=400;
 
y1=dwally1+nreflec*Streetwidth;  w1=dwally2+nreflec*Streetwidth;
miny2=((distance1-Streetwidth)./(y1+Streetwidth)).*(y1+distance2);
minw2=((distance2-Streetwidth)./(w1+Streetwidth)).*(w1+distance1);
   
expectednreflecy2=min(20, ceil((miny2-dwally2)./Streetwidth)); 
expectednreflecw2=min(20, ceil((minw2-dwally1)./Streetwidth));   

y2=dwally2+expectednreflecy2*Streetwidth;
w2=dwally1+expectednreflecw2*Streetwidth;
distancea=sqrt((y1+distance2).^2+(y2+distance2)^2+(heightBS-heightcar)^2); 
proywall=sqrt((y1+distance2).^2+(heightBS-heightcar)^2);
cosbeta=proywall./distancea;
sinbeta=y1./distancea;
refcoef=(sinbeta-(sqrt(real(nground)^2-cosbeta.^2)+%i*imag(nground)))./(sinbeta+(sqrt(real(nground)^2-cosbeta.^2)+%i*imag(nground)));
distanceb=sqrt((w1+distance2).^2+(w2+distance2)^2); 
proywall=sqrt((y2+distance2).^2+(heightBS-heightcar)^2);
cosbeta=proywall./distanceb;
sinbeta=y2./distanceb;
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
   
K=10.^(KdB/10);
scat=grand(1,length(Mpath),'nor',0,sqrt(1/2))+%i*grand(1,length(Mpath),'nor',0,sqrt(1/2));
   
CHAN=Mpath+abs(Mpath).*scat./sqrt(K);
SNR=sum(abs(CHAN(1:5)).^2)*totalpowlinear/noisepowerlin;




