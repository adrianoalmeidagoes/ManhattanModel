clc();
clear

Streetwidth=50;
heightBS=20;
heightcar=2;



speedlight=3e8;
erground=1.1;
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



dloshorx=20;
dwally1=10;
dwally2=22;

Tsymbol=1e-3;  //Subcarrier spacing.

nreflec=1:1:10;
Mpath=zeros(1,30);
KdBprofileLOS=[30 10 0 -3 -6 zeros(1,25)];

dground=sqrt((dloshorx)^2+(dwally1-dwally2)^2+(heightcar+heightBS)^2);  
dlos=sqrt((dloshorx)^2+(dwally1-dwally2)^2+(heightcar-heightBS)^2);
dproyground=sqrt((dloshorx)^2+(dwally1-dwally2)^2);
cosbeta=dproyground./dground;
sinbeta=(heightcar+heightBS)./dground;
refcoefgr=(sinbeta-(sqrt(real(nground)^2-cosbeta.^2)+%i*imag(nground)))./(sinbeta+(sqrt(real(nground)^2-cosbeta.^2)+%i*imag(nground)));

y=dwally1+(nreflec-1)*Streetwidth+dwally2;
distances=sqrt(y.^2+dloshorx.^2+(heightBS-heightcar)^2);
proywall=sqrt(dloshorx.^2+(heightBS-heightcar)^2);
cosbeta=proywall./distances;
sinbeta=y./distances;
refcoef=(sinbeta-(sqrt(real(nwall)^2-cosbeta.^2)+%i*imag(nwall)))./(sinbeta+(sqrt(real(nwall)^2-cosbeta.^2)+%i*imag(nwall)));
    
y2=y-2*dwally1-2*dwally2+2*Streetwidth;
distances2=sqrt(y2.^2+dloshorx.^2+(heightBS-heightcar)^2);
cosbeta=proywall./distances2;
sinbeta=y2./distances2;
refcoef2=(sinbeta-(sqrt(real(nwall)^2-cosbeta.^2)+%i*imag(nwall)))./(sinbeta+(sqrt(real(nwall)^2-cosbeta.^2)+%i*imag(nwall)));
RC1=refcoef.^(nreflec);
RC2=refcoef2.^(nreflec);
    
delaylos=dlos/speedlight;
taplos=floor(delaylos/Tsymbol)+1;
delayground=dground/speedlight;
tapground=floor(delayground/Tsymbol)+1;
delayreflecs1=distances/speedlight;
tapref1=floor(delayreflecs1/Tsymbol)+1;
delayreflecs2=distances2/speedlight;
tapref2=floor(delayreflecs2/Tsymbol)+1;


Mpath(taplos)=exp(-2*%pi*%i*dlos/lambda)/(4*%pi*dlos/lambda);
Mpath(tapground)=Mpath(tapground)+refcoefgr*exp(-2*%pi*%i*dground/lambda)/(4*%pi*dground/lambda);
Mpath(tapref1)=Mpath(tapref1)+RC1.*exp(-2*%pi*%i*distances/lambda)./(4*%pi*distances/lambda);
Mpath(tapref2)=Mpath(tapref2)+RC2.*exp(-2*%pi*%i*distances2/lambda)./(4*%pi*distances2/lambda);
KdB=KdBprofileLOS;
K=10.^(KdB/10);
scat=grand(1,length(Mpath),'nor',0,sqrt(1/2))+%i*grand(1,length(Mpath),'nor',0,sqrt(1/2));
CHAN=Mpath+abs(Mpath).*scat./sqrt(K);
SNR=sum(abs(CHAN(1:5)).^2)*totalpowlinear/noisepowerlin;




   



