clc;
clear all;
%coordinates of ground station bangalore
deg2rad = pi/180;
rad2deg = 180/pi;
latd=13.03; 
lond=77.50;%in degrees
H=0.838;%km
re=6378;%km
ae=6378.1414;
latr=latd*deg2rad; 
lonr=lond*deg2rad;%in radians
Start=2458317.67012; %Start time in julian date 17th july 2018 (04:05:00)UTC
%t2=2458318.67012; %End time in julian date 18th july 2018 (04:05:00)UTC
End=24583.67012;
JD=End; %End Time
JDref=2415020; 
JC=36525; %Julian century from 1900
Startfrac=0.17012; %Fractional part of day (17th july 2018), time at which
%TLE elements or keplerian elements were obtained

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
duration=End-Start; %Enter for how much time simulation has to be run in days.  
days=duration+0.17012; %Enter no of days for simulation from epoch
min=1440;
sec=86400;
noofdays=duration;
%JD=linspace(t1,t2,1440); %updating minute by minute
T=(JD-JDref)./JC;
%UTday=0.041666.*(H+(m/60)+(S/3600));
%UTdeg=UTday.*360;
UT=linspace(0.17012,days,duration*min); %Fractional part of the day
%starting from 4:05:00 UTC
UTD=UT.*360;
UTD=mod(UTD,360);
x=(36000.7689.*T);
X=mod(x,360);
y=0.0004.*(T.^2);
Y=mod(y,360);
gmst=99.6910+X+Y+UTD;
gmst=mod(gmst,360);
lmst=gmst+lond;
lmst=mod(lmst,360);
u=3.986005*10^14;

%IRS P6
in=98.5525;%inclination
ohm0=276.0870;%RAAN
e0=0.00563;%eccentricity
omega0=303.7850;%Argument of perigee
v0=55.928;%True anamoly = Mean anamoly
n0=14.34156; %Mean motion
sat='IRS P6';

% % IRS P5
% in=97.873;
% ohm0=257.9920;
% e0=0.001044;
% omega0=51.2651;
% v0=308.8671;
% n0=14.787;
% sat='IRS P5';

% %CARTOSAT 2AT
% in=97.82;
% ohm0=218.5088;
% e0=0.002763;
% omega0=45.0401;
% v0=315.1040;
% n0=14.834760;
% sat='CARTOSAT 2AT';
% 
% %CARTOSAT 2A
% in=97.8747;
% ohm0=258.1189;
% e0=0.014365;
% omega0=113.7127;
% v0=246.5598;
% n0=14.7869;
% sat='CARTOSAT 2A'; 

% %CARTOSAT 2B
% in=97.4082;
% ohm0=259.6822;
% e0=0.014936;
% omega0=183.8495;
% v0=176.2627;
% n0=14.7869;
% sat='CARTOSAT 2B';
% 
% % CARTOSAT 2C
% in=97.4082;
% ohm0=259.6822;
% e0=0.0014936;
% omega0=129.0233;
% v0=288.5499;
% n0=15.1924;
% sat='CARTOSAT 2C';
% 
% %CARTOSAT 2D
% in=97.4428;
% ohm0=260.1310;
% e0=0.005160;
% omega0=281.3608;
% v0=78.7046;
% n0=15.1925;
% sat='CARTOSAT 2D';

% %CARTOSAT 2E
% in=97.3975;
% ohm0=256.9172;
% e0=0.0080;
% omega0=332.8002;
% v0=27.2809;
% n0=15.1925;
% sat='CARTOSAT 2E';

% %OCEANSAT 2
% in=98.3062;
% ohm0=293.7925;
% e0=0.003537;
% omega0=67.0117;
% v0=293.1451;
% n0=14.50848;
% sat='OCEANSAT 2';

% %RESOURCESAT 2A
% in=98.6703;
% ohm0=273.5977;
% e0=0.000619;
% omega0=109.6629;
% v0=250.4627;
% n0=14.21654;
% sat='RESOURCESAT 2A';

% %RISAT 1
% in=97.5702;
% ohm0=206.8083;
% e0=0.005542;
% omega0=115.9290;
% v0=293.4380;
% n0=15.0955;
% sat='RISAT 1';

% %SARAL
% in=98.5409;
% ohm0=25.3311;
% e0=0.001841;
% omega0=124.3523;
% v0=235.7835;
% n0=14.3201;
% sat='SARAL';

% %IRNSS R1A
% in=29.4959;
% ohm0=108.6179;
% e0=0.01787;
% omega0=178.4468;
% v0=181.5322;
% n0=1.002;
% sat='IRNSS R1A';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrad=2*pi*n0/86400;%Mean motion in rad/s
ndd=nrad*86400*rad2deg;%Mean motion in deg/day
v=v0+(ndd.*duration);%Total motion in terms of days

vrot=linspace(v0,v,duration*min);%True anamoly 
%for total duration updated minute by minute
vrot=mod(vrot,360);%scaled to 360 deg format

a=(u./(nrad.^2))^0.333;%Semimajor axis a
a=a/1000; %km
rapogee=a.*(1+e0);%apogee
rperigee=a.*(1-e0);%perigee
ha=rapogee-re;%Apogee height
hp=rperigee-re;%Perigee height
k1=66063.17;%km^2

%for l=linspace((Startfrac*min),(duration*min),(duration*min))
l=linspace((Startfrac*min),(duration*min),(duration*min)); %Calculating for 
%number of days in minutes
k=(ndd.*k1)./((a.^2).*(1-e0.^2)); %deg/day
dohmdt=cos(in.*deg2rad).*k.*(-1); %Rate of regression of nodes (Related to RAAN)
ohm=ohm0+(dohmdt.*l);
ohm=mod(ohm,360);

domegadt=k.*(2-2.5.*(sin(in.*deg2rad).^2));%deg/day
omega=omega0+(domegadt.*l);
omega=mod(omega,360);
omegar=omega.*deg2rad;
% scatter(l,ohm);
% grid on;
% hold on;
% xlim([0 5000]);
% ylim([0 370]);
% pause(0.01)
% end
gmstr=gmst.*deg2rad;

r=(a.*(1-(e0.^2)))/(1+e0.*cos(v0));
vrotr=vrot.*deg2rad;
rp=r.*(cos(vrotr));
rq=r.*(sin(vrotr));
r1=(cos(ohm).*cos(omega)); 
r2=(sin(ohm).*sin(omega).*cos(in)); 
r3=((-1).*cos(ohm).*sin(omega)); 
r4=(sin(ohm).*cos(omega).*cos(in)); 
r5=(sin(ohm).*cos(omega)); 
r6=cos(ohm).*sin(omega).*cos(in); 
r7=((-1).*sin(ohm).*sin(omega)); 
r8=cos(ohm).*cos(omega).*sin(in);
r9=(sin(omega).*sin(in));
r10=(cos(omega).*sin(in));
R12=r1-r2; 
R34=r3-r4; 
R56=r5+r6; 
R78=r7+r8; 
R9=r9; 
R10=r10;
% alpha=(atan(cos(in).*tan(omegar+vrotr)).*rad2deg)+ohm;
% alpha=mod(alpha,360);
% delta=(asin(sin(in).*sin(omegar+vrotr)).*rad2deg);
% delta=mod(delta,360);
% lons=alpha-gmst;
% lons=mod(lons,360);


R=[R12;R56;R9;R34;R78;R10];
Rm=[R12 R34;R56 R78;R9 R10];%R'
Rd=[rp;rq];%|rp|
           %|rq|
% rI=(R12.*rp)+(R34.*rq);
% rJ=(R56.*rp)+(R78.*rq);
% rK=(R9.*rp)+(R10.*rq);
% S=(rI.^2)+(rJ.^2)+(rK.^2);
% rv=sqrt(S);
for i=1:1:duration*min
R2=R(1:3,i)*rp;%r1 r3 r5 
R3=R(4:6,i)*rq;%r2 r4 r6
rI=R2(1,1:i)+R3(1,1:i);%rI
rJ=R2(2,1:i)+R3(2,1:i);%rJ
rK=R2(3,1:i)+R3(3,1:i);%rK
rIJK=[rI;rJ;rK];
% S=(((rI(1,1:i)).^2)+((rJ(1,1:i)).^2)+((rK(1,1:i)).^2));
% rv=sqrt(S);
end
xd=((rI.*cos(gmstr))+(rJ.*sin(gmstr)));
yd=(((-1).*rI.*sin(gmstr))+rJ.*cos(gmstr));
zd=rK;
lons=atan(yd./xd);
lats=atan(zd./sqrt((xd.^2)+(yd.^2)));
lons=(lons.*rad2deg); 
lats=(lats.*rad2deg);
S=(rI.^2)+(rJ.^2)+(rK.^2);
rv=sqrt(S);
Ee=0.08182;%Eccentricity of the earth
n1=1-((Ee.^2).*(sin(lonr).^2));
N=ae./sqrt(n1);
L=(N+H).*cos(lonr);
lmstr=lmst.*deg2rad;
RI=L.*cosd(lmstr);
RJ=L.*sind(lmstr);
RK=(N.*(1-(Ee.^2))+H).*sind(lonr);%Z
rooI=rI-RI;
rooJ=rJ-RJ;
rooK=rK-RK;
siE=atan(RK./L);
siED=siE.*rad2deg;
w=siE;
s11=sin(w).*cos(lmstr); 
s12=sin(w).*sin(lmstr); 
s13=(-1).*cos(w);
s21=(-1).*sin(lmstr); 
s22=cos(lmstr); 
s23=0;
s31=cos(w).*cos(lmstr); 
s32=cos(w).*sin(lmstr); 
s33=sin(w);
S13=(1:(duration*min)); 
S23=(1:(duration*min)); 
S33=(1:(duration*min));
for i=1:1:(duration*min)
    S13(i)=s13;
    S23(i)=s23;
    S33(i)=s33;
end
%rooSEZ=[s11 s12 S13;s21 s22 S23;s31 s32 S33];
roosez=[s11;s21;s31;s12;s22;s32;S13;S23;S33];
for i=1:1:duration*min
rooA=roosez(1:3,i)*rooI;
rooB=roosez(4:6,i)*rooJ;
rooC=roosez(7:9,i)*rooK;
% rooS=rooA(1,1:i)+rooB(1,1:i)+rooC(1,1:i);
% rooE=rooA(2,1:i)+rooB(2,1:i)+rooC(2,1:i);
% rooZ=rooA(3,1:i)+rooB(3,1:i)+rooC(3,1:i);

% roo=sqrt(de);
% des=(rooZ./roo);
rooS=(s11.*rooI)+(s12.*rooJ)+(S13.*rooK);
rooE=(s21.*rooI)+(s22.*rooJ)+(S23.*rooK);
rooZ=(s31.*rooI)+(s32.*rooJ)+(S33.*rooK);
de=(rooS(1,1:i).^2)+(rooE(1,1:i).^2)+(rooZ(1,1:i).^2);
% de=(rooA(1,1:i).^2)+(rooB(1,1:i).^2)+(rooC(1,1:i).^2);
rooSEZ=[rooS;rooE;rooZ];
% de=(rooS.^2)+(rooE.^2)+(rooZ.^2);
% roo=sqrt(de);
% des=(rooZ./roo);
end
% rooS=(s11.*rooI)+(s12.*rooJ)+(s13.*rooK);
% rooE=(s21.*rooI)+(s22.*rooJ)+(s23.*rooK);
% rooZ=(s31.*rooI)+(s32.*rooJ)+(s33.*rooK);
% rooSEZ=[rooS;rooE;rooZ];
% % de=(rooS.^2)+(rooE.^2)+(rooZ.^2);
% % roo=sqrt(de);
% % des=(rooZ./roo);
% de=(rooS.^2)+(rooE.^2)+(rooZ.^2);
% roo=sqrt(de);
% des=(rooZ./roo);
% Elevation=asind(des);
%Elevation=Elevation.*rad2deg;
%Elevation=(Elevation); 
%Elevation=mod(Elevation,180);
roos=(1:duration*min);
rooe=(1:duration*min);
for i=1:1:duration*min
    roos(i)=rooS(i);
    rooe(i)=rooE(i);
% if(rooE(i)<0)
%     rooE(i)=rooE(i).*(-1);
% 
% end
% end
% for i=1:1:duration*min
% if(rooS(i)<0)
%     rooS(i)=rooS(i).*(-1);
% end
 end
Azimuth=atand(abs(rooE)./abs(rooS));
%Azimuth=mod(Azimuth,360);
% %Azimuth=Azimuth.*rad2deg;
azi=Azimuth;
%Azimuth=mod(Azimuth,360);
roo=sqrt(de);
des=((rooZ)./(roo));
Elevation=asind(des);

for i=1:1:duration*min
if((roos(i)>0)&&(rooe(i)>0))
   Azimuth(i)=180-Azimuth(i);
%Azimuth(i)=mod(180-Azimuth,360);
end
if((roos(i)>0)&&(rooe(i)<0))
        Azimuth(i)=Azimuth(i)+180; 
        %Azimuth=mod(Azimuth,360);
end
if((roos(i)<0)&&(rooe(i)<0))
    Azimuth(i)=360-Azimuth(i);
    %Azimuth=mod(360-Azimuth,360);
end
end
figure('NumberTitle','off','Name',sat);
subplot(4,1,1);
plot(l,Azimuth);
grid on;
xlabel('Time in minutes');
ylabel('Azimuth in degrees');
subplot(4,1,2);
plot(l,Elevation);
grid on;
xlabel('Time in minutes');
ylabel('Elevation in degrees');
subplot(4,1,3);
plot(l,lats);
grid on;
xlabel('Time in minutes');
ylabel('Latitude');
subplot(4,1,4);
plot(l,lons);
grid on;
xlabel('Time in minutes');
ylabel('Longitude');

disp('Simulation Length ');
disp(noofdays);




