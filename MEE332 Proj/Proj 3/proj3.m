syms t;

%time
tau=20;
t1=0:1:tau;

%{

keep the number of points 0 to tau less than 20
the code will be slower

tau=.2;
t1=0:.01:tau;
tau=2;
t1=0:.1:tau;
tau=20;
t1=0:1:tau;
tau=.200;
t1=0:10:tau;
%}

% masses
mpl=.1; %kg per meter
mm=2; %hanging mass

%gravitational accel
ga=9.81;

%lengths
rae=16;
rac=29;
rof=13;
roe=16;
rdc=3;
rbd=35;
rdf=13;
rfe=roe-rof;
rec=rac-rae;
rbc=rbd-rdc;

xo=0;
yo=0;
xa=0;
ya=-17.5-12.5.*cos((pi.*t)./tau);
roa=-ya;
thetao=(acosd((rae.^2-roa.^2-roe.^2)./(-2.*roa.*roe))); %angle AOF
thetaa=(acosd((roe.^2-roa.^2-rae.^2)./(-2.*roa.*rae))); %angle EAO
thetae=180-thetao-thetaa; %angle AEO

%diagonal
rfc=sqrt(rfe.^2+rec.^2-2.*rfe.*rec.*cosd(180-thetae));

%point d
thetad=acosd((rfc.^2-rdf.^2-rdc.^2)./(-2.*rdf.*rdc)); %angle CDF

%point f and c
thetaf=asind((rdc./rfc).*sind(thetad))+acosd(((rfc.^2)+(rfe.^2)-(rec.^2))/(2.*rfc.*rfe)); %angle EFD
thetac=360-thetad-thetaf-(180-thetae); %angle DCE

%angular velocities
wf=diff(thetaf); %angle EFD
we=diff(180-thetae); %angle CEF
wd=diff(thetad); %angle CDF
wc=diff(-(180-thetac)); %angle BCE
wo=diff(thetao); %angle AOF
wa=diff(thetaa); %angle EAO

%angular accelerations
af=diff(wf);
ae=diff(we);
ad=diff(wd);
ac=diff(wc);
ao=diff(wo);
aa=diff(wa);

%trajectory
%f
xf=xo+rof.*sind(thetao);
yf=yo-rof.*cosd(thetao);

%e
xe=xa+rae.*sind(thetaa);
ye=yo+ya+rae.*cosd(thetaa);

%d
xd=xf+rdf.*cosd(thetaf-(90-thetao));
yd=yf+rdf.*sind(thetaf-(90-thetao));

%c
xc=xe+rec.*cosd(180-((180-thetae)+(90-thetao)));
yc=ye+rec.*sind(180-((180-thetae)+(90-thetao)));

%b
xb=rac.*cosd(90-thetaa)+rbc.*sind(180-thetaa-thetac);
yb=ya+rac.*sind(90-thetaa)-rbc.*cosd(180-thetaa-thetac);

%masses of rods
moe=mpl.*roe;
mac=mpl.*rac;
mdf=mpl.*rdf;
mbd=mpl.*rbd;

%center of gravity
coe=roe./2;
cac=rac./2;
cdf=rdf./2;
cbd=rbd./2;

%position of center of gravity
%x
xpcoe=coe.*sind(thetao);
xpcac=cac.*sind(thetaa);
xpcdf=xf+cdf.*cosd(thetaf-(90-thetao));
xpcbd=rac.*cosd(90-thetaa)+(cbd-rdc).*sind(180-thetac-thetaa);

%y
ypcoe=-coe.*cosd(thetao);
ypcac=yo+ya+cac.*cosd(thetaa);
ypcdf=yf+cdf.*sind(thetaf-(90-thetao));
ypcbd=ya+rac.*sind(90-thetaa)-(cbd-rdc).*cosd(180-thetaa-thetac);

%acceleration of center of gravity
%x
acoex=diff(xpcoe,2);
acacx=diff(xpcac,2);
acdfx=diff(xpcdf,2);
acbdx=diff(xpcbd,2);

%y
acoey=diff(ypcoe,2);
acacy=diff(ypcac,2);
acdfy=diff(ypcdf,2);
acbdy=diff(ypcbd,2);

%moment of inertia
iac=((rac.^2).*mac)./(12);
ioe=(((roe.^2).*moe)./(12))+(moe.*coe.^2);
idf=((rdf.^2).*mdf)./(12);
ibd=((rbd.^2).*mbd)./(12);

%pre declare
    Fox=zeros(1,length(t1));
    Foy=zeros(1,length(t1));
    Fax=zeros(1,length(t1));
    Fay=zeros(1,length(t1));
    Fcx=zeros(1,length(t1));
    Fcy=zeros(1,length(t1));
    Fdx=zeros(1,length(t1));
    Fdy=zeros(1,length(t1));
    Fex=zeros(1,length(t1));
    Fey=zeros(1,length(t1));
    Ffx=zeros(1,length(t1));
    Ffy=zeros(1,length(t1));

%matrix eqn
for k=1:length(t1)
    t2=t1(k)
    
    thetao1=(subs(thetao,{t}, {t2}));
    thetaa1=(subs(thetaa,{t}, {t2}));
    thetaf1=(subs(thetaf,{t}, {t2}));
    thetac1=(subs(thetac,{t}, {t2}));

    acoex1=(subs(acoex,{t}, {t2}));
    acacx1=(subs(acacx,{t}, {t2}));
    acdfx1=(subs(acdfx,{t}, {t2}));
    acbdx1=(subs(acbdx,{t}, {t2}));

    acoey1=(subs(acoey,{t}, {t2}));
    acacy1=(subs(acacy,{t}, {t2}));
    acdfy1=(subs(acdfy,{t}, {t2}));
    acbdy1=(subs(acbdy,{t}, {t2}));

    ao1=(subs(ao,{t}, {t2}));
    aa1=(subs(aa,{t}, {t2}));
    ac1=(subs(ac,{t}, {t2}));
    af1=(subs(af,{t}, {t2}));

    A=[-1 0 0 0 0 0 0 0 1 0 1 0;0 -1 0 0 0 0 0 0 0 1 0 1;0 0 0 0 0 0 0 0 coe.*cosd(thetao1) coe.*sind(thetao1) (roe-coe).*cosd(thetao1) (roe-coe).*sind(thetao1);0 0 -1 0 1 0 0 0 -1 0 0 0;0 0 0 -1 0 1 0 0 0 -1 0 0;0 0 -cac.*cosd(thetaa1) cac.*sind(thetaa1) -cac.*cosd(thetaa1) cac.*sind(thetaa1) 0 0 -(rec-cac).*cosd(thetaa1) (rec-cac).*sind(thetaa1) 0 0;0 0 0 0 0 0 1 0 0 0 -1 0;0 0 0 0 0 0 0 1 0 0 0 -1;0 0 0 0 0 0 -cdf.*sind(thetaf1-(90-thetao1)) cdf.*cosd(thetaf1-(90-thetao1)) 0 0 -cdf.*sind(thetaf1-(90-thetao1)) cdf.*cosd(thetaf1-(90-thetao1));0 0 0 0 -1 0 -1 0 0 0 0 0;0 0 0 0 0 -1 0 -1 0 0 0 0;0 0 0 0 (cbd-rdc).*sind(180-thetac1-(90-thetaa1)) (cbd-rdc).*cosd(180-thetac1-(90-thetaa1)) cbd.*sind(180-thetac1-(90-thetaa1)) cbd.*cosd(180-thetac1-(90-thetaa1)) 0 0 0 0];
    b=[moe.*acoex1;moe.*acoey1;ioe.*(-ao1);mac.*acacx1;mac.*acacy1;iac.*(-aa1);mdf.*acdfx1;mdf.*acdfy1;idf.*af1;mbd.*acbdx1;(mbd.*acbdy1+mm.*ga);(ibd.*ac1+(mm.*ga.*cosd(180-thetac1-(90-thetaa1))))];
    x=inv(A)*b;
    
    Fox(k)=x(1);
    Foy(k)=x(2);
    Fax(k)=x(3);
    Fay(k)=x(4);
    Fcx(k)=x(5);
    Fcy(k)=x(6);
    Fdx(k)=x(7);
    Fdy(k)=x(8);
    Fex(k)=x(9);
    Fey(k)=x(10);
    Ffx(k)=x(11);
    Ffy(k)=x(12);
end

figure(1);
hold on;
grid on;
title("Forces of joints O A C D E F x direction");
plot(t1,Fox)
plot(t1,Fax)
plot(t1,Fcx)
plot(t1,Fdx)
plot(t1,Fex)
plot(t1,Ffx)
xlabel('time (s)') ;
ylabel('Force (N)') ;
legend('Fox','Fax','Fcx','Fdx','Fex','Ffx');

figure(2);
hold on;
grid on;
title("Forces of joints O A C D E F y direction");
plot(t1,Foy)
plot(t1,Fay)
plot(t1,Fcy)
plot(t1,Fdy)
plot(t1,Fey)
plot(t1,Ffy)
xlabel('time (s)') ;
ylabel('Force (N)') ;
legend('Foy','Fay','Fcy','Fdy','Fey','Ffy');
hold off;
