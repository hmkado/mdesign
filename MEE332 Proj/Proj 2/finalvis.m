syms t;

%lengths
rae=16;
rac=29;
rof=13;
roe=16;
rdc=3;
rdb=35;
rdf=13;
rfe=roe-rof;
rec=rac-rae;
rbc=rdb-rdc;

%time
tau=2;
t1=0:.1:tau;

xo=0;
yo=0;
xa=0;
ya=-17.5-12.5.*cos((pi.*t)./tau);
roa=-ya;
thetao=(acosd((rae.^2-roa.^2-roe.^2)./(-2.*roa.*roe)));
thetaa=(acosd((roe.^2-roa.^2-rae.^2)./(-2.*roa.*rae)));
thetae=180-thetao-thetaa;

%diagonal
rfc=sqrt(rfe.^2+rec.^2-2.*rfe.*rec.*cosd(180-thetae));

%point d
thetad=acosd((rfc.^2-rdf.^2-rdc.^2)./(-2.*rdf.*rdc));

%point f and c
thetaf=asind((rdc./rfc).*sind(thetad))+acosd(((rfc.^2)+(rfe.^2)-(rec.^2))/(2.*rfc.*rfe));
thetac=360-thetad-thetaf-(180-thetae);

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

%plot of angular velocities
figure(1);
title("Angular velocity at joints");
grid on;
hold on;
plot(t1,(subs(wo, {t},{t1})),'-');
plot(t1,(subs(wa, {t},{t1})),'-');
plot(t1,(subs(wc, {t},{t1})),'-');
plot(t1,(subs(wd, {t},{t1})),'-');
plot(t1,(subs(we, {t},{t1})),'-');
plot(t1,(subs(wf, {t},{t1})),'-');
xlabel('time (s)') ;
ylabel('Angular Velocity (rad/s)') ;
legend('Angle AOF','Angle EAO','Angle BCE','Angle CDF','Angle CEF','Angle EFD');
hold off;

%plot of angular accelerations
figure(2);
title("Angular Acceleration at joints");
grid on;
hold on;
plot(t1,(subs(ao, {t},{t1})),'-');
plot(t1,(subs(aa, {t},{t1})),'-');
plot(t1,(subs(ac, {t},{t1})),'-');
plot(t1,(subs(ad, {t},{t1})),'-');
plot(t1,(subs(ae, {t},{t1})),'-');
plot(t1,(subs(af, {t},{t1})),'-');
xlabel('time (s)') ;
ylabel('Angular Accelration (rad/s^2)') ;
legend('Angle AOF','Angle EAO','Angle BCE','Angle CDF','Angle CEF','Angle EFD');
hold off;

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

%linear velocities in x direction
vxa=0; %no x direction
vxb=sqrt(diff(xb).^2+diff(yb).^2);
vxc=sqrt(diff(xc).^2+diff(yc).^2);
vxd=sqrt(diff(xd).^2+diff(yd).^2);
vxe=sqrt(diff(xe).^2+diff(ye).^2);
vxf=sqrt(diff(xf).^2+diff(yf).^2);

%y direction
vya=diff(ya);
vyb=diff(yb);
vyc=diff(yc);
vyd=diff(yd);
vye=diff(ye);
vyf=diff(yf);

%linear acceleration in x direction
axa=0; %no x direction
axb=diff(xb,2);
axc=diff(xc,2);
axd=diff(xd,2);
axe=diff(xe,2);
axf=diff(xf,2);

%y direction
aya=diff(ya,2); %no x direction
ayb=diff(yb,2);
ayc=diff(yc,2);
ayd=diff(yd,2);
aye=diff(ye,2);
ayf=diff(yf,2);

%linear velocity plot x direction
figure(3);
grid on;
hold on;
plot(t1,(subs(vxa, {t},{t1})),'-');
plot(t1,(subs(vxb, {t},{t1})),'-');
plot(t1,(subs(vxc, {t},{t1})),'-');
plot(t1,(subs(vxd, {t},{t1})),'-');
plot(t1,(subs(vxe, {t},{t1})),'-');
plot(t1,(subs(vxf, {t},{t1})),'-');
title("Linear Velocity of Joints x direction");
xlabel('time (s)') ;
ylabel('Velocity (cm/s)') ;
legend('A','B','C','D','E','F');
hold off;

%linear velocity plot x direction
figure(4);
grid on;
hold on;
plot(t1,(subs(vya, {t},{t1})),'-');
plot(t1,(subs(vyb, {t},{t1})),'-');
plot(t1,(subs(vyc, {t},{t1})),'-');
plot(t1,(subs(vyd, {t},{t1})),'-');
plot(t1,(subs(vye, {t},{t1})),'-');
plot(t1,(subs(vyf, {t},{t1})),'-');
title("Linear Velocity of Joints y direction");
xlabel('time (s)') ;
ylabel('Velocity (cm/s)') ;
legend('A','B','C','D','E','F');
hold off;

%linear acceleration plot x direction
figure(5);
grid on;
hold on;
plot(t1,(subs(axa, {t},{t1})),'-');
plot(t1,(subs(axb, {t},{t1})),'-');
plot(t1,(subs(axc, {t},{t1})),'-');
plot(t1,(subs(axd, {t},{t1})),'-');
plot(t1,(subs(axe, {t},{t1})),'-');
plot(t1,(subs(axf, {t},{t1})),'-');
title("Linear Acceleration of Joints x direction");
xlabel('time (s)') ;
ylabel('Acceleration (cm/s^2)') ;
legend('A','B','C','D','E','F');
hold off;

%linear acceleration plot y direction
figure(6);
grid on;
hold on;
plot(t1,(subs(aya, {t},{t1})),'-');
plot(t1,(subs(ayb, {t},{t1})),'-');
plot(t1,(subs(ayc, {t},{t1})),'-');
plot(t1,(subs(ayd, {t},{t1})),'-');
plot(t1,(subs(aye, {t},{t1})),'-');
plot(t1,(subs(ayf, {t},{t1})),'-');
title("Linear Acceleration of Joints y direction");
xlabel('time (s)') ;
ylabel('Acceleration (cm/s^2)') ;
legend('A','B','C','D','E','F');
hold off;

%plot trajectory of points A,B,C,D,E,F
figure(7);
grid on;
hold on;
axis equal;
title("Trajectories in XY plane");
p1=plot((subs(xa, {t},{t1})),(subs(ya, {t},{t1})),'-');
p2=plot((subs(xb, {t},{t1})),(subs(yb, {t},{t1})),'-');
p3=plot((subs(xc, {t},{t1})),(subs(yc, {t},{t1})),'-');
p4=plot((subs(xd, {t},{t1})),(subs(yd, {t},{t1})),'-');
p5=plot((subs(xe, {t},{t1})),(subs(ye, {t},{t1})),'-');
p6=plot((subs(xf, {t},{t1})),(subs(yf, {t},{t1})),'-');
xlabel('x position (cm)') ;
ylabel('y position (cm)') ;
legend('A','B','C','D','E','F');
axis([-5 65 -40 5]);
hold off;
