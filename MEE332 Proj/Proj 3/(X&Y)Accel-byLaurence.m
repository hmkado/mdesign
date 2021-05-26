% - Accelerations in the X & Y directions for force analysis

clc
clear variables
close all


% Link Lengths (cm):
 r_ae = 16 ;
 r_ac = 29 ;
 r_of = 13 ;
 r_oe = 16 ;
 r_dc = 3  ;
 r_db = 35 ; 
 r_df = 13 ;  
% 
% Symbol Assignment:
 syms t
% 
% Angle Equations:
 tau = 2            ; % (Opening times (s))
 ax  = zeros(1,100) ;
 ay  = (-17.5-12.5.*cos((pi.*t)./tau))+30 ;
 ang_a = (acos(1 - (((30-ay).^2)./(2.*(r_ae).^2))))./2 ;
%     
% *Position Equations:
 ex = r_ae.*cos(ang_a) ;
 ey = 30 - r_ae.*sin(ang_a) ;
 fx = r_of.*cos(ang_a) ;
 fy = 30 - r_of.*sin(ang_a) ;
 dx = r_of.*cos(ang_a) + r_df.*cos(ang_a) ;
 dy = fy + r_df.*sin(ang_a) ;
 cx = dx + r_dc.*cos(ang_a) ;
 cy = dy - r_dc.*sin(ang_a) ;
 bx = dx + r_db.*cos(ang_a) ;
 by = dy - r_db.*sin(ang_a) ;
% 
% _Velocity Values (x-dir):
 vex = diff(ex) ;
 vfx = diff(fx) ;
 vdx = diff(dx) ;
 vcx = diff(cx) ;
 vbx = diff(bx) ;
 vax = diff(ax) ;
% 
% Velocity Values (y-dir):
 vay = diff(ay) ;
 vby = diff(by) ;
 vcy = diff(cy) ;
 vdy = diff(dy) ;
 vey = diff(ey) ;
 vfy = diff(fy) ;
% 
% Acceleration Values (x-dir):
 aex = diff(vex) ;
 afx = diff(vfx) ;
 adx = diff(vdx) ;
 acx = diff(vcx) ;
 abx = diff(vbx) ;
% 
% Acceleration Values (y-dir):
 aay = diff(vay) ;
 aby = diff(vby) ;
 acy = diff(vcy) ;
 ady = diff(vdy) ;
 aey = diff(vey) ;
 afy = diff(vfy) ;
% 
% Evaluating Position Equations:
 t     = linspace(0,2,100) ; % time array 
 ay    = eval(subs(ay)) ;
 ang_a = eval(subs(ang_a)) ;
 ex    = eval(subs(ex)) ;
 ey    = eval(subs(ey)) ;
 fx    = eval(subs(vex)) ;
 fy    = eval(subs(fy)) ;
 dx    = eval(subs(dx)) ;
 dy    = eval(subs(dy)).*ones(1,100) ; 
 cx    = eval(subs(cx)) ;
 ay    = eval(subs(cy)) ;
 bx    = eval(subs(bx)) ;
 by    = eval(subs(by)) ;
% 
% Evaluating Velocity Equations:
 vax = eval(subs(vax)) ;
 vay = eval(subs(vay)) ;
 vbx = eval(subs(vbx)) ;
 vby = eval(subs(vby)) ;
 vcx = eval(subs(vcx)) ;
 vcy = eval(subs(vcy)) ;
 vdx = eval(subs(vdx)) ;
 vdy = eval(subs(vdy)) ;
 vex = eval(subs(vex)) ;
 vey = eval(subs(vey)) ;
 vfx = eval(subs(vfx)) ;
 vfy = eval(subs(vfy)) ;
% 
% Evaluating Acceleration Equations:
 aay = eval(subs(aay)) ;
 abx = eval(subs(abx)) ;
 aby = eval(subs(aby)) ;
 acx = eval(subs(acx)) ;
 acy = eval(subs(acy)) ;
 adx = eval(subs(adx)) ;
 ady = eval(subs(ady)) ;
 aex = eval(subs(aex)) ;
 aey = eval(subs(aey)) ;
 afx = eval(subs(afx)) ;
 afy = eval(subs(afy)) ;
 
 %linear acceleration plot x direction
figure(1);
grid on;
hold on;
plot(t,abx,'-');
plot(t,acx,'-');
plot(t,adx,'-');
plot(t,aex,'-');
plot(t,afx,'-');
title("Linear Acceleration of Joints x direction");
xlabel('time (s)') ;
ylabel('Acceleration (cm/s^2)') ;
legend('B','C','D','E','F');
hold off;

%linear acceleration plot y direction
figure(2);
grid on;
hold on;
plot(t,aay,'-');
plot(t,aby,'-');
plot(t,acy,'-');
plot(t,ady,'-');
plot(t,aey,'-');
plot(t,afy,'-');
title("Linear Acceleration of Joints y direction");
xlabel('time (s)') ;
ylabel('Acceleration (cm/s^2)') ;
legend('A','B','C','D','E','F');
hold off;
