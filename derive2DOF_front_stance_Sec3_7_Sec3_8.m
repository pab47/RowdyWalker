%derive2DOF_front_stance_Sec3_7_Sec3_8
%Files used to derive equations given in Sec 3.7 and 3.8
clc;
clear all;

syms q1 q2 real %Angles as defined in figures 
syms u1 u2 real %Angular velocity
syms ud1 ud2 real%Angular Acceleration
syms gam g real %slope of ramp,, gravity
syms c1 w1 c2 w2 r real % Distances as defined in figures
syms v1 v2 r1 r2 real 
syms m1 m2 real%m1 is body and front fixed foot; m2 is moving leg
syms I1 I2 real %inertia about COM for front and rear leg.
syms xp1 real

i=[1 0 0];
j=[0 1 0];
k=[0 0 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Reference Frames                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X1 = sin(q1)*i-cos(q1)*j;
Y1 = cos(q1)*i+sin(q1)*j;

X2 = sin(q1-q2)*i - cos(q1-q2)*j; 
Y2 = cos(q1-q2)*i + sin(q1-q2)*j;

J  = -sin(gam)*i+cos(gam)*j; %Fictitious ramp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Position Vectors                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Original position vectors (see Fig 6 a)
r_H_G1 = c1*X1 - w1*Y1; % m1 is fixed body and this is the stance leg
r_H_G2 = c2*X2 - w2*Y2;  % m2 is the free leg and this is the swing leg
r_P1_H = r*j;

%Derived position vectors
r_P1_G1 = r_P1_H + r_H_G1;
r_P1_G2 = r_P1_H + r_H_G2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Angular Velocities and Accelerations       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

om1 = u1*k;
om2 = (u1-u2)*k;
 
al1 = ud1*k;
al2 = (ud1-ud2)*k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Constraints, Linear Velocities and Accelerations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q = [q1; q2];
u = [u1; u2];
ud = [ud1; ud2];

xh =  -r*q1 + xp1;
yh =   r;
    
    xhdot  = jacobian(xh,q)*u;
    xhddot = jacobian(xhdot,q)*u + jacobian(xhdot,u)*ud;
    yhdot  = jacobian(yh,q)*u;
    yhddot = jacobian(yhdot,q)*u + jacobian(yhdot,u)*ud;    

v_H  = xhdot*i+yhdot*j; 

a_H  = xhddot*i+yhddot*j;

v_G1 = v_H +cross(om1,r_H_G1);
v_G2 = v_H +cross(om2,r_H_G2);

a_G1 = a_H +cross(om1,cross(om1,r_H_G1))+cross(al1,r_H_G1);
a_G2 = a_H +cross(om2,cross(om2,r_H_G2))+cross(al2,r_H_G2);

%%%%% for two legs %%%%%%
M_P1    = dot(m1*g*cross(r_P1_G1,-J)+m2*g*cross(r_P1_G2,-J),k);
Hdot_P1 = dot(I1*al1+m1*cross(r_P1_G1,a_G1)+I2*al2+m2*cross(r_P1_G2,a_G2),k);

M_H    = dot(m2*g*cross(r_H_G2,-J),k); 
Hdot_H = dot(I2*al2+m2*cross(r_H_G2,a_G2),k);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Angular Momentum, (Before and After Heelstrike  % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H_P1p = dot(I1*om1+m1*cross(r_P1_G1,v_G1)+I2*om2+m2*cross(r_P1_G2,v_G2),k);
H_P1n = subs(H_P1p,{u1,u2,q1,q2},{v1,v2,r1,r2});

AM_P1 = H_P1p;
Mh11  =  subs(AM_P1,[u1 u2],[1 0]);
RHSh1 = H_P1n;
disp('Copy paste in rear_strike');
disp('-----------');
disp(['Mh11 = ',char(simplify(Mh11)),';']);
disp(['RHSh1 = ',char(simplify(RHSh1)),';']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   EOM:[M-Hdot=0](Single Stance and Double Stance)  % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AMB_P1 = M_P1-Hdot_P1;
AMB_H  = M_H-Hdot_H;


%%%% two bodies %%%%%%
   RHSs1 = -subs(AMB_P1,[ud1 ud2],[0 0]); 
   Ms11  =  subs(AMB_P1,[ud1 ud2],[1 0]) + RHSs1;
   Ms12  =  subs(AMB_P1,[ud1 ud2],[0 1]) + RHSs1;


   RHSs2 = -subs(AMB_H,[ud1 ud2],[0 0]); 
   Ms21  =  subs(AMB_H,[ud1 ud2],[1 0]) + RHSs2;
   Ms22  =  subs(AMB_H,[ud1 ud2],[0 1]) + RHSs2;
   
disp(' ');
disp('Copy paste in single_stance2_front_stance');
disp('-----------');
disp(['M11 = ',char(simplify(Ms11)),';']);
disp(['M12 = ',char(simplify(Ms12)),';']);
disp(['M21 = ',char(simplify(Ms21)),';']);
disp(['M22 = ',char(simplify(Ms22)),';']);
disp(['RHS1 = ',char(simplify(RHSs1)),';']);
disp(['RHS2 = ',char(simplify(RHSs2)),';']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Energies                                   %  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get required positions of masses
%wrt ramp frame
xH = xh - xp1;
yH = yh;
xG1 = xH+dot(r_H_G1,i);
yG1 = yH+dot(r_H_G1,j);
xG2 = xH+dot(r_H_G2,i);
yG2 = yH+dot(r_H_G2,j);
% 
% %wrt to global frame
Y_H = yH*cos(gam) - xH*sin(gam);
Y_G1 = yG1*cos(gam) - xG1*sin(gam); 
Y_G2 = yG2*cos(gam) - xG2*sin(gam);
X_H = xH*cos(gam) + yH*sin(gam);
X_G1 = xG1*cos(gam) + yG1*sin(gam);
X_G2 = xG2*cos(gam) + yG2*sin(gam);

KE = 0.5*(simplify(I1*u1*u1+m1*dot(v_G1,v_G1) + I2*(u1-u2)*(u1-u2)+m2*dot(v_G2,v_G2)));
PE = simplify(m1*g*Y_G1+m2*g*Y_G2);
disp(' ');
disp('total_energy2_front_stance');
TE = simplify(KE + PE);
DTE = jacobian(TE,q)*u + jacobian(TE,u)*ud;
disp(['TE = ',char(simplify(TE)),';']);
disp(['DTE = ',char(simplify(DTE)),';']);
