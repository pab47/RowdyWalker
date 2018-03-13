function RowdyWalkerSim
% MATLAB files for the paper: 
% A 3D Printed, Non-Assembly, Passive Dynamic Walking Toy: Design and Analysis
% By Christian Trevino, Joseph Galloway, Pranav A. Bhounsule
% Submitted to ASME Journal of Mechanisms and Robotics
% More information on this link: tiny.cc/pranavb
%

clc
clear all
close all

format long

robot.f = 0.9; %fraction in the body           
robot.r =  5.5/100; %fixed                     
robot.c1 = 0.35*robot.r; %along length.        
robot.w1 = 0.75/100; %below hip and backward     
robot.c2 = 0.7*robot.r;                       
robot.w2 = robot.w1*.95;                          
robot.gam = atan(2.1/15.2);                  
robot.front_stop = -1.2/5.5;                 
robot.I1 = 0.00013;
robot.I2 = 0.0000030; 

robot.m1 = robot.f*(157/1000);
robot.m2 = (1-robot.f)*(157/1000);
robot.g = 9.81;
robot.fps = 30;
robot.C1 = 0.005; %damping friction in the rear leg

robot.width = 0.15;
robot.front2rear = robot.gam;
robot.rear_stop = 0;
robot.rear2front = robot.gam;
robot.Tol = 1e-12; %integration tolerance
robot.Fsolve_Tol = 1e-10; %fsolve tolerance
robot.animate_x_adjust = 0.007;
robot.animate_y_adjust = -0.007;
Fontsize = 12; %for font size on plots

%%%%%% options that can be set %%%
robot.makeAVI = 0; %writes the file rowdy_walker_sim.avi 
robot.animate_face = 1; %put patch object on the face
robot.animate_border = 0; %animates only the border of rowdy
robot.verbose = 0; %set to 1 to print results
steps = 4; %number of steps


%%% guess for the initial state (see Fig 5. Eqn 1 in the attached paper).
q1 = robot.gam; 
u1 = -1;
z0 = [q1 u1];
%%%% fixed point %%%
%z0=[-0.019156109860067  -1.900074694074746]; %fsolve gives this root

fixed_points = 1;
if (fixed_points)
    %%%% Root finding, Period one gait %%%%
    options = optimset('TolFun',robot.Fsolve_Tol,'TolX',robot.Fsolve_Tol,'Display','off');
    [zstar,fval,exitflag,output,jacob] = fsolve(@fixedpt,z0,options,robot);
    if exitflag == 1
        disp('Fixed points are');
        %fval
        z0=zstar;
        disp(z0);
        %fixed point z0 = [-0.019467427570676  -1.897577473434229];
    else
        error('Root finder not converged, change guess or change system parameters')
    end

    %%% Stability, using linearised eigenvalue %%%
    J=partialder(@onestep,zstar,robot);
    % disp('Jacobian is ');
    % J
    disp(' ');
    disp('EigenValues for linearized map are');
    disp(eig(J))

end

%%%%% Get data for all the steps %%%
[z,t,zz,TE] = onestep(z0,robot,steps);

Distance_sim = zz(end,5)-zz(1,5);
Time_sim = t(end);
Velocity_sim = Distance_sim/Time_sim;
disp('Parameter:  Values  ');
disp('-------------------------------------');
disp(['Distance:   ',num2str(Distance_sim)]);
disp(['Time:   ',num2str(Time_sim)]);
disp(['Velocity:   ',num2str(Velocity_sim)]);

%%%% Animate result %%%
figure(1)
animate(t,zz,robot);

%%% Plot data %%%
h=figure(2);
subplot(3,1,1);
plot(t,zz(:,1),'r','LineWidth',2); ylabel('absolute angle (rad)'); hold on;
plot(t,zz(:,2),'b','LineWidth',1); 
ylabel('Absolute angle (rad)','FontSize',Fontsize);
legend('front','rear','Location','Best');
grid on

subplot(3,1,2);
plot(t,zz(:,3),'r','LineWidth',2); ylabel('absolute rate (rad/s)'); hold on
plot(t,zz(:,4),'b','LineWidth',1); 
ylabel('Absolute angular rate (rad/s)','FontSize',Fontsize);
xlabel('Time (s)','FontSize',Fontsize);
legend('front','rear','Location','Best');
grid on
%%% NOTE: absolute value + CCW and - is CW.
%h = figure(3)
subplot(3,1,3);
plot(t,TE,'k');
xlabel('time');
ylabel('Total Energy');
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FUNCTIONS START HERE %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================================
function zdiff=fixedpt(z0,robot)
%===================================================================
zdiff=onestep(z0,robot)-z0; 

%===================================================================
function J=partialder(FUN,z,robot)
%===================================================================
pert=1e-5;

%%%% Using forward difference, accuracy linear %%%
% y0=feval(FUN,z,GL_DIM); 
% for i=1:length(z)
%     ztemp=z;
%     ztemp(i)=ztemp(i)+pert; 
%     J(:,i)=(feval(FUN,ztemp,GL_DIM)-y0) ;
% end
% J=(J/pert);

%%% Using central difference, accuracy quadratic %%%
for i=1:length(z)
    ztemp1=z; ztemp2=z;
    ztemp1(i)=ztemp1(i)+pert; 
    ztemp2(i)=ztemp2(i)-pert; 
    J(:,i)=(feval(FUN,ztemp1,robot)-feval(FUN,ztemp2,robot)) ;
end
J=J/(2*pert);

%===================================================================
function [z_temp,t_temp,zz,TE] = onestep(z0,robot,steps)
%===================================================================

flag = 1;
if nargin<2
    error('need more inputs to onestep');
elseif nargin<3
    flag = 0; %send only last state, for root finder and jacobian
    steps = 1;
end

q1 = z0(1);
u1 = z0(2);
q2 = 0;
u2 = 0;

r = robot.r;
xh = 0;%-r*q1 + xp1;
vxh = -r*u1; 
yh =  r;
vyh = 0; 

t0 = 0;
z0 = [q1 u1 q2 u2 xh vxh yh vyh];

E0 = total_energy1_front_stance(t0,z0,robot);
TE = E0;
t_temp = t0; 
z_temp = [z0 E0]; 
zz = [z0(1) z0(1)-z0(3) z0(2) z0(2)-z0(4) z0(5) z0(7)];
if (robot.verbose && nargin==3)
    disp(['start state is z = ',num2str(zz(1:4))]);
end

for i=1:steps
    %%% starts in leg together and front tilted forward %%%%
    %%%%%% one dof, front 2 rear %%%%%%%
    E0 = total_energy1_front_stance(t0,z0,robot);
    z0 = [z0,E0];
    t0 = t_temp(end,:);
    dt = 3;
    options=odeset('abstol',robot.Tol,'reltol',robot.Tol,'events',@front2rear);
    tspan = linspace(t0,t0+dt,1000);
    [t_temp1, z_temp1] = ode113(@single_stance1_front,tspan,z0,options,robot);
    zz1 = [z_temp1(:,1) z_temp1(:,1)-z_temp1(:,3) ...
           z_temp1(:,2) z_temp1(:,2)-z_temp1(:,4)...
           z_temp1(:,5) z_temp1(:,7)];

    t_temp = [t_temp; t_temp1(2:end)];
    if (robot.verbose && nargin==3)
        disp(['front2rear, t = ',num2str(t_temp1(end))]);
        disp(['Absolute state, zz = ',num2str(zz1(end,1:4))]);
    end
    
    z_temp = [z_temp; z_temp1(2:end,:)];
    TE = [TE; z_temp1(2:end,9)];
    
    zz = [zz; zz1(2:end,:)];

    %%%%% two dof, rear support %%%%%%
    q1 = z_temp1(end,1)-z_temp1(end,3); 
    u1 = z_temp1(end,2)-z_temp1(end,4);
    q2 = z_temp1(end,3); %=0
    u2 = z_temp1(end,4);%=0
    xh = z_temp1(end,5);
    vxh = z_temp1(end,6); 
    yh =  z_temp1(end,7);
    vyh = z_temp1(end,8);
    z0 = [q1 u1 q2 u2 xh vxh yh vyh];
    t0 = t_temp1(end);
    E0 = total_energy2_rear_stance(t0,z0,robot);
    z0 = [z0 E0];

    dt = 3;
    options=odeset('abstol',robot.Tol,'reltol',robot.Tol,'events',@front_stop);
    tspan = linspace(t0,t0+dt,1000);
    [t_temp2, z_temp2] = ode113(@single_stance2_rear_stance,tspan,z0,options,robot);
    zz2 = [z_temp2(:,1)-z_temp2(:,3) z_temp2(:,1) ...
           z_temp2(:,2)-z_temp2(:,4) z_temp2(:,2) ...
           z_temp2(:,5) z_temp2(:,7)]; 
    t_temp = [t_temp; t_temp2(2:end)];
    z_temp = [z_temp; z_temp2(2:end,:)];
    zz = [zz; zz2(2:end,:)];
    TE = [TE; z_temp2(2:end,9)];
    if (robot.verbose && nargin==3)
        disp(['front_stop, t = ',num2str(t_temp2(end))]);
        disp(['Absolute state, zz = ',num2str(zz2(end,1:4))]);
    end
 
   %%%%%%% one dof, rear 2 front %%
    r1 = z_temp2(end,1); 
    v1 = z_temp2(end,2);
    r2 = z_temp2(end,3);
    v2 = z_temp2(end,4);

    [q1,u1,q2,u2] = front_strike(r1,v1,r2,v2,robot);
    
    xh = z_temp2(end,5);
    vxh = -r*u1; 
    yh =  z_temp2(end,7);
    vyh = z_temp2(end,8);
    z0 = [q1 u1 q2 u2 xh vxh yh vyh]; 
    t0 = t_temp(end,:);
    E0 = total_energy2_rear_stance(t0,z0,robot);
    
    z0 = [z0 E0];
    dt = 3;
    options=odeset('abstol',robot.Tol,'reltol',robot.Tol,'events',@rear2front);
    tspan = linspace(t0,t0+dt,1000);
    [t_temp3, z_temp3] = ode113(@single_stance1_rear,tspan,z0,options,robot);
    zz3 = [z_temp3(:,1)-z_temp3(:,3) z_temp3(:,1) ...
           z_temp3(:,2)-z_temp3(:,4) z_temp3(:,2) ...
           z_temp3(:,5) z_temp3(:,7)];
    t_temp = [t_temp; t_temp3(2:end)];
    z_temp = [z_temp; z_temp3(2:end,:)];
    zz = [zz; zz3(2:end,:)];

    TE = [TE; z_temp3(2:end,9)];
    if (robot.verbose && nargin==3)
        disp(['rear2front, t = ',num2str(t_temp3(end))]);
        disp(['Absolute state, zz = ',num2str(zz3(end,1:4))]);
    end
    
    %%%%%% two dof, front support %%%%%%%
    q1 = z_temp3(end,1)-z_temp3(end,3);  
    u1 = z_temp3(end,2);
    q2 = -z_temp3(end,3);
    u2 = 0;
    xh = z_temp3(end,5);
    vxh = z_temp3(end,6); 
    yh =  z_temp3(end,7);
    vyh = z_temp3(end,8);
    z0 = [q1 u1 q2 u2 xh vxh yh vyh];
    E0 = total_energy2_front_stance(t0,z0,robot);
    z0 = [z0 E0];
    t0 = t_temp(end,:);
    dt = 3;
    options=odeset('abstol',robot.Tol,'reltol',robot.Tol,'events',@rear_stop);
    tspan = linspace(t0,t0+dt,1000);
    [t_temp4, z_temp4] = ode113(@single_stance2_front_stance,tspan,z0,options,robot);
    zz4 = [z_temp4(:,1) z_temp4(:,1)-z_temp4(:,3) ...
           z_temp4(:,2) z_temp4(:,2)-z_temp4(:,4) ...
           z_temp4(:,5) z_temp4(:,7)];
    t_temp = [t_temp; t_temp4(2:end)];
    z_temp = [z_temp; z_temp4(2:end,:)];
    zz = [zz; zz4(2:end,:)];
     TE = [TE; z_temp4(2:end,9)];
    if (robot.verbose && nargin==3)
        disp(['rear_stop, t = ',num2str(t_temp4(end))]);
        disp(['Absolute state, zz = ',num2str(zz4(end,1:4))]);
    end
    
    %%%%%%% init condition for front rocking %%%
    r1 = z_temp4(end,1); 
    v1 = z_temp4(end,2);
    r2 = z_temp4(end,3);
    v2 = z_temp4(end,4);
    [q1,u1,q2,u2] = rear_strike(r1,v1,r2,v2,robot);
      
    xh = z_temp4(end,5);
    vxh = -r*u1; 
    yh =  z_temp4(end,7);
    vyh = z_temp4(end,8);
    z0 = [q1 u1 q2 u2 xh vxh yh vyh];

end

if flag==0
   z_temp=[q1 u1];
   t_temp = [];
   zz = [];
end

%===================================================================
function animate(t_all,z_all,robot)
%===================================================================
fps = robot.fps;
r = robot.r;
m1 = robot.m1; c1 = robot.c1; w1 = robot.w1; 
m2 = robot.m2; c2 = robot.c2; w2 = robot.w2; 
gam = robot.gam;

[mm,nn] = size(z_all);
t = linspace(0,t_all(end),fps*ceil(t_all(end)));
for i=1:nn
    z(:,i) = interp1(t_all,z_all(:,i),t);
end

if (robot.animate_face)
    [X,Y] = rowdy_data;
    scale = robot.width/(max(X) - min(X));
    X = robot.animate_x_adjust+scale*X; Y = robot.animate_y_adjust+scale*Y;
    N = length(X);
end

qq1 = z(:,1); qq2 = z(:,2);
xxh = z(:,5); yyh = z(:,6);

if (robot.makeAVI)
    mov = VideoWriter('rowdy_walker_sim.avi'); 
    open(mov);
end

for i=1:length(t)
    q1 = qq1(i); q2 = qq2(i);
    xh = xxh(i); yh = yyh(i);
        
    x_H_com = (m1*(c1*sin(q1) - w1*cos(q1)) + m2*(c2*sin(q1 - q2) - w2*cos(q1 - q2)))/(m1 + m2);
    y_H_com = -(m1*(c1*cos(q1) + w1*sin(q1)) + m2*(c2*cos(q1 - q2) + w2*sin(q1 - q2)))/(m1 + m2);
    [x_H_com,y_H_com] = world_frame(gam,x_H_com,y_H_com);
     
    com_x = xh+x_H_com; 
    com_y = yh+y_H_com;
    
    
    [hinge_x,hinge_y] = world_frame(gam,xh,yh);
    plot(hinge_x,hinge_y,'ro','Markersize',5,'Markerfacecolor','g'); hold on;
    
     [com_x,com_y] = world_frame(gam,com_x,com_y);
     plot(com_x,com_y,'ko','Markersize',20,'Markerfacecolor','k');
     
    ramp_x_min = -1; ramp_x_max = 1;
    ramp_y_min = 0; ramp_y_max = 0;
    [ramp_x_min,ramp_y_min] = world_frame(gam,ramp_x_min,ramp_y_min);
    [ramp_x_max,ramp_y_max] = world_frame(gam,ramp_x_max,ramp_y_max);
    line([ramp_x_min ramp_x_max],[ramp_y_min ramp_y_max],'Linewidth',2,'Color','k'); %ramp
     
    center_x = xh; center_y = yh;
    tip_x = xh+r*sin(q1); tip_y = yh-r*cos(q1);
    tip_x2 = xh+r*sin(q2); tip_y2 = yh-r*cos(q2); 
    
    [tip_x,tip_y] = world_frame(gam,tip_x,tip_y);
    [tip_x2,tip_y2] = world_frame(gam,tip_x2,tip_y2);

    th = 0.013/robot.r;
    curve_x_fwd = center_x + r*sin(q1+th); curve_y_fwd = center_y - r*cos(q1+th);
    curve_x2_fwd = center_x + r*sin(q1+2*th); curve_y2_fwd = center_y - r*cos(q1+2*th);
    [curve_x_fwd, curve_y_fwd] = world_frame(gam,curve_x_fwd, curve_y_fwd);
    [curve_x2_fwd, curve_y2_fwd] = world_frame(gam,curve_x2_fwd, curve_y2_fwd);

    th = -0.013/robot.r;
    curve_x_bwd = center_x + r*sin(q2+th); curve_y_bwd = center_y - r*cos(q2+th);
    curve_x2_bwd = center_x + r*sin(q2+2*th); curve_y2_bwd = center_y - r*cos(q2+2*th);
    [curve_x_bwd, curve_y_bwd] = world_frame(gam,curve_x_bwd, curve_y_bwd);
    [curve_x2_bwd, curve_y2_bwd] = world_frame(gam,curve_x2_bwd, curve_y2_bwd);
    
    if (robot.animate_border)
           patch([hinge_x tip_x curve_x_fwd curve_x2_fwd ],...
          [hinge_y tip_y curve_y_fwd curve_y2_fwd ],...
          [255,255,255]/255,'EdgeColor','k'); %front leg
          patch([hinge_x tip_x2 curve_x_bwd curve_x2_bwd ],...
          [hinge_y tip_y2 curve_y_bwd curve_y2_bwd ],...
          [255,255,255]/255,'EdgeColor','k');  
          %[255,140,0]/255,'EdgeColor',[255,140,0]/255);       
    else
    patch([hinge_x tip_x curve_x_fwd curve_x2_fwd ],...
          [hinge_y tip_y curve_y_fwd curve_y2_fwd ],...
          [241, 90, 34]/255,'EdgeColor',[241, 90, 34]/255); %front leg
    patch([hinge_x tip_x2 curve_x_bwd curve_x2_bwd ],...
          [hinge_y tip_y2 curve_y_bwd curve_y2_bwd ],...
          [255,140,0]/255,'EdgeColor',[255,140,0]/255);  
          %[255,140,0]/255,'EdgeColor',[255,140,0]/255);  
    end
    
    if (robot.animate_face)  
        R = rotation(q1);
        for ii=1:N
            pos_face = [R(1,1); R(1,2)]*X(ii) +  [R(2,1); R(2,2)]*Y(ii);
            x_face(ii) = xh + pos_face(1); y_face(ii) = yh+pos_face(2);
            [x_face(ii),y_face(ii)] = world_frame(gam,x_face(ii),y_face(ii));
        end
        if (robot.animate_border)
            patch(x_face,y_face,[255,255,255]/255,'EdgeColor','k'); 
        else
            patch(x_face,y_face,[241, 90, 34]/255,'EdgeColor',[255,255,255]/255); 
        end
    end
    
%    %%%%% draw hinge, tip_x, and tip_x2
     plot(hinge_x,hinge_y,'ko','Markersize',3,'MarkerFacecolor','k');
     plot(tip_x,tip_y,'ko','Markersize',3,'MarkerFacecolor','k');
     plot(tip_x2,tip_y2,'ko','Markersize',3,'MarkerFacecolor','k');
    
    axis('equal'); 
    axis([-0.08 0.21 -0.05 0.1]); 
    if (robot.makeAVI)
        axis off
        set(gcf,'Color',[1,1,1])
        currFrame = getframe;
        writeVideo(mov,currFrame);
    end
    pause(0.05);
    hold off;
    
end

if (robot.makeAVI)
    close(mov);
end

%===================================================================
function R = rotation(theta)
%===================================================================
c = cos(theta); s = sin(theta);
R = [c s; -s c];

%===================================================================
function [x,y] = world_frame(theta,x,y)
%===================================================================
c = cos(theta); s = sin(theta);
R = [c s; -s c];
orig_vec = [x;y];
new_vec = R*orig_vec;
x = new_vec(1); y = new_vec(2);


%===================================================================
function [gstop, isterminal,direction]=front2rear(t,z,robot)
%===================================================================
  
q1 = z(1); 

gstop = q1 - robot.front2rear;% gam;
isterminal = 1; %terminate integration
%q1 goes from - to 0 to +
%q1-gam goes from - to +
direction = 1;

%===================================================================
function [gstop, isterminal,direction]=rear2front(t,z,robot)
%===================================================================
  
q1 = z(1); 
q2 = z(3);

gstop = (q1-q2) - robot.rear2front;% gam;
isterminal = 1; %terminate integration
% % q1 goes from + to 0 to -
% % q1-gam goes from + to -
direction = -1;

%===================================================================
function [gstop, isterminal,direction]=front_stop(t,z,robot)
%===================================================================
 
q2 = z(3); 

gstop = q2 - robot.front_stop;%
isterminal = 1; %terminate integration
%q3 goes from 0 to -
%q1-robot.rear2front goes from + to -
direction = -1;

%===================================================================
function [gstop, isterminal,direction]=rear_stop(t,z,robot)
%===================================================================
 
q2 = z(3); 

gstop = q2 - robot.rear_stop;%
isterminal = 1; %terminate integration
%q2 goes from 0 to -
%q1-robot.rear2front goes from + to -
direction = -1;


%===================================================================
function zdot=single_stance1_front(t,z,robot)   %PAB (DONE)
%===================================================================

q1 = z(1);   u1 = z(2);                         
q2 = z(3);   u2 = z(4);                         
xh = z(5);  vxh = z(6);                       
yh = z(7);  vyh = z(8);                     

m1 = robot.m1;  m2 = robot.m2; g = robot.g;   gam = robot.gam;  
r = robot.r;    c2 = robot.c2; w2 = robot.w2; c1 = robot.c1; w1 = robot.w1;  
I1 = robot.I1; I2 = robot.I2;

MM = 2*c1*m1*r*cos(q1) - I2 - c1^2*m1 - c2^2*m2 - m1*r^2 - m2*r^2 - m1*w1^2 - m2*w2^2 - I1 + 2*m1*r*w1*sin(q1) + 2*c2*m2*r*cos(q1 - q2) + 2*m2*r*w2*sin(q1 - q2);
RHS = g*m1*r*sin(gam) - g*m2*w2*cos(gam - q1 + q2) - c2*g*m2*sin(gam - q1 + q2) + g*m2*r*sin(gam) - c1*g*m1*sin(gam - q1) - g*m1*w1*cos(gam - q1) + c1*m1*r*u1^2*sin(q1) - m1*r*u1^2*w1*cos(q1) + c2*m2*r*u1^2*sin(q1 - q2) - m2*r*u1^2*w2*cos(q1 - q2);

X = MM \ RHS;                                    

ud1 = X(1);                                       
ud2 = 0;                                       

axh = -r*ud1; 
ayh = 0; 

DTE = u1*(I1*ud1 + I2*ud1 + c1^2*m1*ud1 + c2^2*m2*ud1 + m1*r^2*ud1 + m2*r^2*ud1 + m1*ud1*w1^2 + m2*ud1*w2^2 - c2*g*m2*sin(gam - q1 + q2) - g*m2*w2*cos(gam - q1 + q2) + g*m1*r*sin(gam) + g*m2*r*sin(gam) - c1*g*m1*sin(gam - q1) - g*m1*w1*cos(gam - q1) - 2*m1*r*ud1*w1*sin(q1) - 2*c2*m2*r*ud1*cos(q1 - q2) - 2*m2*r*ud1*w2*sin(q1 - q2) + c1*m1*r*u1^2*sin(q1) - m1*r*u1^2*w1*cos(q1) + c2*m2*r*u1^2*sin(q1 - q2) - m2*r*u1^2*w2*cos(q1 - q2) - 2*c1*m1*r*ud1*cos(q1));


zdot = [u1 ud1 u2 ud2 ...          
        vxh axh vyh ayh DTE]';                   

    
%===================================================================
function zdot=single_stance1_rear(t,z,robot)   %PAB (DONE)
%===================================================================

q1 = z(1);   u1 = z(2);                         
q2 = z(3);   u2 = z(4);                         
xh = z(5);  vxh = z(6);                       
yh = z(7);  vyh = z(8);                     

m1 = robot.m1;  m2 = robot.m2; g = robot.g;   gam = robot.gam;  
r = robot.r;    c2 = robot.c2; w2 = robot.w2; c1 = robot.c1; w1 = robot.w1;   
I1 = robot.I1; I2 = robot.I2;

MM = 2*c2*m2*r*cos(q1) - I2 - c1^2*m1 - c2^2*m2 - m1*r^2 - m2*r^2 - m1*w1^2 - m2*w2^2 - I1 + 2*m2*r*w2*sin(q1) + 2*c1*m1*r*cos(q1 - q2) + 2*m1*r*w1*sin(q1 - q2);
RHS = g*m1*r*sin(gam) - g*m1*w1*cos(gam - q1 + q2) - c1*g*m1*sin(gam - q1 + q2) + g*m2*r*sin(gam) - c2*g*m2*sin(gam - q1) - g*m2*w2*cos(gam - q1) + c2*m2*r*u1^2*sin(q1) - m2*r*u1^2*w2*cos(q1) + c1*m1*r*u1^2*sin(q1 - q2) - m1*r*u1^2*w1*cos(q1 - q2);

X = MM \ RHS;                                    

ud1 = X(1);                                       
ud2 = 0;                                       

axh = -r*ud1; 
ayh = 0; 

DTE = u1*(I1*ud1 + I2*ud1 + c1^2*m1*ud1 + c2^2*m2*ud1 + m1*r^2*ud1 + m2*r^2*ud1 + m1*ud1*w1^2 + m2*ud1*w2^2 - c1*g*m1*sin(gam - q1 + q2) - g*m1*w1*cos(gam - q1 + q2) + g*m1*r*sin(gam) + g*m2*r*sin(gam) - c2*g*m2*sin(gam - q1) - g*m2*w2*cos(gam - q1) - 2*m2*r*ud1*w2*sin(q1) - 2*c1*m1*r*ud1*cos(q1 - q2) - 2*m1*r*ud1*w1*sin(q1 - q2) + c2*m2*r*u1^2*sin(q1) - m2*r*u1^2*w2*cos(q1) + c1*m1*r*u1^2*sin(q1 - q2) - m1*r*u1^2*w1*cos(q1 - q2) - 2*c2*m2*r*ud1*cos(q1));


zdot = [u1 ud1 u2 ud2 ...          
        vxh axh vyh ayh DTE]';                   

%===================================================================
function zdot=single_stance2_rear_stance(t,z,robot) %PAB (DONE)
%===================================================================

q1 = z(1);   u1 = z(2);                         
q2 = z(3);   u2 = z(4);                         
xh = z(5);  vxh = z(6);                       
yh = z(7);  vyh = z(8);                     

m1 = robot.m1;  m2 = robot.m2; g = robot.g;   gam = robot.gam;  
r = robot.r;    c2 = robot.c2; w2 = robot.w2; c1 = robot.c1; w1 = robot.w1; 
C1 = robot.C1;
I1 = robot.I1; I2 = robot.I2;


M11 = 2*c2*m2*r*cos(q1) - I2 - c1^2*m1 - c2^2*m2 - m1*r^2 - m2*r^2 - m1*w1^2 - m2*w2^2 - I1 + 2*m2*r*w2*sin(q1) + 2*c1*m1*r*cos(q1 - q2) + 2*m1*r*w1*sin(q1 - q2);
M12 = I1 + c1^2*m1 + m1*w1^2 - c1*m1*r*cos(q1 - q2) - m1*r*w1*sin(q1 - q2);
M21 = c1*m1*r*cos(q1 - q2) - c1^2*m1 - m1*w1^2 - I1 + m1*r*w1*sin(q1 - q2);
M22 = I1 + c1^2*m1 + m1*w1^2;
RHS1 = C1*u1-m1*((c1*sin(q1 - q2) - w1*cos(q1 - q2))*(u1 - u2)^2*(c1*cos(q1 - q2) - r + w1*sin(q1 - q2)) - (c1*cos(q1 - q2) + w1*sin(q1 - q2))*(c1*sin(q1 - q2) - w1*cos(q1 - q2))*(u1 - u2)^2) - m2*(u1^2*(c2*sin(q1) - w2*cos(q1))*(c2*cos(q1) - r + w2*sin(q1)) - u1^2*(c2*cos(q1) + w2*sin(q1))*(c2*sin(q1) - w2*cos(q1))) - g*m2*(sin(gam)*(c2*cos(q1) - r + w2*sin(q1)) - cos(gam)*(c2*sin(q1) - w2*cos(q1))) - g*m1*(sin(gam)*(c1*cos(q1 - q2) - r + w1*sin(q1 - q2)) - cos(gam)*(c1*sin(q1 - q2) - w1*cos(q1 - q2)));
RHS2 = 0*u2+g*m1*(cos(gam)*(c1*sin(q1 - q2) - w1*cos(q1 - q2)) - sin(gam)*(c1*cos(q1 - q2) + w1*sin(q1 - q2)));
 
MM = [M11 M12;                               
     M21 M22];                               

RHS = [RHS1; RHS2];                       

X = MM \ RHS;                                    

ud1 = X(1);                                       
ud2 = X(2);                                       

axh = -r*ud1; 
ayh = 0; 

DTE = I1*u1*ud1 - I1*u1*ud2 - I1*u2*ud1 + I2*u1*ud1 + I1*u2*ud2 + c1^2*m1*u1*ud1 - c1^2*m1*u1*ud2 - c1^2*m1*u2*ud1 + c1^2*m1*u2*ud2 + c2^2*m2*u1*ud1 + m1*r^2*u1*ud1 + m2*r^2*u1*ud1 + m1*u1*ud1*w1^2 - m1*u1*ud2*w1^2 - m1*u2*ud1*w1^2 + m1*u2*ud2*w1^2 + m2*u1*ud1*w2^2 - c2*g*m2*u1*sin(gam - q1) - g*m2*u1*w2*cos(gam - q1) + c2*m2*r*u1^3*sin(q1) - m2*r*u1^3*w2*cos(q1) - c1*g*m1*u1*sin(gam - q1 + q2) + c1*g*m1*u2*sin(gam - q1 + q2) + c1*m1*r*u1^3*sin(q1 - q2) - g*m1*u1*w1*cos(gam - q1 + q2) + g*m1*u2*w1*cos(gam - q1 + q2) - m1*r*u1^3*w1*cos(q1 - q2) + g*m1*r*u1*sin(gam) + g*m2*r*u1*sin(gam) - m1*r*u1*u2^2*w1*cos(q1 - q2) + 2*m1*r*u1^2*u2*w1*cos(q1 - q2) - 2*c2*m2*r*u1*ud1*cos(q1) - 2*m2*r*u1*ud1*w2*sin(q1) - 2*c1*m1*r*u1*ud1*cos(q1 - q2) + c1*m1*r*u1*ud2*cos(q1 - q2) + c1*m1*r*u2*ud1*cos(q1 - q2) - 2*m1*r*u1*ud1*w1*sin(q1 - q2) + m1*r*u1*ud2*w1*sin(q1 - q2) + m1*r*u2*ud1*w1*sin(q1 - q2) + c1*m1*r*u1*u2^2*sin(q1 - q2) - 2*c1*m1*r*u1^2*u2*sin(q1 - q2);

zdot = [u1 ud1 u2 ud2 ...           
        vxh axh vyh ayh DTE]';                   
 
    
%===================================================================
function zdot=single_stance2_front_stance(t,z,robot)  %PAB (DONE)
%===================================================================

q1 = z(1);   u1 = z(2);                         
q2 = z(3);   u2 = z(4);                         
xh = z(5);  vxh = z(6);                       
yh = z(7);  vyh = z(8);                     

m1 = robot.m1;  m2 = robot.m2; g = robot.g;   gam = robot.gam;  
r = robot.r;    c2 = robot.c2; w2 = robot.w2; c1 = robot.c1; w1 = robot.w1;   
I1 = robot.I1; I2 = robot.I2;

M11 = 2*c1*m1*r*cos(q1) - I2 - c1^2*m1 - c2^2*m2 - m1*r^2 - m2*r^2 - m1*w1^2 - m2*w2^2 - I1 + 2*m1*r*w1*sin(q1) + 2*c2*m2*r*cos(q1 - q2) + 2*m2*r*w2*sin(q1 - q2);
M12 = I2 + c2^2*m2 + m2*w2^2 - c2*m2*r*cos(q1 - q2) - m2*r*w2*sin(q1 - q2);
M21 = c2*m2*r*cos(q1 - q2) - c2^2*m2 - m2*w2^2 - I2 + m2*r*w2*sin(q1 - q2);
M22 = I2 + c2^2*m2 + m2*w2^2;
RHS1 = - m2*((c2*sin(q1 - q2) - w2*cos(q1 - q2))*(u1 - u2)^2*(c2*cos(q1 - q2) - r + w2*sin(q1 - q2)) - (c2*cos(q1 - q2) + w2*sin(q1 - q2))*(c2*sin(q1 - q2) - w2*cos(q1 - q2))*(u1 - u2)^2) - m1*(u1^2*(c1*sin(q1) - w1*cos(q1))*(c1*cos(q1) - r + w1*sin(q1)) - u1^2*(c1*cos(q1) + w1*sin(q1))*(c1*sin(q1) - w1*cos(q1))) - g*m1*(sin(gam)*(c1*cos(q1) - r + w1*sin(q1)) - cos(gam)*(c1*sin(q1) - w1*cos(q1))) - g*m2*(sin(gam)*(c2*cos(q1 - q2) - r + w2*sin(q1 - q2)) - cos(gam)*(c2*sin(q1 - q2) - w2*cos(q1 - q2)));
RHS2 = g*m2*(cos(gam)*(c2*sin(q1 - q2) - w2*cos(q1 - q2)) - sin(gam)*(c2*cos(q1 - q2) + w2*sin(q1 - q2)));
  
MM = [M11 M12;                               
     M21 M22];                               

RHS = [RHS1; RHS2];                       

X = MM \ RHS;                                    

ud1 = X(1);                                       
ud2 = X(2);                                       

axh = -r*ud1; 
ayh = 0; 

DTE = I1*u1*ud1 + I2*u1*ud1 - I2*u1*ud2 - I2*u2*ud1 + I2*u2*ud2 + c1^2*m1*u1*ud1 + c2^2*m2*u1*ud1 - c2^2*m2*u1*ud2 - c2^2*m2*u2*ud1 + c2^2*m2*u2*ud2 + m1*r^2*u1*ud1 + m2*r^2*u1*ud1 + m1*u1*ud1*w1^2 + m2*u1*ud1*w2^2 - m2*u1*ud2*w2^2 - m2*u2*ud1*w2^2 + m2*u2*ud2*w2^2 - c1*g*m1*u1*sin(gam - q1) - g*m1*u1*w1*cos(gam - q1) + c1*m1*r*u1^3*sin(q1) - m1*r*u1^3*w1*cos(q1) - c2*g*m2*u1*sin(gam - q1 + q2) + c2*g*m2*u2*sin(gam - q1 + q2) + c2*m2*r*u1^3*sin(q1 - q2) - g*m2*u1*w2*cos(gam - q1 + q2) + g*m2*u2*w2*cos(gam - q1 + q2) - m2*r*u1^3*w2*cos(q1 - q2) + g*m1*r*u1*sin(gam) + g*m2*r*u1*sin(gam) - m2*r*u1*u2^2*w2*cos(q1 - q2) + 2*m2*r*u1^2*u2*w2*cos(q1 - q2) - 2*c1*m1*r*u1*ud1*cos(q1) - 2*m1*r*u1*ud1*w1*sin(q1) - 2*c2*m2*r*u1*ud1*cos(q1 - q2) + c2*m2*r*u1*ud2*cos(q1 - q2) + c2*m2*r*u2*ud1*cos(q1 - q2) - 2*m2*r*u1*ud1*w2*sin(q1 - q2) + m2*r*u1*ud2*w2*sin(q1 - q2) + m2*r*u2*ud1*w2*sin(q1 - q2) + c2*m2*r*u1*u2^2*sin(q1 - q2) - 2*c2*m2*r*u1^2*u2*sin(q1 - q2);

zdot = [u1 ud1 u2 ud2 ...           
        vxh axh vyh ayh DTE]';    
    
%===================================================================
function TE = total_energy1_front_stance(t,z,robot) %PAB
%===================================================================

m1 = robot.m1; c1 = robot.c1; w1 = robot.w1;
m2 = robot.m2; c2 = robot.c2; w2 = robot.w2;
r = robot.r; g = robot.g; gam = robot.gam;
I1 = robot.I1; I2 = robot.I2;

for i=1:length(t)
    q1 = z(i,1); u1 = z(i,2);
    q2 = z(i,3); %u2 = z(i,4);
    TE(i,1) = (I1*u1^2)/2 + (I2*u1^2)/2 + (m2*(u1^2*(c2*sin(q1 - q2) - w2*cos(q1 - q2))^2 + (r*u1 - u1*(c2*cos(q1 - q2) + w2*sin(q1 - q2)))^2))/2 + (m1*(u1^2*(c1*sin(q1) - w1*cos(q1))^2 + (r*u1 - u1*(c1*cos(q1) + w1*sin(q1)))^2))/2 - g*m2*(cos(gam)*(c2*cos(q1 - q2) - r + w2*sin(q1 - q2)) - sin(gam)*(q1*r - c2*sin(q1 - q2) + w2*cos(q1 - q2))) + g*m1*(sin(gam)*(q1*r - c1*sin(q1) + w1*cos(q1)) - cos(gam)*(c1*cos(q1) - r + w1*sin(q1)));
end


%===================================================================
function TE = total_energy2_front_stance(t,z,robot) %PAB (DONE)
%===================================================================

m1 = robot.m1; c1 = robot.c1; w1 = robot.w1;
m2 = robot.m2; c2 = robot.c2; w2 = robot.w2;
r = robot.r; g = robot.g; gam = robot.gam;
I1 = robot.I1; I2 = robot.I2;

for i=1:length(t)
    q1 = z(i,1); u1 = z(i,2);
    q2 = z(i,3); u2 = z(i,4);
    TE(i,1) = (m2*(((c2*cos(q1 - q2) + w2*sin(q1 - q2))*(u1 - u2) - r*u1)^2 + (c2*sin(q1 - q2) - w2*cos(q1 - q2))^2*(u1 - u2)^2))/2 + (I1*u1^2)/2 + (I2*(u1 - u2)^2)/2 + (m1*(u1^2*(c1*sin(q1) - w1*cos(q1))^2 + (r*u1 - u1*(c1*cos(q1) + w1*sin(q1)))^2))/2 - g*m2*(cos(gam)*(c2*cos(q1 - q2) - r + w2*sin(q1 - q2)) - sin(gam)*(q1*r - c2*sin(q1 - q2) + w2*cos(q1 - q2))) + g*m1*(sin(gam)*(q1*r - c1*sin(q1) + w1*cos(q1)) - cos(gam)*(c1*cos(q1) - r + w1*sin(q1)));
end

%===================================================================
function [TE,KE_m1,KE_m2] = total_energy2_rear_stance(t,z,robot) %PAB (DONE)
%===================================================================

m1 = robot.m1; c1 = robot.c1; w1 = robot.w1;
m2 = robot.m2; c2 = robot.c2; w2 = robot.w2;
r = robot.r; g = robot.g; gam = robot.gam;
I1 = robot.I1; I2 = robot.I2;

for i=1:length(t)
    q1 = z(i,1); u1 = z(i,2);
    q2 = z(i,3); u2 = z(i,4);
    TE(i,1) = (m1*(((c1*cos(q1 - q2) + w1*sin(q1 - q2))*(u1 - u2) - r*u1)^2 + (c1*sin(q1 - q2) - w1*cos(q1 - q2))^2*(u1 - u2)^2))/2 + (I2*u1^2)/2 + (I1*(u1 - u2)^2)/2 + (m2*(u1^2*(c2*sin(q1) - w2*cos(q1))^2 + (r*u1 - u1*(c2*cos(q1) + w2*sin(q1)))^2))/2 - g*m1*(cos(gam)*(c1*cos(q1 - q2) - r + w1*sin(q1 - q2)) - sin(gam)*(q1*r - c1*sin(q1 - q2) + w1*cos(q1 - q2))) + g*m2*(sin(gam)*(q1*r - c2*sin(q1) + w2*cos(q1)) - cos(gam)*(c2*cos(q1) - r + w2*sin(q1)));
    
end  

KE_m2(i,1) = (I2*u1^2)/2 + (m2*(u1^2*(c2*sin(q1) - w2*cos(q1))^2 + (r*u1 - u1*(c2*cos(q1) + w2*sin(q1)))^2))/2;
KE_m1(i,1) = (I1*(u1 - u2)^2)/2 + (m1*(((c1*cos(q1 - q2) + w1*sin(q1 - q2))*(u1 - u2) - r*u1)^2 + (c1*sin(q1 - q2) - w1*cos(q1 - q2))^2*(u1 - u2)^2))/2;

%===================================================================
function [q1,u1,q2,u2] = front_strike(r1,v1,r2,v2,robot) %PAB (DONE)
%===================================================================

m1 = robot.m1; c1 = robot.c1; w1 = robot.w1;
m2 = robot.m2; c2 = robot.c2; w2 = robot.w2;
r = robot.r; 
I1 = robot.I1; I2 = robot.I2;
%g = robot.g; gam = robot.gam;

q1 = r1;
q2 = r2;
u2 = 0;

Mh11 = I1 + I2 + c1^2*m1 + c2^2*m2 + m1*r^2 + m2*r^2 + m1*w1^2 + m2*w2^2 - 2*c2*m2*r*cos(q1) - 2*m2*r*w2*sin(q1) - 2*c1*m1*r*cos(q1 - q2) - 2*m1*r*w1*sin(q1 - q2);
RHSh1 = I2*v1 + I1*(v1 - v2) + m1*(((c1*cos(r1 - r2) + w1*sin(r1 - r2))*(v1 - v2) - r*v1)*(c1*cos(r1 - r2) - r + w1*sin(r1 - r2)) + (c1*sin(r1 - r2) - w1*cos(r1 - r2))^2*(v1 - v2)) + m2*(v1*(c2*cos(r1) - r + w2*sin(r1))^2 + v1*(c2*sin(r1) - w2*cos(r1))^2);

u1 = Mh11 \ RHSh1; 

%===================================================================
function [q1,u1,q2,u2] = rear_strike(r1,v1,r2,v2,robot) %PAB (DONE)
%===================================================================

m1 = robot.m1; c1 = robot.c1; w1 = robot.w1;
m2 = robot.m2; c2 = robot.c2; w2 = robot.w2;
r = robot.r; 
I1 = robot.I1; I2 = robot.I2;

q1 = r1;
q2 = r2;
u2 = 0;
  
Mh11 = I1 + I2 + c1^2*m1 + c2^2*m2 + m1*r^2 + m2*r^2 + m1*w1^2 + m2*w2^2 - 2*c1*m1*r*cos(q1) - 2*m1*r*w1*sin(q1) - 2*c2*m2*r*cos(q1 - q2) - 2*m2*r*w2*sin(q1 - q2);
RHSh1 = I1*v1 + I2*(v1 - v2) + m2*(((c2*cos(r1 - r2) + w2*sin(r1 - r2))*(v1 - v2) - r*v1)*(c2*cos(r1 - r2) - r + w2*sin(r1 - r2)) + (c2*sin(r1 - r2) - w2*cos(r1 - r2))^2*(v1 - v2)) + m1*(v1*(c1*cos(r1) - r + w1*sin(r1))^2 + v1*(c1*sin(r1) - w1*cos(r1))^2);

u1 = Mh11 \ RHSh1;    
