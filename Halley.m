%% Halley comet orbit simulation

clc;
clear;
close all;

%adding an additional comment

addpath("Additional Functions\Ephemerides\");
addpath("Additional Functions\Time converstion\");

miu_sun = 1.327124e11; % km3/s2
miu_earth = 398600; % km3/s2
R_sun= 696340; % km
R_earth = 6378.14; % km
AU=1.496e8; %km

scale=5; % For animation purposes, scales the tail of the comet

%% HalLey Orbital Elements
a=2667950017;
e=0.96714291;
i=deg2rad(162.26269058);
raan=deg2rad(58.42);
aop=deg2rad(111.33248510452);
ta0=-6*pi/8;
date_0=[1986 2 5 21 29 15];
date_0_mjd2000=date2mjd2000(date_0);

pars0=[a,e,i,raan,aop,ta0];

y0=op2vec(pars0,miu_sun);

%% Additional Calculations
orbital_period= 2*pi*sqrt(a^3/miu_sun);
orbital_period_y=orbital_period/(86400*365);
distance_range=4*AU;
ta_range=acos((a*(1-e^2)/(distance_range)-1)/e);
t_range=ta2t(ta_range,pars0,miu_sun);

%% Closest approach
%{
min_d = 4*AU;
for rev=1:100
    date_perihel_mjd2000=date_0_mjd2000+rev*orbital_period_y;
    time=linspace(0,2*t_range,1000);
    ta_temp_0=-ta_range;
    pars0_temp=[a,e,i,raan,aop,ta_temp_0];
    y0_temp=op2vec(pars0_temp,miu_sun);
    options=odeset('RelTol',1e-12,'AbsTol',1e-12);
    [~,y_temp]=ode45(@OrbitFunc,time,y0_temp,options,miu_sun);
    
end
%}
%% Integration

orbital_period= 2*pi*sqrt(a^3/miu_sun);
disp(['Orbital period: ' num2str(orbital_period/(86400*365)) ' Years']);

time=linspace(0,orbital_period,10000);
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
[~,y]=ode45(@OrbitFunc,time,y0,options,miu_sun);


% Earth's orbit
a_e=AU;
e_e=0.0167;
i_e=0;
raan_e=0;
aop_e=0;
ta_e=0;
pars_e=[a_e,e_e,i_e,raan_e,aop_e,ta_e];
y0_e=op2vec(pars_e,miu_sun);
orbital_period_e= 2*pi*sqrt(a_e^3/miu_sun);
time_e=linspace(0,orbital_period_e,1000);
[~,y_e]=ode45(@OrbitFunc,time_e,y0_e,options,miu_sun);
%% Plots

figure();
title('Orbital view');
plot3(y(:,1),y(:,2),y(:,3));
hold on;
drawSphere(R_sun);
plot3(y_e(:,1),y_e(:,2),y_e(:,3),'r');
axis equal;
grid on;

%% Animation

figure()
plot3(y(:,1),y(:,2),y(:,3));
hold on;
drawSphere(10*R_sun);
plot3(y_e(:,1),y_e(:,2),y_e(:,3),'r');
axis equal;
grid on;
xlim(5*AU*[-1 1]);
ylim(5*AU*[-1 1]);
zlim(5*AU*[-1 1]);
for k=1:length(time)
    comet=plot3(y(k,1),y(k,2),y(k,3),'or');
    tail=drawTail(y(k,1:3),scale);
    drawnow;
    pause(1e-1);
    delete(comet);
    delete(tail);
    
    % Break statement
    op_temp = vec2op (y(k,:),miu_sun);
    if and(op_temp(6)>3*pi/4,op_temp(6)<(2*pi-3*pi/4)) break; end
end
hold off;

%% Functions
function p = drawSphere (r,o,n)
% r is radius
% o is the coordination of the center point
% n is the number of faces (default is 20 by 20)

switch nargin
    case 1
        [x,y,z]=sphere;
        p=surf(x*r,y*r,z*r);
    case 2
        [x,y,z]=sphere;
        p=surf(x*r+o(1),y*r+o(2),z*r+o(3));
    case 3
        [x,y,z]=sphere(n);
        p=surf(x*r+o(1),y*r+o(2),z*r+o(3));
end
end
function p = drawTail (r,scale)
if nargin==1
    scale=1;
end
AU=1.496e8; %km
MaxTail=50e6*scale; %km
Halley_perhel=0.586*AU;
if norm(r)<(MaxTail+Halley_perhel)
    tip=r;
    tail= r + (MaxTail+Halley_perhel-norm(r))*r/norm(r);
    p = line([tip(1) tail(1)],[tip(2) tail(2)],[tip(3) tail(3)],'Color','red');
else
    p = plot3(NaN,NaN,NaN);
end
end
function T = TransMat (theta,axis)
switch axis
    case 1
        T = [1 0 0 ; 0 cos(theta) sin(theta) ; 0 -sin(theta) cos(theta)];
    case 2
        T = [cos(theta) 0 -sin(theta) ; 0 1 0 ; sin(theta) 0 cos(theta)];
    case 3
        T = [cos(theta) sin(theta) 0 ; -sin(theta) cos(theta) 0 ; 0 0 1];
end
end
function ta = t2ta (t,pars,miu)
a = pars(1);
e = pars(2);
func = @(E) E-e*sin(E)-sqrt(miu/a^3)*t;
E = fzero(func,pi);
ta = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
end
function t = ta2t (ta,pars,miu)
a = pars(1);
e = pars(2);
E = 2*atan( tan(ta/2)*sqrt((1-e)/(1+e)) );
t = (E-e*sin(E))*sqrt(a^3/miu);
end
function vecs = op2vec (pars,miu)
a=pars(1);
e=pars(2);
i=pars(3);
raan=pars(4);
aop=pars(5);
ta=pars(6);

h       = sqrt(miu*a*(1-e^2));
r_size  = a*(1-e^2)/(1+e*cos(ta));
vr_size = miu/h*e*sin(ta);
vs_size = miu/h*(1+e*cos(ta));

rsw2eci=TransMat(-raan,3)*TransMat(-i,1)*TransMat(-(ta+aop),3);

r  = rsw2eci*[1;0;0]*r_size;
vr = rsw2eci*[1;0;0]*vr_size;
vs = rsw2eci*[0;1;0]*vs_size;
v  = vr+vs;

vecs = [r;v];
end
function pars = vec2op (vecs,miu)
% Position and velocity vector to orbital parameter converter
r=vecs(1:3);
v=vecs(4:6);

a = -miu/2*(1/(norm(v)^2/2-miu/norm(r)));

h = cross(r,v);

e = cross(v,h)/miu-r/norm(r);

i = real(acos(dot(h,[0;0;1])/norm(h)));

an = cross([0;0;1],h)/norm(cross([0;0;1],h));

raan = real(acos(dot(an,[1;0;0])));
if an(2)<0 
    raan=2*pi-raan; 
end

aop = real(acos(dot(e,an)/norm(e)));
if e(3)<0 
    aop=2*pi-aop; 
end

ta = real(acos(dot(r,e)/(norm(e)*norm(r))));
if dot(v,r)<=0
    ta=2*pi-ta; 
end

pars = [a,norm(e),i,raan,aop,ta];
end
function [dy,parout] = OrbitFunc (t,y,miu)
if ~iscolumn(y)
    y=y';
end
r_vec=y(1:3);
v_vec=y(4:6);
r_size=norm(r_vec);
v_size=norm(v_vec);
dv_vec=-miu/r_size^3*r_vec;
dy = [v_vec;dv_vec];

orbital_pars = vec2op (y,miu);
parout=[r_size,v_size,orbital_pars];

end