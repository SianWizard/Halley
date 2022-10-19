%% Analytic way to calculate the closest approach of Earth and Halley

clc;
clear;
close all;


addpath("Additional Functions\Ephemerides\");
addpath("Additional Functions\Time converstion\");
addpath("Additional Functions\Other\");

%% Predefined values
miu_sun = 1.327124e11; % km3/s2
AU=1.496e8; %km


%% Earth
a_earth=AU;
e_earth=0.0167;
i_earth=0;

tp_earth=date2mjd2000([2022 1 4 5 52 0]); % Latest time of Earth at perihelion

x_earth=[0 1 0]';
y_earth=[-1 0 0]';
z_earth=[0 0 1]';

perifocal_earth=[x_earth y_earth z_earth];
% Earth_year_secs=2*pi*sqrt(AU^3/miu_sun);
% Earth_year_days=Earth_year_secs/86400;
% 
% t=date2mjd2000(tp_earth)-Earth_year_days*1000+365.2648/2;
% f = t2f (t);



%% Halley
a_haley=2667950017;
e_haley=0.96714291;
i_haley=deg2rad(162.26269058);
raan_haley=deg2rad(58.42);
aop_haley=deg2rad(111.33248510452);
tp_haley=date2mjd2000([1986 2 5 21 29 15]);% Latest time of Haley at perihelion
f0_haley=0;


[r0_haley, v0_haley] = par2car([a_haley e_haley i_haley raan_haley aop_haley f0_haley], miu_sun);
x_haley=r0_haley/norm(r0_haley);
z_haley=cross(r0_haley,v0_haley)/norm(cross(r0_haley,v0_haley));
y_haley=cross(z_haley,r0_haley)/norm(cross(z_haley,r0_haley));
perifocal_haley=[x_haley y_haley z_haley];
%% Main calculation

initial_time=[1980 1 1 0 0 0];
final_time=[2080 1 1 0 0 0];
n=10000;

t_initial=date2mjd2000(initial_time);
t_final=date2mjd2000(final_time);

data1=[a_earth e_earth tp_earth];
data2=[a_haley e_haley tp_haley];
A1=perifocal_earth;
A2=perifocal_haley;

t_vect=linspace(t_initial,t_final,n);
distances=zeros(1,n);

for i=1:n
    distances(i)=MainFunc (t_vect(i),data1,data2,A1,A2);
end

[ClosestApproach,i_CA]=min(distances); % Find the minimum

%% Plotting
CA_AU=ClosestApproach/AU
Time_CA=mjd20002date(t_vect(i_CA))

% Transforming date numbers
DateNumber=zeros(1,n);
for j=1:n
    timeArray=mjd20002date(t_vect(j));
    DateNumber(j)=datenum(timeArray(1),timeArray(2),timeArray(3),timeArray(4),timeArray(5),timeArray(6));
end

figure()
hold on;
plot(DateNumber,distances/AU);
plot(DateNumber(i_CA),CA_AU,'or','LineWidth',2);
datetick('x',26);
xlabel('Date');
ylabel('Distance [AU]');
xlim([DateNumber(1) DateNumber(end)]);
grid on;
grid minor;
hold off;

%% Functions

function f = t2f (t,tp,a,e,miu)
% t [mjd2000 days]
% default values are for earth, orbiting around sun
switch nargin
    case 1
        tp=[2022 1 4 5 52 0];
        tp=date2mjd2000(tp);
        miu=1.327124e11;
        a=1.496e8;
        e=0.0167;
    case 2
        miu=1.327124e11;
        a=1.496e8;
        e=0.0167;
    case 3
        miu=1.327124e11;
        e=0.0167;
    case 4
        miu=1.327124e11;
    case 5

end

options=optimset('TolX',1e-16);
func=@(E) E-e*sin(E)-sqrt(miu/a^3)*(t-tp)*86400;
try
    E=fzero(func,pi,options);
catch
    E=fsolve(func,pi,options);
    disp('Fsolve used');
end
f=2*atan(sqrt((1+e)/(1-e))*tan(E/2));
end
function r_vec = f2r_perifocal (a,e,f)
r_norm=a*(1-e^2)/(1+e*cos(f));
r_x=r_norm*cos(f);
r_y=r_norm*sin(f);
r_vec=[r_x r_y 0]';
end
function A_1to2 = RotMat1to2 (A1,A2)
% rotating rf 1 to rf 2
A_1to2=A2/A1;
end
function dist = MainFunc (t,data1,data2,A1,A2,miu)
% time = [mjd2000]
% data = [a e tp]
if nargin<6
    miu=1.327124e11;
end
data=[data1;data2];
r_vec=zeros(3,2);
for i=1:2
    a=data(i,1);
    e=data(i,2);
    tp=data(i,3);
    f= t2f (t,tp,a,e,miu);
    r_vec(:,i) = f2r_perifocal (a,e,f);
end
r1_perifocal=r_vec(:,1);
r2_perifocal=r_vec(:,2);
A_1to2 = RotMat1to2 (A1,A2);
r1_newRF2=A_1to2*r1_perifocal;
dist=norm(r2_perifocal-r1_newRF2);
end