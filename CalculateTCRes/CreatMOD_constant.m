%生成常量MOD_constant.mat文件供各部分子函数调用
clc
clear

MaxPRN=96; % used in other modules
SatNum=0 ;  % Satellite Number
ParaNum=4 ;  % Coordinate, Tropsphere parameter, and clock
c=299792458.0;   % velocity of light (m/s)
% Elements of WGS84 ellipsolide
SemiAxis=6378137;  %major semi-axis(meter), GGSP standard
ecc2=0.0066943800229d0;  %square of eccentricity
omg=7.2921151467d-5; %  WGS 84: Earth rotation angular velocity (rad/sec)
GM=3.986005d14 ;       %  WGS 84: universal gravitational param (m^3/s^2)
BDSSatClass=[1,1,1,1,1,2,2,2,2,2,3,3,2,3,0,0,1];%BDS卫星类型1为GEO 2为IGSO 3为MEO
save('MOD_constant.mat')


