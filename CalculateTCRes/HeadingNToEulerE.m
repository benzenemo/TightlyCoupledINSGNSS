function Euler_E = HeadingNToEulerE( heading_n,est_L_b,est_lambda )
%将当地水平坐标系（NED）下的航向角转换为相对于ECEF系三轴做旋转，用于将双GNSS天线测向新息（航向角修正）转化为对ECEF系下的系统变量的直接观测
%输入：
% heading_n N系下的航向角 rad
% est_L_b 纬度 rad
% est_lambda 经度 rad
%输出：
% Euler_E=[roll;pitch;yaw] 均为ECEF系 rad

Delta_C_n_b=Euler_to_CTM([0,0,heading_n]);
cos_lat = cos(est_L_b);
sin_lat = sin(est_L_b);
cos_long = cos(est_lambda);
sin_long = sin(est_lambda);
C_n_e=[-sin_lat * cos_long, -sin_lat * sin_long,  cos_lat;...
               -sin_long,            cos_long,        0;...
     -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat]';
Delta_C_e_b= C_n_e*Delta_C_n_b*C_n_e';
Euler_E=CTM_to_Euler(Delta_C_e_b);

end

