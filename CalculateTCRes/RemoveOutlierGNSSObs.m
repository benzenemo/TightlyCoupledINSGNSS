function [ index_elimi,RecClockPre,RecClockBias,Innovation,RecClockRateBias] = RemoveOutlierGNSSObs( RecPosi,RecVelo,GNSSObs,est_clock_previous,tor_s,RecClockPreprocOptions,GNSS_config )
%剔除仰角、载噪比不符合门限要求的观测
%在执行定位解算之前，通过先验的接收机位置反算接收机的钟差，每个伪距观测都能够得到一个钟差值
%以各星钟差的中位数作为预估的接收机钟差RecClockPre
%若某颗卫星解算得到的钟差与RecClockPre差异(该差异即使用RecClockPre作为钟差计算到的新息)大于RecClockBiasThreshold，剔除
%RecClockPreprocOptions=2时RecClockBias与Innovation相等，RecClockPreprocOptions=0时解算Innovation使用的接收机钟差是通过一步递推得到的
% Inputs:
% RecPosi           先验的GNSS接收机天线位置
% RecVelo           先验接收机速度
% GNSSObs           预处理后的GNSS观测数组
%  Column 1: epoch
%  Column 2: Obsweek (week)
%  Column 3: Obssec (s)
%  Column 4: PRN
%  Column 5: Ionospheric Free pseudorange linear combination (m)
%  Column 6: slant tropospheric delay (m)
%  Column 7: Satellite clock error (s)
%  Column 8: relativity correction (m)
%  Column 9-11: Satellite position in ECEF(m)
%  Column 12: range rate(m/s)
%  Column 13: Rate of Satellite clock  (s/s)
%  Column 14-16: Satellite velocity in ECEF(m/s)
%  Column 17: flag= 0 means this PRN was Removed in  GNSS Single point Navigation Solution
%  Column 18: Elevation angle (deg)
%  Column 19: Azimuth (deg)
%  Column 20: user ranging error  (m)
%  Column 21: residual (m)
%  Column 22: a priori Pseudo-distance noise standard deviation (m)
% est_clock_previous            上一次Kalman滤波得到的接收机钟差钟速
%  Column 1: 接收机钟差引起的测距误差 （m）
%  Column 2: 接收机钟速 （m/s）
% tor_s  Time interval（s）
% RecClockPreprocOptions        0表示全程使用传统模型 2表示全程使用钟差修预处理正策略 
% GNSS_config
%     .epoch_interval     Interval between GNSS epochs (s)
%     .init_est_r_ea_e    Initial estimated position (m; ECEF)
%     .mask_angle         Mask angle (deg)
%     .mask_SignalStrenth Mask Carrier to noise ratio(dbhz)
%     .rx_clock_offset    Receiver clock offset at time=0 (m)
%     .rx_clock_drift     Receiver clock drift at time=0 (m/s)
%     .intend_no_GNSS_meas       可处理的GNSS卫星数上限，观测数量大于该值触发剔星操作
% Outputs:
% RecClockPre 先验接收机钟差引起的测距误差（m）
% index_elimi 列向量，0表示保留，1表示剔除卫星，卫星序列与GNSSObs的序列对应
% RecClockBias 各卫星计算得到的接收机钟差与RecClockPre的差异
% RecClockRateBias 各卫星计算得到的接收机钟速与RecClockPre的差异
% Innovation 新息 列数为观测数目的2倍，上半部分为伪距新息，下半部分为伪距率新息
c=299792458.0;   % velocity of light (m/s)
RecClockBiasThreshold=30;
InnovationThreshold=30;
[no_GNSS_meas,~]=size(GNSSObs);
% 计算预估的伪距
RecClockPre=median(GNSSObs(:,5)-GNSSObs(:,6)+GNSSObs(:,7)*c+GNSSObs(:,8)-...
     sqrt(dot(GNSSObs(:,9:11)'-RecPosi*ones(1,no_GNSS_meas),GNSSObs(:,9:11)'-RecPosi*ones(1,no_GNSS_meas)))');
%  计算预估的伪距率 忽略Sagnac
u_as_e_T = zeros(no_GNSS_meas,3);
RecClockRateAll =zeros(no_GNSS_meas,1);
for i=1:no_GNSS_meas
    delta_r =  GNSSObs(i,9:11)' - RecPosi; 
    range = sqrt(delta_r' * delta_r);
    u_as_e_T(i,1:3) = delta_r' / range;
    RecClockRateAll(i,1)= GNSSObs(i,12)-u_as_e_T(i,1:3) * (GNSSObs(i,14:16)'- RecVelo);
end
RecClockRatePre=median(RecClockRateAll); 
RecClockBias=GNSSObs(:,5)-GNSSObs(:,6)+GNSSObs(:,7)*c+GNSSObs(:,8)-...
     sqrt(dot(GNSSObs(:,9:11)'-RecPosi*ones(1,no_GNSS_meas),GNSSObs(:,9:11)'-RecPosi*ones(1,no_GNSS_meas)))'-RecClockPre;
RecClockRateBias= RecClockRateAll-RecClockRatePre;
if RecClockPreprocOptions==0
    RecClockUpdate=est_clock_previous(1)+est_clock_previous(2)*tor_s;
    Innovation=[GNSSObs(:,5)-GNSSObs(:,6)+GNSSObs(:,7)*c+GNSSObs(:,8)-...
        sqrt(dot(GNSSObs(:,9:11)'-RecPosi*ones(1,no_GNSS_meas),GNSSObs(:,9:11)'-RecPosi*ones(1,no_GNSS_meas)))'-RecClockUpdate*ones(no_GNSS_meas,1);...
        RecClockRateAll-est_clock_previous(2)*ones(no_GNSS_meas,1)];
    index_elimi=abs (Innovation(1:no_GNSS_meas))>InnovationThreshold;
else
    Innovation=[RecClockBias;RecClockRateBias];
    index_elimi=abs (RecClockBias)>RecClockBiasThreshold;
end
index_lowele=GNSSObs(:,18)<GNSS_config.mask_angle;
index_lowSignalStrenth=GNSSObs(:,23)<GNSS_config.mask_SignalStrenth|GNSSObs(:,24)<GNSS_config.mask_SignalStrenth;
index_elimi=index_lowele|index_lowSignalStrenth|index_elimi;
end

