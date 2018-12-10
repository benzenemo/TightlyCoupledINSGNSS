function [out_profile,out_IMU_bias_est,out_clock,out_KF_SD,PickSubsetRes,RecClockBiasRes,out_MeasurementNoise_SD,out_Resi,InnovationRes,StdInnovationRes] =...
    Tightly_coupled_INS_GNSS(FilePath,old_time,old_est_r_eb_e,old_est_v_eb_e,...
    est_clock,attitude_ini,GNSS_config,TC_KF_config,L_ba_b,Total_GNSS_epoch)
%Tightly_coupled_INS_GNSS - Simulates inertial navigation using ECEF
% navigation equations and kinematic model, GNSS and tightly coupled
% INS/GNSS integration. 
%
% Software for use with "Principles of GNSS, Inertial, and Multisensor
% Integrated Navigation Systems," Second Edition.
%
% This function created 12/4/2012 by Paul Groves
% Inputs:
%   FilePath.GNSSFile   GNSS Obs File For Tightly_coupled_INS_GNSS
%       % Format of GNSS observation array:
%           Column 1: epoch
%           Column 2: Obsweek (GPS week)
%           Column 3: Obssec (GPS second)
%           Column 4: PRN (PRN=GPS_PRN=GLONASS_PRN+32=BDS_PRN+56)
%           Column 5: Ionospheric Free pseudorange linear combination (m)
%           Column 6: slant tropospheric delay (m)
%           Column 7: Satellite clock error (code bias have been correted) (s)
%           Column 8: relativity correction (m)
%           Column 9-11: Satellite position in ECEF(Earth rotation correction have been made,ECEF at Receiving moment)(m)
%           Column 12: range rate(m/s)
%           Column 13: Rate of Satellite clock  (s/s)
%           Column 14-16: Satellite velocity in ECEF(Earth rotation correction have been made,ECEF at Receiving moment)(m/s)
%           Column 17: flag= 0 means this PRN was Removed in  GNSS Single point Least squares Solution
%           Column 18: Elevation angle (deg)
%           Column 19: Azimuth (deg)
%           Column 20: user ranging error (By RTK solution) (m)
%           Column 21: GNSS Single point Least squares Solution Residual (m)
%           Column 22: a priori Pseudo-distance noise standard deviation (m)
%           Column 23-24:signal strength of each frequency(dbHz)
%           Column 25：Ionosphere delay by KLOBUCHAR(if Available)(m)
%   FilePath.INSFile    INS Obs File For Tightly_coupled_INS_GNSS
%       % Format of INS observation array:
%           Column 1: Obssec (GPS second)
%           Column 2: gyroscope observation (Forward) (rad/s)
%           Column 3: gyroscope observation (Right) (rad/s)
%           Column 4: gyroscope observation (Down) (rad/s)
%           Column 5: accelerometer observation (Forward) (m*s-2)
%           Column 6: accelerometer observation (Right) (m*s-2)
%           Column 7: accelerometer observation (Down) (m*s-2)
%   old_time,old_est_r_eb_e,old_est_v_eb_e,est_clock,attitude_ini 
%       Initial alignment result
%   GNSS_config 
%       defined in INS_GNSS_Demo_7
%   TC_KF_config
%       defined in INS_GNSS_Demo_7
%   L_ba_b
%       The vector from the IMU center to the phase center of the GNSS antenna (forward, right, down)
%   Total_GNSS_epoch
%       Number of GNSS observations
% Outputs:
%   out_profile        Navigation solution as a motion profile array
%        Format of output IMU biases array:
%          Column 1: time (sec)
%          Column 2: estimated ECEF-X position (m)
%          Column 3: estimated ECEF-Y position (m)
%          Column 4: estimated ECEF-Z position (m)
%          Column 5: estimated ECEF-X velocity (m/s)
%          Column 6: estimated ECEF-Y velocity (m/s)
%          Column 7: estimated ECEF-Z velocity (m/s)
%          Column 8: estimated roll angle (rad)
%          Column 9: estimated pitch angle (rad)
%          Column 10: estimated yaw angle (rad)
%          Column 11: estimated Local geodetic frame (origin defined by STA ) North position (m)
%          Column 12: estimated Local geodetic frame (origin defined by STA ) East position (m)
%          Column 13: estimated Local geodetic frame (origin defined by STA ) Up position (m)
%   out_IMU_bias_est   Kalman filter IMU bias estimate array
        % Format of output IMU biases array:
        %  Column 1: time (sec)
        %  Column 2: estimated X accelerometer bias (m/s^2)
        %  Column 3: estimated Y accelerometer bias (m/s^2)
        %  Column 4: estimated Z accelerometer bias (m/s^2)
        %  Column 5: estimated X gyro bias (rad/s)
        %  Column 6: estimated Y gyro bias (rad/s)
        %  Column 7: estimated Z gyro bias (rad/s)
%   out_clock          GNSS Receiver clock estimate array
        % Format of receiver clock array:
        %  Column 1: time (sec)
        %  Column 2: estimated clock offset (m)
        %  Column 3: estimated clock drift (m/s)
%   out_KF_SD          Output Kalman filter state uncertainties
        % Format of KF state uncertainties array:
        %  Column 1: time (sec)
        %  Column 2: X attitude error uncertainty (rad)
        %  Column 3: Y attitude error uncertainty (rad)
        %  Column 4: Z attitude error uncertainty (rad)
        %  Column 5: X velocity error uncertainty (m/s)
        %  Column 6: Y velocity error uncertainty (m/s)
        %  Column 7: Z velocity error uncertainty (m/s)
        %  Column 8: X position error uncertainty (m)
        %  Column 9: Y position error uncertainty (m)
        %  Column 10: Z position error uncertainty (m)
        %  Column 11: X accelerometer bias uncertainty (m/s^2)
        %  Column 12: Y accelerometer bias uncertainty (m/s^2)
        %  Column 13: Z accelerometer bias uncertainty (m/s^2)
        %  Column 14: X gyro bias uncertainty (rad/s)
        %  Column 15: Y gyro bias uncertainty (rad/s)
        %  Column 16: Z gyro bias uncertainty (rad/s)
        %  Column 17: clock offset uncertainty (m)
        %  Column 18: clock drift uncertainty (m/s)
%   PickSubsetRes      PDOP 、Num of Sat、time cost of Satellite Selection  
%   RecClockBiasRes    存储各卫星反算的接收机钟差同RecClockPre的差异
%   InnovationRes      新息矩阵， RecClockPreprocOptions==2时与RecClockBiasRes相同
%   out_MeasurementNoise_SD     measurement noise SD output record(m)
%   out_Resi            Residual record 
%   StdInnovationRes  标准化后的新息
% Copyright 2012, Paul Groves
% License: BSD; see license.txt for details
% Modified 2017/8 by LiuXiao BUAA benzenemo@buaa.edu.cn % 20170311B104ZXY
% Begins
global STA
load(FilePath.GNSSFile);
load(FilePath.INSFile);
if TC_KF_config.DoubleGNSSFlag==1
    load(FilePath.DoubleGNSSHeadingFile);
end
[no_epochs,~]=size(IMU);
load('..\SharedMat\MOD_constant')

% Initialize true navigation solution

% Determine satellite positions and velocities

% Initialize the GNSS biases. Note that these are assumed constant throughout 
% the simulation and are based on the initial elevation angles. Therefore, 
% this function is unsuited to simulations longer than about 30 min.

% Generate GNSS measurements

% Determine Least-squares GNSS position solution
[old_est_L_b,old_est_lambda_b,old_est_h_b,old_est_v_eb_n] =...
    pv_ECEF_to_NED(old_est_r_eb_e,old_est_v_eb_e);
est_L_b = old_est_L_b;
est_lambda_b=old_est_lambda_b;
% Initialize estimated attitude solution
old_est_C_b_n=Euler_to_CTM(attitude_ini)';%转置符号是我加的，函数Euler_to_CTM中的欧拉角并非我们定义的姿态角
[~,~,old_est_C_b_e] = NED_to_ECEF(old_est_L_b,...
    old_est_lambda_b,old_est_h_b,old_est_v_eb_n,old_est_C_b_n);

% Initialize output profile record and errors record
out_profile = zeros(no_epochs,13);
out_IMU_bias_est=zeros(Total_GNSS_epoch,7);
out_clock=zeros(Total_GNSS_epoch,4);
out_KF_SD=zeros(Total_GNSS_epoch,18);
PickSubsetRes=zeros(Total_GNSS_epoch,4);
[GNSSObsNum,~]=size(GNSSObs);
RecClockBiasRes=zeros(GNSSObsNum,5);
if strcmp(TC_KF_config.KFMethod,'M-LSKF')
    MLS.MAXloop=10;%迭代次数大于该值，则停止迭代
    MLS.DeltaPosi=0.01;%迭代时三维位置变化小于该值，则停止迭代
    MLS.OutlierCount=0;%残差超限的观测个数
%     References:杨元喜. 自适应动态导航定位[M]. 测绘出版社, 2017. P100
    MLS.WeightFunctionCategory=3;%2表示使用不归零的2段权函数 3表示使用3段权函数
    MLS.PseudorangeC=2.5;%伪距不归零2段权函数标准化残差分段点
    MLS.PseudorangeRateC=2.5;%伪距率不归零2段权函数标准化残差分段点
    MLS.K0=2.0;%三段权函数标准化残差分段点K0 
    MLS.K1=5.0;%三段权函数标准化残差分段点K1
elseif strcmp(TC_KF_config.KFMethod,'IAE-KF')
    IAE.WeightFunctionCategory=3;%2表示使用不归零的2段权函数 3表示使用3段权函数
    IAE.K0=1.5;%三段权函数标准化新息分段点K0 
    IAE.K1=5.0;%三段权函数标准化新息分段点K1 
end


InnovationRes=zeros(GNSSObsNum,5);
StdInnovationRes=zeros(GNSSObsNum,5);

out_MeasurementNoise_SD=zeros(GNSSObsNum,4);
out_Resi=zeros(GNSSObsNum,7);
% Generate output profile record
out_profile(1,1) = old_time;
out_profile(1,2:4)=(old_est_r_eb_e+old_est_C_b_e*L_ba_b)';
out_profile(1,5:7) = old_est_v_eb_e';
out_profile(1,8:10) = CTM_to_Euler(old_est_C_b_n')';
Rotation=[   -sin(est_L_b)*cos(est_lambda_b),-sin(est_L_b)*sin(est_lambda_b),cos(est_L_b);
                    -sin(est_lambda_b),cos(est_lambda_b),0;
                    cos(est_L_b)*cos(est_lambda_b),cos(est_L_b)*sin(est_lambda_b),sin(est_L_b)];
NEU=Rotation*(out_profile(1,2:4)-STA.STA(1).Coor(1:3))';
out_profile(1,11:13)=NEU';

% Determine errors and generate output record
% Initialize Kalman filter P matrix and IMU bias states
P_matrix = Initialize_TC_P_matrix(TC_KF_config);
est_IMU_bias = zeros(6,1);


% Generate IMU bias and clock output records
out_IMU_bias_est(1,1) = old_time;
out_IMU_bias_est(1,2:7) = est_IMU_bias';
out_clock(1,1) = old_time;
out_clock(1,2:3) = est_clock;
out_clock(1,4)=0;
out_clock(1,5)=0;
% Generate KF uncertainty record
out_KF_SD(1,1) = old_time;
for i =1:17
    out_KF_SD(1,i+1) = sqrt(P_matrix(i,i));
end % for i
% Initialize GNSS model timing
time_last_GNSS = floor(old_time);%GNSS观测更新时刻
time_last_ZVT = floor(old_time);%静止检测时刻
TC_KF_config.StationarityFlag = 0;%静止标识，默认为0标识非静止
GNSS_epoch = 1;
% Progress bar
dots = '....................';
bars = '||||||||||||||||||||';
rewind = '\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b';
fprintf(strcat('Processing: ',dots));
progress_mark = 0;
progress_epoch = 0;
% Main loop
for epoch = 2:no_epochs

    % Update progress bar
    if (epoch - progress_epoch) > (no_epochs/20)
        progress_mark = progress_mark + 1;
        progress_epoch = epoch;
        fprintf(strcat(rewind,bars(1:progress_mark),...
            dots(1:(20 - progress_mark))));
    end % if epoch    
    
    % Input data from motion profile
    time=IMU(epoch,1);

    % Time interval
    tor_i = time - old_time;
    
    % Calculate specific force and angular rate
    meas_f_ib_b=IMU(epoch,5:7)';
    meas_omega_ib_b=IMU(epoch,2:4)';
    % Correct IMU errors
    meas_f_ib_b = meas_f_ib_b - est_IMU_bias(1:3);
    meas_omega_ib_b = meas_omega_ib_b - est_IMU_bias(4:6);
    
    if sum(TC_KF_config.StationarityDetectionFlag)~=0&&time - time_last_ZVT >= TC_KF_config.StationarityDetectionInterval
        % 静止检测
        TC_KF_config = StationarityDetection(TC_KF_config,time);
        time_last_ZVT=time;
    end
    % Update estimated navigation solution
    [est_r_eb_e,est_v_eb_e,est_C_b_e] = Nav_equations_ECEF(tor_i,...
        old_est_r_eb_e,old_est_v_eb_e,old_est_C_b_e,meas_f_ib_b,...
        meas_omega_ib_b);

    % Determine whether to update GNSS simulation and run Kalman
    % filter%此处改为是否更新GNSS观测
    if (time - time_last_GNSS) >= GNSS_config.epoch_interval
        %依据观测时刻
        index_GNSS=find(GNSSObs(:,3)==floor(time)&~ismember(GNSSObs(:,4),GNSS_config.omit));
        no_GNSS_meas=length(index_GNSS);
        if no_GNSS_meas < 1
            disp(['Time' num2str(time) 'GNSS观测缺失'])
        else
            GNSS_epoch = GNSS_epoch + 1;
            tor_s = floor(time) - time_last_GNSS;  % KF time interval
            time_last_GNSS = floor(time);
            est_r_ea_e_for_RecClockPre=est_r_eb_e+est_C_b_e*L_ba_b; %将INS位置换算为GNSS天线位置
            est_v_ea_e_for_RecClockPre=est_v_eb_e;
            %扣除ISB
            index_BDS=ismember(GNSSObs(index_GNSS,4),GNSS_config.BDSPRNList);
            if ~isempty(index_BDS)
                GNSSObs(index_GNSS(index_BDS),5)=GNSSObs(index_GNSS(index_BDS),5)-GNSS_config.ISBBDS_GPS;
            end
            %筛选异常卫星(高度截止角、SNR、新息筛选卫星)
            [index_elimi,RecClockPre,RecClockBias,Innovation,RecClockRateBias]=RemoveOutlierGNSSObs(est_r_ea_e_for_RecClockPre,...
                est_v_ea_e_for_RecClockPre,GNSSObs(index_GNSS,:),est_clock,tor_s,TC_KF_config.RecClockPreprocOptions,GNSS_config);
                RecClockBiasRes(index_GNSS,:)=[GNSSObs(index_GNSS,2:4),RecClockBias,RecClockRateBias];
                if TC_KF_config.RecClockPreprocOptions==0
                    InnovationRes(index_GNSS,:)=[GNSSObs(index_GNSS,2:4),Innovation(1:no_GNSS_meas),Innovation(no_GNSS_meas+1:end)];
                end
            index_GNSS(index_elimi)=[];%剔除异常卫星
            no_GNSS_meas=length(index_GNSS);    %重新统计卫星数  
            if no_GNSS_meas < 1
                disp(['Time' num2str(time) 'GNSS观测在RemoveOutlierGNSSObs环节被全部剔除'])
            else
                %依据指定的卫星数目选星
                PDOP=999;
                if GNSS_config.intend_no_GNSS_meas<no_GNSS_meas            %判断可用卫星个数与待选卫星个数大小关系
                    %返回剔除卫星的Mask
                    [index_elimi,PDOP,RunningTime]=PickSubsetOfGNSS('RedundancyMatrix',GNSS_config.intend_no_GNSS_meas,est_r_ea_e_for_RecClockPre,GNSSObs(index_GNSS,9:11));
                    index_GNSS(index_elimi)=[];%完成剔星
                    no_GNSS_meas=length(index_GNSS);    %重新统计卫星数  
                elseif no_GNSS_meas>=4
                    RunningTime=0;
                    LOS=ones(no_GNSS_meas,1)*est_r_ea_e_for_RecClockPre'-GNSSObs(index_GNSS,9:11);
                    for i=1:no_GNSS_meas
                        LOS(i,1:3)=LOS(i,1:3)/norm(LOS(i,1:3));
                    end
                    PDOP=sqrt(trace(inv(LOS'*LOS)));
                end
                PickSubsetRes(GNSS_epoch,1)=time;
                PickSubsetRes(GNSS_epoch,2)=PDOP;
                PickSubsetRes(GNSS_epoch,3)=RunningTime;
                PickSubsetRes(GNSS_epoch,4)=no_GNSS_meas;
                % Generate GNSS measurements%此处改为读取伪距伪距率反推发射时刻求卫星位置，反推需要知道接收机的时间的GNSST，如何得到？对原始伪距做修正，要扣除电离层、对流层、卫星钟的影响？原始伪距率的预处理也需要考虑
                GNSS_measurements=[GNSSObs(index_GNSS,5)-GNSSObs(index_GNSS,6)+GNSSObs(index_GNSS,7)*c+GNSSObs(index_GNSS,8),...
                    GNSSObs(index_GNSS,12)+GNSSObs(index_GNSS,13)*c,GNSSObs(index_GNSS,9:11),GNSSObs(index_GNSS,14:16),GNSSObs(index_GNSS,18)];   
                if TC_KF_config.RecClockPreprocOptions==2%强制全程使用钟差预处理策略
                    Clock_Reset_Flag=1;%1表示GNSS数据做了钟跳修复，状态更新时钟差钟漂置零
                    GNSS_measurements(:,1)=GNSS_measurements(:,1)-RecClockPre;
                else
                    RecClockPre=0;
                    Clock_Reset_Flag=0;%1表示GNSS数据做了钟跳修复，状态更新时钟差钟漂置零,默认为0
                end
                %读取双天线测向数据
                TC_KF_config.DoubleGNSSDataFlag=0;%默认无双天线数据
                DoubleGNSSHeadingEpoch=zeros(1,2);%双天线测量值默认为0
                if TC_KF_config.DoubleGNSSFlag==1
                    index_DoubleGNSS=find(DoubleGNSSHeading.sow==floor(time)&DoubleGNSSHeading.HdgSD<20);
                    if ~isempty(index_DoubleGNSS)
                        TC_KF_config.DoubleGNSSDataFlag=1;
                        DoubleGNSSHeadingEpoch=[DoubleGNSSHeading.Heading(index_DoubleGNSS),DoubleGNSSHeading.HdgSD(index_DoubleGNSS)]*pi/180;
                    end
                end
                if strcmp(TC_KF_config.KFMethod,'ClassicKF')
                    % Run Integration Kalman filter
                    [est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,R_matrix] =...
                        TC_KF_Epoch(GNSS_measurements,no_GNSS_meas,tor_s,est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,...
                        meas_f_ib_b,TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,DoubleGNSSHeadingEpoch);
                elseif strcmp(TC_KF_config.KFMethod,'M-LSKF')
                    MLS_loopflag=true;
                    MLS_loopCount=0;
                    MLS.OutlierCount=0;
                    index_LowEle=GNSS_measurements(:,9)<30;
                    ElevationGain=ones(no_GNSS_meas,1);
                    ElevationGain(index_LowEle)=1./(2*sind(GNSS_measurements(index_LowEle,9)));
                    if TC_KF_config.StationarityFlag == 0 && TC_KF_config.DoubleGNSSDataFlag == 0%静止检测结果为动态且无双天线观测
                        R_matrix=eye(no_GNSS_meas*2);
                    elseif TC_KF_config.StationarityFlag == 1 && TC_KF_config.DoubleGNSSDataFlag == 0%静止检测结果为静态且无双天线观测
                        R_matrix=eye(no_GNSS_meas*2+1);
                    elseif TC_KF_config.StationarityFlag == 0 && TC_KF_config.DoubleGNSSDataFlag == 1%静止检测结果为动态且有双天线观测
                        R_matrix=eye(no_GNSS_meas*2+3);
                    elseif TC_KF_config.StationarityFlag == 1 && TC_KF_config.DoubleGNSSDataFlag == 1%静止检测结果为静态且有双天线观测
                        R_matrix=eye(no_GNSS_meas*2+4);
                    end

                    R_matrix(1:no_GNSS_meas,1:no_GNSS_meas) = diag(ElevationGain)*...
                    TC_KF_config.pseudo_range_SD^2;

                    R_matrix((no_GNSS_meas + 1):2*no_GNSS_meas,(no_GNSS_meas + 1):2*no_GNSS_meas) =...
                        diag(ElevationGain.*ElevationGain)...
                        * TC_KF_config.range_rate_SD^2;

                    if  TC_KF_config.StationarityFlag == 1 && TC_KF_config.DoubleGNSSDataFlag == 0%静止检测结果为静态且无双天线观测
                        R_matrix(2*no_GNSS_meas+1,2*no_GNSS_meas+1)=TC_KF_config.ZARU_DGrySD^2;
                    elseif TC_KF_config.StationarityFlag == 0 && TC_KF_config.DoubleGNSSDataFlag == 1%静止检测结果为动态且有双天线观测
                        R_matrix(2*no_GNSS_meas+1:2*no_GNSS_meas+3,2*no_GNSS_meas+1:2*no_GNSS_meas+3)=eye(3)*DoubleGNSSHeadingEpoch(1,2)^2;
                    elseif TC_KF_config.StationarityFlag == 1 && TC_KF_config.DoubleGNSSDataFlag == 1%静止检测结果为静态且有双天线观测 
                        R_matrix(2*no_GNSS_meas+1,2*no_GNSS_meas+1)=TC_KF_config.ZARU_DGrySD^2;
                        R_matrix(2*no_GNSS_meas+2:2*no_GNSS_meas+4,2*no_GNSS_meas+2:2*no_GNSS_meas+4)=eye(3)*DoubleGNSSHeadingEpoch(1,2)^2;
                    end
                    while MLS_loopflag
                        if MLS_loopCount==0
                            est_r_eb_e_temp_old=est_r_eb_e;
                        else
                            R_matrix=R_matrix_NextLoop;
                        end
                        [est_C_b_e_temp,est_v_eb_e_temp,est_r_eb_e_temp,est_IMU_bias_temp,est_clock_temp,P_matrix_temp,R_matrix_NextLoop,MLS] =...
                            MLS_KF_Epoch(GNSS_measurements,no_GNSS_meas,tor_s,est_C_b_e,...
                            est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,...
                            meas_f_ib_b,TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,R_matrix,MLS,DoubleGNSSHeadingEpoch);
                        MLS_loopCount=MLS_loopCount+1;
                        if MLS_loopCount==1
                            Resi=MLS.Resi;
                            StdResi=MLS.StdResi;
                        end                        
                        if MLS.OutlierCount==0 %残差全部在门限之内 不进行迭代重新赋权
                            break
                        end
                        if norm(est_r_eb_e_temp_old-est_r_eb_e_temp,2)<MLS.DeltaPosi||MLS_loopCount>=MLS.MAXloop
                            MLS_loopflag=false;
                        end
                        est_r_eb_e_temp_old=est_r_eb_e_temp;
                    end
                    est_C_b_e=est_C_b_e_temp;
                    est_v_eb_e=est_v_eb_e_temp;
                    est_r_eb_e=est_r_eb_e_temp;
                    est_IMU_bias=est_IMU_bias_temp;
                    est_clock=est_clock_temp;
                    P_matrix=P_matrix_temp;
                elseif strcmp(TC_KF_config.KFMethod,'IAE-KF')
                    index_LowEle=GNSS_measurements(:,9)<30;
                    ElevationGain=ones(no_GNSS_meas,1);
                    ElevationGain(index_LowEle)=1./(2*sind(GNSS_measurements(index_LowEle,9)));
                    R_matrix=eye(no_GNSS_meas*2);
                    R_matrix(1:no_GNSS_meas,1:no_GNSS_meas) = diag(ElevationGain)*...
                    TC_KF_config.pseudo_range_SD^2;
                    R_matrix((no_GNSS_meas + 1):end,(no_GNSS_meas + 1):end) =...
                        diag(ElevationGain)...
                        * TC_KF_config.range_rate_SD^2;
                    [est_C_b_e,est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,IAE,R_matrix] =...
                        IAE_KF_Epoch(GNSS_measurements,no_GNSS_meas,tor_s,est_C_b_e,...
                        est_v_eb_e,est_r_eb_e,est_IMU_bias,est_clock,P_matrix,...
                        meas_f_ib_b,est_L_b,TC_KF_config,L_ba_b,meas_omega_ib_b,Clock_Reset_Flag,R_matrix,IAE);
                end
                % Generate IMU bias and clock output records
                out_IMU_bias_est(GNSS_epoch,1) = time;
                out_IMU_bias_est(GNSS_epoch,2:7) = est_IMU_bias';
                out_clock(GNSS_epoch,1) = time;
                out_clock(GNSS_epoch,2:3) = est_clock;
                out_clock(GNSS_epoch,4)=RecClockPre;
                % Generate KF uncertainty output record
                out_KF_SD(GNSS_epoch,1) = time;
                for i =1:17
                    out_KF_SD(GNSS_epoch,i+1) = sqrt(P_matrix(i,i));
                end % for i
                %measurement noise SD output record(m)
                if ~strcmp(TC_KF_config.KFMethod,'ClassicKF')
                    out_MeasurementNoise_SD(index_GNSS,:)=[GNSSObs(index_GNSS,2:4),sqrt(diag(R_matrix(1:no_GNSS_meas,1:no_GNSS_meas)))];
                end
                if strcmp(TC_KF_config.KFMethod,'M-LSKF')
                    out_Resi(index_GNSS,:)=[GNSSObs(index_GNSS,2:4),Resi(1:no_GNSS_meas),StdResi(1:no_GNSS_meas),Resi(no_GNSS_meas+1:2*no_GNSS_meas),StdResi(no_GNSS_meas+1:2*no_GNSS_meas)];
                end
                if strcmp(TC_KF_config.KFMethod,'IAE-KF')
                    StdInnovationRes(index_GNSS,:)=[GNSSObs(index_GNSS,2:4),IAE.StdInno(1:no_GNSS_meas),IAE.StdInno(no_GNSS_meas+1:2*no_GNSS_meas)];
                end
            end
        end
    end % if time

    % Convert navigation solution to NED
    [est_L_b,est_lambda_b,~,~,est_C_b_n] =...
        ECEF_to_NED(est_r_eb_e,est_v_eb_e,est_C_b_e);
    
    % Generate output profile record
    out_profile(epoch,1) = time;
    out_profile(epoch,2:4)=(est_r_eb_e+est_C_b_e*L_ba_b)';%为了便于比对，最后输出的是GNSS天线的位置
    out_profile(epoch,5:7) = (est_v_eb_e+est_C_b_e*( Skew_symmetric(meas_omega_ib_b)*L_ba_b))';%同理，输出GNSS天线速度，忽略地球自传相关的杆臂效应
    out_profile(epoch,8:10) = CTM_to_Euler(est_C_b_n')';
    Rotation=[   -sin(est_L_b)*cos(est_lambda_b),-sin(est_L_b)*sin(est_lambda_b),cos(est_L_b);
                        -sin(est_lambda_b),cos(est_lambda_b),0;
                        cos(est_L_b)*cos(est_lambda_b),cos(est_L_b)*sin(est_lambda_b),sin(est_L_b)];
    NEU=Rotation*(out_profile(epoch,2:4)-STA.STA(1).Coor(1:3))';
    out_profile(epoch,11:13)=NEU';
    % Determine errors and generate output record%若没有真实位置速度姿态做参考，本区域代码作废

    
    % Reset old values
    old_time = time;
    old_est_r_eb_e = est_r_eb_e;
    old_est_v_eb_e = est_v_eb_e;
    old_est_C_b_e = est_C_b_e;

end %epoch

% Complete progress bar
fprintf(strcat(rewind,bars,'\n'));

% Ends