%CreatMOD_STA生成测站配置结构体STA
clc
clear
load('ControlNet');
STA.Num=1;%total station number
STA.FixNum=0;% number of fixed stations 
STA.STA(1).Name='BUAAStaticBeiDouPPPRes2017Day320';%Station Name
STA.STA(1).SKD='K';%Station coordinate status, fix or kinemetic
STA.STA(1).Coor=[-2171981.63614967,4386001.44447237,4076182.47907183];
[L_b,lambda_b,h_b,~,~] = ECEF_to_NED(STA.STA(1).Coor',zeros(3,1),eye(3));
STA.STA(1).BLH=[L_b*180/pi,lambda_b*180/pi,h_b];
% STA.STA(1).BLH=[34.207029018 , 117.135303736  ,  45.4034];%Station BLH  latitude(deg) longitude(deg)  height(m) 
STA.STA(1).NEU=[0,0,0];%天线相对于坐标参考点的北东天位置
% [STA.STA(1).Coor,~,~]=NED_to_ECEF(STA.STA(1).BLH(1)*pi/180,STA.STA(1).BLH(2)*pi/180,STA.STA(1).BLH(3),zeros(3,1),eye(3));
% STA.STA(1).Coor=STA.STA(1).Coor';%Station Coordination
STA.STA(1).Rotation=[   -sind(STA.STA(1).BLH(1))*cosd(STA.STA(1).BLH(2)),-sind(STA.STA(1).BLH(1))*sind(STA.STA(1).BLH(2)),cosd(STA.STA(1).BLH(1));
                        -sind(STA.STA(1).BLH(2)),cosd(STA.STA(1).BLH(2)),0;
                        cosd(STA.STA(1).BLH(1))*cosd(STA.STA(1).BLH(2)),cosd(STA.STA(1).BLH(1))*sind(STA.STA(1).BLH(2)),sind(STA.STA(1).BLH(1))];
STA.STA(1).Coor=STA.STA(1).Coor+STA.STA(1).NEU*STA.STA(1).Rotation;%将坐标转到天线参考点
[L_b,lambda_b,h_b,~,~] = ECEF_to_NED(STA.STA(1).Coor',zeros(3,1),eye(3));
STA.STA(1).BLH=[L_b*180/pi,lambda_b*180/pi,h_b];
MJD=YMDHMS2Mjd(int_year,1,int_doy, 12, 0, 0.d0);
if strcmp(cdattype,'GPT')==1
    [STA.STA(1).Trop.press ,STA.STA(1).Trop.temp, STA.STA(1).Trop.rhumity, STA.STA(1).Trop.undu]...
        = gpt( MJD,STA.STA(1).BLH(1)*pi/180,STA.STA(1).BLH(2)*pi/180,STA.STA(1).BLH(3) );
elseif strcmp(cdattype,'GPT2')==1
    disp('气象函数GPT2尚未完成，请使用配置GPT')
    pause
end
if strcmp(cztd,'SAAS')==1
    [STA.STA(1).Trop.ZHD,STA.STA(1).Trop.ZWD]=ZTD_SAAS(STA.STA(1).Trop.press ,...
        STA.STA(1).Trop.temp, STA.STA(1).Trop.rhumity,STA.STA(1).BLH(1)*pi/180,...
        STA.STA(1).BLH(3));
elseif strcmp(cztd,'EGNOS')==1
    disp('天顶模型EGNOS尚未完成，请使用配置SAAS')
    pause    
end
save('MOD_STA_BUAAStaticBeiDouPPPRes2017Day320','STA')



    



