function [ EliminationMask,PDOP, RunningTime] = PickSubsetOfGNSS( method,IntendedSatNum,RecPosi,SatPosi )
%本程序作为选星的第二个环节（第一个环节为剔除异常卫星），剔除部分卫星使得剩余卫星数满足要求
%输入：
% method            选星方法标识
%                       'RedundancyMatrix' 参考文献：Park C W, How J P. Method and apparatus for selecting optimal satelittes in global positioning system: U.S. Patent 6,727,850[P]. 2004-4-27.
%                       'Optimal' 遍历所有卫星组合，找最小PDOP
% IntendedSatNum    希望保留的卫星个数
% RecPosi           先验的GNSS接收机天线位置
% SatPosi           卫星位置ECEF，每行为一颗卫星的位置，依次为XYZ，单位：m
% 输出：
% EliminationMask   列向量，0表示保留，1表示剔除，卫星序列与SatPosi的序列对应
% PDOP              剔星后的PDOP值
% RunningTime       选星程序运行时间(不含PDOP计算)
tic %开始计时
[SatNum,~]=size(SatPosi);
EliminationMask=false(SatNum,1);
if strcmpi(method, 'RedundancyMatrix')
    LOS=SatPosi-ones(SatNum,1)*RecPosi';
    MatrixD=eye(SatNum);
    for SatIndexRow=2:SatNum
        for SatIndexCol=1:SatIndexRow-1
            MatrixD(SatIndexRow,SatIndexCol)=dot(LOS(SatIndexRow,1:3),LOS(SatIndexCol,1:3))/...
                norm(LOS(SatIndexRow,1:3))/norm(LOS(SatIndexCol,1:3));
        end
    end
    MatrixRedundancy=MatrixD.*MatrixD;
    MatrixRedundancy=MatrixRedundancy+tril(MatrixRedundancy,-1)';
    while sum(EliminationMask)<SatNum-IntendedSatNum
        Redundancy=sum(MatrixRedundancy);
        [~,EliminationFlag]=max(Redundancy);
        EliminationMask(EliminationFlag)=true;
        MatrixRedundancy(EliminationFlag,:)=0;
        MatrixRedundancy(:,EliminationFlag)=0;
    end
elseif strcmpi(method, 'Optimal')
    AllCase=nchoosek(1:SatNum,IntendedSatNum);
    PDOPAllCase=zeros(length(AllCase),1);
    for CaseIndex=1:length(AllCase)
        EliminationMask=true(SatNum,1);%reset
        EliminationMask(AllCase(CaseIndex,:))=false;
        LOS=ones(sum(~EliminationMask),1)*RecPosi'-SatPosi(~EliminationMask,1:3);
        for i=1:sum(~EliminationMask)
            LOS(i,1:3)=LOS(i,1:3)/norm(LOS(i,1:3));
        end
        PDOPAllCase(CaseIndex)=sqrt(trace(inv(LOS'*LOS)));
    end
    [PDOP,IndexMinPDOP]=min(PDOPAllCase);
    EliminationMask=true(SatNum,1);%reset
    EliminationMask(AllCase(IndexMinPDOP,:))=false;
    RunningTime=toc;
    return
else
    disp(['GNSS 选星算法:' method '未定义'])
end
RunningTime=toc;
LOS=ones(sum(~EliminationMask),1)*RecPosi'-SatPosi(~EliminationMask,1:3);
for i=1:sum(~EliminationMask)
    LOS(i,1:3)=LOS(i,1:3)/norm(LOS(i,1:3));
end
PDOP=sqrt(trace(inv(LOS'*LOS)));
return


