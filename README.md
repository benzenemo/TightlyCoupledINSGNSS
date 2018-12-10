# Tightly Coupled INS GNSS
本程序基于教材：
Paul D. Groves, 练军想, 唐康华, 潘献飞. GNSS与惯性及多传感器组合导航系统原理 : Principles of GNSS, inertial, and multisensor integrated navigation systems. 北京:国防工业出版社, 2015.
附赠的INS/GNSS紧组合仿真代码进行修改得到，能够实现基于伪距、伪距率、INS数据的紧组合解算。
文件夹CalculateTCRes为紧组合解算程序，INS_GNSS_Demo_7为主脚本
文件夹SharedMat为运行所需常量数组
文件夹TCdata为一组手推车实验数据，包含InertialExplorer软件INS/RTK模式输出的参考导航解（DGNSSRES文件夹），预处理后的GNSS观测（GNSSObsForCouple文件夹），预处理后的INS观测（IMU文件夹），双GNSS天线测向数据（SPANE1文件夹）
