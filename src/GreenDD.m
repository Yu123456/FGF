function [xx,RG,IG,time,optional]=GreenDD(internal,r,z,h,flag)
% 本程序编写于2014年8月24日，通过计算频域格林函数微分方程形式计算频域格林函数值
% 直接调用函数 [xx,OG,ODG,OGr,ODGr,OGz,ODGz,Optional]=GreenD(internal,h,r,z,flag)

% r,z  场点位置
% internal 频率区间
% h 步长
% flag 标号，flag=1,计算 G,DG, flag=2 计算G,DG,Gr,DGr,Gz,DGz

% RG 返回 [G,DG,Gr,DGr,Gz,DGz]实部
% IG 返回 [G,DG,Gr,DGr,Gz,DGz]虚部
% time 返回计算时间

switch flag
    case 1
        tic;
        [xx,OG,ODG,~,~,~,~,optional]=GreenD(internal,h,r,z,flag);
        RG=[OG(:,1),ODG(:,1)];
        IG=[OG(:,2),ODG(:,2)];
        time=toc;
    case 2
        tic;
        [xx,OG,ODG,OGr,ODGr,OGz,ODGz,optional]=GreenD(internal,h,r,z,flag);
        RG=[OG(:,1),ODG(:,1),OGr(:,1),ODGr(:,1),OGz(:,1),ODGz(:,1)];
        IG=[OG(:,2),ODG(:,2),OGr(:,2),ODGr(:,2),OGz(:,2),ODGz(:,2)];
        time=toc;
end

end

