function [RG,IG,time]=GreenFI(r,z,omega,flag)
% 本程序编写于2014年8月24日，通过直接计算频域格林函数积分方程形式计算频域格林函数值
% 直接调用函数 RG=GreenF(r,z,omega,flag)

% r,z  场点位置
% omega 频率
% flag 标号，flag=1,计算 G,DG, flag=2 计算G,DG,Gr,DGr,Gz,DGz

% RG 返回 [G,DG,Gr,DGr,Gz,DGz]实部
% IG 返回 [G,DG,Gr,DGr,Gz,DGz]虚部
% time 返回计算时间

switch flag
    case 1
        tic;
        n=length(omega);
        RG=zeros(n,2);
        IG=zeros(n,2);
        for i=1:n
            RGF=GreenF(r,z,omega(i),flag);
            RG(i,:)=RGF(1,1:2);
            IG(i,:)=RGF(2,1:2);
        end
        time=toc;
    case 2
        tic;
        n=length(omega);
        RG=zeros(n,6);
        IG=zeros(n,6);
        for i=1:n
            RGF=GreenF(r,z,omega(i),flag);
            RG(i,:)=RGF(1,:);
            IG(i,:)=RGF(2,:);
        end
        time=toc;
end

end

