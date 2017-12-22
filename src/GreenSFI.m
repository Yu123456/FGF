function [RG,time]=GreenSFI(r,z,omega,flag1,flag2)
% 本程序编写于2014年8月24日，通过直接计算频域格林函数积分方程形式计算频域格林函数值
% 分别计算 G DG; Gr DGr; Gz DGz 实部与虚部
% 直接调用函数 RG=GreenSF(r,z,omega,flag1,flag2)

% r,z  场点位置
% omega 频率
% flag1 标号，flag=1,计算实部，flag1=2,计算虚部
% flag2 标号，flag2=1,计算 G,DG, flag2=2 计算Gr,DGr,flag2=3,计算 Gz,DGz

% RG 返回 [原始值，导数值]
% time 返回计算时间

tic;
switch flag1 % 实部或虚部
    case 1 % 实部
        switch flag2 % G DG, Gr DGr, Gz DGz
            case 1 % G DG
                n=length(omega);
                RG=zeros(n,2);
                for i=1:n
                    RGF=GreenSF(r,z,omega(i),1,1);
                    RG(i,:)=RGF(:);
                end
            case 2 % Gr DGr
                n=length(omega);
                RG=zeros(n,2);
                for i=1:n
                    RGF=GreenSF(r,z,omega(i),1,2);
                    RG(i,:)=RGF(:);
                end
            case 3 % Gz DGz
                n=length(omega);
                RG=zeros(n,2);
                for i=1:n
                    RGF=GreenSF(r,z,omega(i),1,3);
                    RG(i,:)=RGF(:);
                end
        end
    case 2 % 虚部
        switch flag2 % G DG, Gr DGr, Gz DGz
            case 1 % G DG
                n=length(omega);
                RG=zeros(n,2);
                for i=1:n
                    RGF=GreenSF(r,z,omega(i),2,1);
                    RG(i,:)=RGF(:);
                end
            case 2 % Gr DGr
                n=length(omega);
                RG=zeros(n,2);
                for i=1:n
                    RGF=GreenSF(r,z,omega(i),2,2);
                    RG(i,:)=RGF(:);
                end
            case 3 % Gz DGz
                n=length(omega);
                RG=zeros(n,2);
                for i=1:n
                    RGF=GreenSF(r,z,omega(i),2,3);
                    RG(i,:)=RGF(:);
                end
        end
end
time=toc;

end

