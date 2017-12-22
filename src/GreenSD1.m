function [xx,RIG,Optional,time]=GreenSD1(internal,h,r,z,flag1,flag2)
% 本程序编写于2014年11月19日，硕士论文，格林函数微分方程求解,分别求解 G,DG;Gr,DGr;Gz,DGz;
% 主要实现变步长，当 omega>=1 时， h'=10*h
% 注意，此程序需要保证 omega>1
% internal 频率求解区间，程序中以t代替omega
% h 微分方程计算步长，初始步长
% r,z 位置值

% xx 频率值
%RIG 分别代表格林函数及其一阶导数，格林函数梯度方程及其一阶导数值,第一列为原始值，第二列为导数值
% flag1 标号，flag1=1,计算实部，flag1=2,计算虚部
% flag2 标号，flag2=1,计算 G,DG, flag2=2 计算Gr,DGr,flag2=3,计算 Gz,DGz
% Optional TBNM求解中优化参数
% time 运行时间

if internal(2)<=1
    error('频率必须大于 1 才能进行变步长！');
    return;
end

tic;

% 将区间分开
h1=10*h;
hx=internal(1):h:1;
hn=length(hx);
internal1=[internal(1),hx(hn)];
internal2=[hx(hn), internal(2)];

if r>0 && z<0
    switch flag1 % 实部或虚部
        case 1 % 实部
            switch flag2 % 计算 G DG,Gr DGr,Gz DGz
                case 1 % G DG
                    % 计算初值
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[原始值，导数值] 输出为格林
                    [x1,yR1,OptR1]=GreenTBNM([RG(1),RG(2)], internal1, h, r, z, 1, 1);  % 计算实部
                    [x2,yR2,OptR2]=GreenTBNM([yR1(hn,1),yR1(hn,2)], internal2, h1, r, z, 1, 1);  % 计算实部
                    RIG=[yR1(1:hn-1,:);yR2];  % 两列，第一列为原始值，第二列为导数值
                    Optional=[OptR1,OptR2];
                    xx=[x1(1:hn-1);x2];
                case 2 % Gr DGr
                    % 计算初值
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[原始值，导数值] 输出为格林
                    [x,yR,OptR]=GreenTBNM([RG(1),RG(2)], internal, h, r, z, 2, 1);  % 计算实部
                    RIG=yR;  % 两列，第一列为原始值，第二列为导数值
                    Optional=OptR;
                    xx=x;
                case 3 % Gz DGz
                    % 计算初值
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[原始值，导数值] 输出为格林
                    [x,yR,OptR]=GreenTBNM([RG(1),RG(2)], internal, h, r, z, 3, 1);  % 计算实部
                    RIG=yR;  % 两列，第一列为原始值，第二列为导数值
                    Optional=OptR;
                    xx=x;
            end
        case 2  % 虚部
            switch flag2 % 计算 G DG,Gr DGr,Gz DGz
                case 1 % G DG
                    % 计算初值
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[原始值，导数值] 输出为格林
                    [x,yR,OptR]=GreenTBNM([RG(1),RG(2)], internal, h, r, z, 1, 2);  % 计算虚部
                    RIG=yR;  % 两列，第一列为原始值，第二列为导数值
                    Optional=OptR;
                    xx=x;
                case 2 % Gr DGr
                    % 计算初值
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[原始值，导数值] 输出为格林
                    [x,yR,OptR]=GreenTBNM([RG(1),RG(2)], internal, h, r, z, 2, 2);  % 计算实部
                    RIG=yR;  % 两列，第一列为原始值，第二列为导数值
                    Optional=OptR;
                    xx=x;
                case 3 % Gz DGz
                    % 计算初值
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[原始值，导数值] 输出为格林
                    [x,yR,OptR]=GreenTBNM([RG(1),RG(2)], internal, h, r, z, 3, 2);  % 计算实部
                    RIG=yR;  % 两列，第一列为原始值，第二列为导数值
                    Optional=OptR;
                    xx=x;
            end
    end
else
    warning('r,z 输入有误！');
    return;
end
time=toc;

end



