function [xx,OG,ODG,OGr,ODGr,OGz,ODGz, Optional]=GreenD(internal,h,r,z,flag)
% 本程序编写于2014年8月13日，硕士论文，格林函数微分方程求解
% internal 频率求解区间，程序中以t代替omega
% h 微分方程计算步长
% r,z 位置值

% xx 频率值
% OG,ODG,OGr,ODGr,OGz,ODGz 分别代表格林函数及其一阶导数，格林函数梯度方程及其
% 一阶导数值，第一列为实部，第二列为虚部
% flag 标号，flag=1 表示计算 OG, ODG, flag=2 表示计算OG,ODG,OGr,ODGr,OGz,ODGz
% Optional TBNM求解中优化参数

Optional=zeros(2,3);
if r>0 && z<0
    switch flag
        case 1 % 计算OG, ODG
            % 计算初值
            RG=GreenF(r,z,internal(1),flag); % RG=[G,DG,Gr,DGr,Gz,DGz] 输出为格林
                                             % 函数及其梯度方程值，及对 Omega 一阶
                                             % 导数值，第一行为实部，第二行为虚部
            G=RG(:,1); 
            DG=RG(:,2);
            [x,yR,OptR]=GreenTBNM([G(1),DG(1)], internal, h, r, z, 1, 1);  % 计算实部
            [x,yI,OptI]=GreenTBNM([G(2),DG(2)], internal, h, r, z, 1, 2);  % 计算虚部
            OG=[yR(:,1),yI(:,1)];
            ODG=[yR(:,2),yI(:,2)];
            Optional(:,1)=[OptR; OptI];
            OGr=0;
            ODGr=0;
            OGz=0;
            ODGz=0;
            xx=x;
        case 2 % 计算OG,ODG,OGr,ODGr,OGz,ODGz
            % 计算初值
            RG=GreenF(r,z,internal(1),flag); % RG=[G,DG,Gr,DGr,Gz,DGz] 输出为格林
                                             % 函数及其梯度方程值，及对 Omega 一阶
                                             % 导数值，第一行为实部，第二行为虚部
            G=RG(:,1); 
            DG=RG(:,2);
            [x,yR,OptR]=GreenTBNM([G(1),DG(1)], internal, h, r, z, 1, 1);  % G计算实部
            [x,yI,OptI]=GreenTBNM([G(2),DG(2)], internal, h, r, z, 1, 2);  % G计算虚部
            OG=[yR(:,1),yI(:,1)];
            ODG=[yR(:,2),yI(:,2)];
            Optional(:,1)=[OptR; OptI];
            Gr=RG(:,3); 
            DGr=RG(:,4);
            [x,yR,OptR]=GreenTBNM([Gr(1),DGr(1)], internal, h, r, z, 2, 1);  % Gr计算实部
            [x,yI,OptI]=GreenTBNM([Gr(2),DGr(2)], internal, h, r, z, 2, 2);  % Gr计算虚部
            OGr=[yR(:,1),yI(:,1)];
            ODGr=[yR(:,2),yI(:,2)];
            Optional(:,2)=[OptR; OptI];
            Gz=RG(:,5); 
            DGz=RG(:,6);
            [x,yR,OptR]=GreenTBNM([Gz(1),DGz(1)], internal, h, r, z, 3, 1);  % Gz计算实部
            [x,yI,OptI]=GreenTBNM([Gz(2),DGz(2)], internal, h, r, z, 3, 2);  % Gz计算虚部
            OGz=[yR(:,1),yI(:,1)];
            ODGz=[yR(:,2),yI(:,2)];
            Optional(:,3)=[OptR; OptI];
            xx=x;
    end
else
    warning('r,z 输入有误！');
    return;
end

end

