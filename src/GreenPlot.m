function GreenPlot(internal,h,r,z,flag)
% 本程序编写于2014年8月13日，Green函数作图，硕士论文
% internal 频率求解区间
% h 微分方程计算步长
% r,z 位置值
% flag 标号，flag=1 表示计算 OG, ODG, flag=2 表示计算OG,ODG,OGr,ODGr,OGz,ODGz

switch flag
    case 1
        [xx,OG,ODG,OGr,ODGr,OGz,ODGz, Optional]=GreenD(internal,h,r,z,flag);
        % OG,ODG,OGr,ODGr,OGz,ODGz 分别代表格林函数及其一阶导数，格林函数梯度方程及其
        % 一阶导数值，第一列为实部，第二列为虚部
        % Optional TBNM求解中优化参数
        figure;
        plot(xx,OG(:,1));
        title(['G 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(1,1))]);
        legend('Real G');
        xlabel('频率\omega');
        figure;
        plot(xx,OG(:,2));
        title(['G 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(2,1))]);
        legend('Image G');
        xlabel('频率\omega');
        figure;
        plot(xx,ODG(:,1));
        title(['DG 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(1,1))]);
        legend('Real DG');
        xlabel('频率\omega');
        figure;
        plot(xx,ODG(:,2));
        title(['DG 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(2,1))]);
        legend('Image DG');
        xlabel('频率\omega');
        figure;
        plot(xx,OG(:,1),'--r',xx,OG(:,2),'-k');
        title(['G 函数频率图，(r,z)=(',num2str(r),',',num2str(z),')']);
        legend('Real G','Image G');
        xlabel('频率\omega');
        figure;
        plot(xx,ODG(:,1),'--r',xx,ODG(:,2),'-k');
        title(['DG 函数频率图，(r,z)=(',num2str(r),',',num2str(z),')']);
        legend('Real DG','Image DG');
        xlabel('频率\omega');
    case 2
        [xx,OG,ODG,OGr,ODGr,OGz,ODGz, Optional]=GreenD(internal,h,r,z,flag);
        % OG,ODG,OGr,ODGr,OGz,ODGz 分别代表格林函数及其一阶导数，格林函数梯度方程及其
        % 一阶导数值，第一列为实部，第二列为虚部
        % Optional TBNM求解中优化参数
        figure;
        plot(xx,OG(:,1));
        title(['G 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(1,1))]);
        legend('Real G');
        xlabel('频率\omega');
        figure;
        plot(xx,OG(:,2));
        title(['G 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(2,1))]);
        legend('Image G');
        xlabel('频率\omega');
        figure;
        plot(xx,OGr(:,1));
        title(['Gr 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(1,2))]);
        legend('Real Gr');
        xlabel('频率\omega');
        figure;
        plot(xx,OGr(:,2));
        title(['Gr 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(2,2))]);
        legend('Image Gr');
        xlabel('频率\omega');
        figure;
        plot(xx,OGz(:,1));
        title(['Gz 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(1,3))]);
        legend('Real Gz');
        xlabel('频率\omega');
        figure;
        plot(xx,OGz(:,2));
        title(['Gz 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(2,3))]);
        legend('Image Gz');
        xlabel('频率\omega');
        figure;
        plot(xx,ODG(:,1));
        title(['DG 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(1,1))]);
        legend('Real DG');
        xlabel('频率\omega');
        figure;
        plot(xx,ODG(:,2));
        title(['DG 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(2,1))]);
        legend('Image DG');
        xlabel('频率\omega');
        figure;
        plot(xx,ODGr(:,1));
        title(['DGr 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(1,2))]);
        legend('Real DGr');
        xlabel('频率\omega');
        figure;
        plot(xx,ODGr(:,2));
        title(['DGr 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(2,2))]);
        legend('Image DGr');
        xlabel('频率\omega');
        figure;
        plot(xx,ODGz(:,1));
        title(['DGz 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(1,3))]);
        legend('Real DGz');
        xlabel('频率\omega');
        figure;
        plot(xx,ODGz(:,2));
        title(['DGz 函数频率图，(r,z)=(',num2str(r),',',num2str(z),'),微分优化参数\omega=',num2str(Optional(2,3))]);
        legend('Image DGz');
        xlabel('频率\omega');
        figure;
        plot(xx,OG(:,1),'--r',xx,OGr(:,1),'-k',xx,OGz(:,1),'-.b');
        title(['G,Gr,Gz 函数频率图，(r,z)=(',num2str(r),',',num2str(z),')']);
        legend('Real G','Real Gr','Real Gz');
        xlabel('频率\omega');
        figure;
        plot(xx,OG(:,2),'--r',xx,OGr(:,2),'-k',xx,OGz(:,2),'-.b');
        title(['G,Gr,Gz 函数频率图，(r,z)=(',num2str(r),',',num2str(z),')']);
        legend('Image G','Image Gr','Image Gz');
        xlabel('频率\omega');
        figure;
        plot(xx,ODG(:,1),'--r',xx,ODGr(:,1),'-k',xx,ODGz(:,1),'-.b');
        title(['DG,DGr,DGz 函数频率图，(r,z)=(',num2str(r),',',num2str(z),')']);
        legend('Real DG','Real DGr','Real DGz');
        xlabel('频率\omega');
        figure;
        plot(xx,ODG(:,2),'--r',xx,ODGr(:,2),'-k',xx,ODGz(:,2),'-.b');
        title(['DG,DGr,DGz 函数频率图，(r,z)=(',num2str(r),',',num2str(z),')']);
        legend('Image DG','Image DGr','Image DGz');
        xlabel('频率\omega');
end
end

