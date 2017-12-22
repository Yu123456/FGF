function GreenCom(internal,r,z,h,flag)
% 本程序编写于2014年8月24日，Green函数直接求解积分方程与微分方程精度及运行时间的比较
% 调用函数 [xx,RG,IG,time]=GreenDD(internal,r,z,h,flag)
% 调用函数 [RG,IG,time]=GreenFI(r,z,omega,flag)

% r,z  场点位置
% internal 频率区间
% h 步长
% flag 标号，flag=1,计算 G,DG, flag=2 计算G,DG,Gr,DGr,Gz,DGz，flag=3，这是一个
% 测试情形，此时只计算 G DG 实部

switch flag
    case 1
        [xx,RGD,IGD,timeD,optional]=GreenDD(internal,r,z,h,flag);
        [RGI,IGI,timeI]=GreenFI(r,z,xx,flag);
        errorR=RGD-RGI;
        errorI=IGD-IGI;
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,1));
        hl=legend('G 实部');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,1));
        hl=legend('G 虚部');
        set(hl,'FontSize',10);
        h1=suptitle({['Green 函数误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,1))]});
        set(h1,'FontSize',10);
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,2));
        hl=legend('G 导数实部');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,2));
        hl=legend('G 导数虚部');
        set(hl,'FontSize',10);
        h2=suptitle({['Green 函数导数误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,1))]});
        set(h2,'FontSize',10);
    case 2
        [xx,RGD,IGD,timeD,optional]=GreenDD(internal,r,z,h,flag);
        [RGI,IGI,timeI]=GreenFI(r,z,xx,flag);
        errorR=RGD-RGI;
        errorI=IGD-IGI;
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,1));
        hl=legend('G 实部');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,1));
        hl=legend('G 虚部');
        set(hl,'FontSize',10);
        h1=suptitle({['Green 函数误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,1))]});
        set(h1,'FontSize',10);
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,2));
        hl=legend('G 导数实部');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,2));
        hl=legend('G 导数虚部');
        set(hl,'FontSize',10);
        h2=suptitle({['Green 函数导数误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,1))]});
        set(h2,'FontSize',10);
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,3));
        hl=legend('Gr 实部');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,3));
        hl=legend('Gr 虚部');
        set(hl,'FontSize',10);
        h1=suptitle({['Green r 梯度误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,2))]});
        set(h1,'FontSize',10);
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,4));
        hl=legend('G r 梯度导数实部');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,4));
        hl=legend('G r 梯度导数虚部');
        set(hl,'FontSize',10);
        h2=suptitle({['Green r 梯度导数误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,2))]});
        set(h2,'FontSize',10);
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,5));
        hl=legend('Gz 实部');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,5));
        hl=legend('Gz 虚部');
        set(hl,'FontSize',10);
        h1=suptitle({['Green z 梯度误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,3))]});
        set(h1,'FontSize',10);
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,6));
        hl=legend('G z 梯度导数实部');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,6));
        hl=legend('G z 梯度导数虚部');
        set(hl,'FontSize',10);
        h2=suptitle({['Green z 梯度导数误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,3))]});
        set(h2,'FontSize',10);
end

end

