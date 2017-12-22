function GreenSCom(internal,r,z,h,flag1,flag2)
% 本程序编写于2014年8月28日，Green函数直接求解积分方程与微分方程精度及运行时间的比较,分别对比
% 调用函数 [xx,RIG,Optional,time]=GreenSD(internal,h,r,z,flag1,flag2)
% 调用函数 [RG,time]=GreenSFI(r,z,omega,flag1,flag2)

% r,z  场点位置
% internal 频率区间
% h 步长
% flag1 标号，flag=1,计算实部，flag1=2,计算虚部
% flag2 标号，flag2=1,计算 G,DG, flag2=2 计算Gr,DGr,flag2=3,计算 Gz,DGz
% 测试情形，此时只计算 G DG 实部

switch flag1 % 实部或虚部
    case 1 % 实部
        switch flag2 % G DG,Gr DGr,Gz DGz
            case 1 % G DG
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                [RG,timeI]=GreenSFI(r,z,xx,flag1,flag2);
                error=RIG-RG;
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('G','DG');
                set(hl,'FontSize',10);
                h1=title(['G 实部，微分法时间: ',num2str(timeD),'优化参数 \omega=',num2str(Optional)]);
                set(h1,'FontSize',10);
                figure;
                subplot(2,1,1);
                plot(xx,error(:,1));
                hl=legend('G');
                set(hl,'FontSize',10);
                subplot(2,1,2);
                plot(xx,error(:,2));
                hl=legend('DG');
                set(hl,'FontSize',10);
                h1=suptitle({['Green 函数实部误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
            case 2 % Gr DGr
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                [RG,timeI]=GreenSFI(r,z,xx,flag1,flag2);
                error=RIG-RG;
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('Gr','DGr');
                set(hl,'FontSize',10);
                h1=title(['Gr 实部，微分法时间: ',num2str(timeD),'优化参数 \omega=',num2str(Optional)]);
                set(h1,'FontSize',10);
                figure;
                subplot(2,1,1);
                plot(xx,error(:,1));
                hl=legend('Gr');
                set(hl,'FontSize',10);
                subplot(2,1,2);
                plot(xx,error(:,2));
                hl=legend('DGr');
                set(hl,'FontSize',10);
                h1=suptitle({['Green 函数r梯度实部误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
            case 3 % Gz DGz
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                [RG,timeI]=GreenSFI(r,z,xx,flag1,flag2);
                error=RIG-RG;
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('Gz','DGz');
                set(hl,'FontSize',10);
                h1=title(['Gz 实部，微分法时间: ',num2str(timeD),'优化参数 \omega=',num2str(Optional)]);
                set(h1,'FontSize',10);
                figure;
                subplot(2,1,1);
                plot(xx,error(:,1));
                hl=legend('Gr');
                set(hl,'FontSize',10);
                subplot(2,1,2);
                plot(xx,error(:,2));
                hl=legend('DGr');
                set(hl,'FontSize',10);
                h1=suptitle({['Green 函数z梯度实部误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
        end
    case 2 % 虚部
        switch flag2 % G DG,Gr DGr,Gz DGz
            case 1 % G DG
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                [RG,timeI]=GreenSFI(r,z,xx,flag1,flag2);
                error=RIG-RG;
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('G','DG');
                set(hl,'FontSize',10);
                h1=title(['G 虚部，微分法时间: ',num2str(timeD),'优化参数 \omega=',num2str(Optional)]);
                set(h1,'FontSize',10);
                figure;
                subplot(2,1,1);
                plot(xx,error(:,1));
                hl=legend('G');
                set(hl,'FontSize',10);
                subplot(2,1,2);
                plot(xx,error(:,2));
                hl=legend('DG');
                set(hl,'FontSize',10);
                h1=suptitle({['Green 函数虚部误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
            case 2 % Gr DGr
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                [RG,timeI]=GreenSFI(r,z,xx,flag1,flag2);
                error=RIG-RG;
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('Gr','DGr');
                set(hl,'FontSize',10);
                h1=title(['Gr 虚部，微分法时间: ',num2str(timeD),'优化参数 \omega=',num2str(Optional)]);
                set(h1,'FontSize',10);
                figure;
                subplot(2,1,1);
                plot(xx,error(:,1));
                hl=legend('Gr');
                set(hl,'FontSize',10);
                subplot(2,1,2);
                plot(xx,error(:,2));
                hl=legend('DGr');
                set(hl,'FontSize',10);
                h1=suptitle({['Green 函数r梯度虚部误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
            case 3 % Gz DGz
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                [RG,timeI]=GreenSFI(r,z,xx,flag1,flag2);
                error=RIG-RG;
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('Gz','DGz');
                set(hl,'FontSize',10);
                h1=title(['Gz 虚部，微分法时间: ',num2str(timeD),'优化参数 \omega=',num2str(Optional)]);
                set(h1,'FontSize',10);
                figure;
                subplot(2,1,1);
                plot(xx,error(:,1));
                hl=legend('Gr');
                set(hl,'FontSize',10);
                subplot(2,1,2);
                plot(xx,error(:,2));
                hl=legend('DGr');
                set(hl,'FontSize',10);
                h1=suptitle({['Green 函数z梯度虚部误差，积分法时间: ',num2str(timeI),',微分法时间: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
        end
end


end

