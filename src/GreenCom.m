function GreenCom(internal,r,z,h,flag)
% �������д��2014��8��24�գ�Green����ֱ�������ַ�����΢�ַ��̾��ȼ�����ʱ��ıȽ�
% ���ú��� [xx,RG,IG,time]=GreenDD(internal,r,z,h,flag)
% ���ú��� [RG,IG,time]=GreenFI(r,z,omega,flag)

% r,z  ����λ��
% internal Ƶ������
% h ����
% flag ��ţ�flag=1,���� G,DG, flag=2 ����G,DG,Gr,DGr,Gz,DGz��flag=3������һ��
% �������Σ���ʱֻ���� G DG ʵ��

switch flag
    case 1
        [xx,RGD,IGD,timeD,optional]=GreenDD(internal,r,z,h,flag);
        [RGI,IGI,timeI]=GreenFI(r,z,xx,flag);
        errorR=RGD-RGI;
        errorI=IGD-IGI;
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,1));
        hl=legend('G ʵ��');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,1));
        hl=legend('G �鲿');
        set(hl,'FontSize',10);
        h1=suptitle({['Green ���������ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,1))]});
        set(h1,'FontSize',10);
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,2));
        hl=legend('G ����ʵ��');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,2));
        hl=legend('G �����鲿');
        set(hl,'FontSize',10);
        h2=suptitle({['Green �������������ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
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
        hl=legend('G ʵ��');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,1));
        hl=legend('G �鲿');
        set(hl,'FontSize',10);
        h1=suptitle({['Green ���������ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,1))]});
        set(h1,'FontSize',10);
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,2));
        hl=legend('G ����ʵ��');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,2));
        hl=legend('G �����鲿');
        set(hl,'FontSize',10);
        h2=suptitle({['Green �������������ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,1))]});
        set(h2,'FontSize',10);
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,3));
        hl=legend('Gr ʵ��');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,3));
        hl=legend('Gr �鲿');
        set(hl,'FontSize',10);
        h1=suptitle({['Green r �ݶ������ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,2))]});
        set(h1,'FontSize',10);
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,4));
        hl=legend('G r �ݶȵ���ʵ��');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,4));
        hl=legend('G r �ݶȵ����鲿');
        set(hl,'FontSize',10);
        h2=suptitle({['Green r �ݶȵ��������ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,2))]});
        set(h2,'FontSize',10);
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,5));
        hl=legend('Gz ʵ��');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,5));
        hl=legend('Gz �鲿');
        set(hl,'FontSize',10);
        h1=suptitle({['Green z �ݶ������ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,3))]});
        set(h1,'FontSize',10);
        figure;
        subplot(2,1,1);
        plot(xx,errorR(:,6));
        hl=legend('G z �ݶȵ���ʵ��');
        set(hl,'FontSize',10);
        subplot(2,1,2);
        plot(xx,errorI(:,6));
        hl=legend('G z �ݶȵ����鲿');
        set(hl,'FontSize',10);
        h2=suptitle({['Green z �ݶȵ��������ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
            ',\omega =',num2str(optional(1,3))]});
        set(h2,'FontSize',10);
end

end

