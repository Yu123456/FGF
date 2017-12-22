function GreenSCom(internal,r,z,h,flag1,flag2)
% �������д��2014��8��28�գ�Green����ֱ�������ַ�����΢�ַ��̾��ȼ�����ʱ��ıȽ�,�ֱ�Ա�
% ���ú��� [xx,RIG,Optional,time]=GreenSD(internal,h,r,z,flag1,flag2)
% ���ú��� [RG,time]=GreenSFI(r,z,omega,flag1,flag2)

% r,z  ����λ��
% internal Ƶ������
% h ����
% flag1 ��ţ�flag=1,����ʵ����flag1=2,�����鲿
% flag2 ��ţ�flag2=1,���� G,DG, flag2=2 ����Gr,DGr,flag2=3,���� Gz,DGz
% �������Σ���ʱֻ���� G DG ʵ��

switch flag1 % ʵ�����鲿
    case 1 % ʵ��
        switch flag2 % G DG,Gr DGr,Gz DGz
            case 1 % G DG
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                [RG,timeI]=GreenSFI(r,z,xx,flag1,flag2);
                error=RIG-RG;
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('G','DG');
                set(hl,'FontSize',10);
                h1=title(['G ʵ����΢�ַ�ʱ��: ',num2str(timeD),'�Ż����� \omega=',num2str(Optional)]);
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
                h1=suptitle({['Green ����ʵ�������ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
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
                h1=title(['Gr ʵ����΢�ַ�ʱ��: ',num2str(timeD),'�Ż����� \omega=',num2str(Optional)]);
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
                h1=suptitle({['Green ����r�ݶ�ʵ�������ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
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
                h1=title(['Gz ʵ����΢�ַ�ʱ��: ',num2str(timeD),'�Ż����� \omega=',num2str(Optional)]);
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
                h1=suptitle({['Green ����z�ݶ�ʵ�������ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
        end
    case 2 % �鲿
        switch flag2 % G DG,Gr DGr,Gz DGz
            case 1 % G DG
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                [RG,timeI]=GreenSFI(r,z,xx,flag1,flag2);
                error=RIG-RG;
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('G','DG');
                set(hl,'FontSize',10);
                h1=title(['G �鲿��΢�ַ�ʱ��: ',num2str(timeD),'�Ż����� \omega=',num2str(Optional)]);
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
                h1=suptitle({['Green �����鲿�����ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
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
                h1=title(['Gr �鲿��΢�ַ�ʱ��: ',num2str(timeD),'�Ż����� \omega=',num2str(Optional)]);
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
                h1=suptitle({['Green ����r�ݶ��鲿�����ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
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
                h1=title(['Gz �鲿��΢�ַ�ʱ��: ',num2str(timeD),'�Ż����� \omega=',num2str(Optional)]);
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
                h1=suptitle({['Green ����z�ݶ��鲿�����ַ�ʱ��: ',num2str(timeI),',΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
        end
end


end

