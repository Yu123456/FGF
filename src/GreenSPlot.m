function GreenSPlot(internal,r,z,h,flag1,flag2)
% �������д��2014��8��30�գ�Green����΢�ַ��̷ֱ���⣬��ͼ
% ���ú��� [xx,RIG,Optional,time]=GreenSD(internal,h,r,z,flag1,flag2)


% r,z  ����λ��
% internal Ƶ������
% h ����
% flag1 ��ţ�flag=1,����ʵ����flag1=2,�����鲿
% flag2 ��ţ�flag2=1,���� G,DG, flag2=2 ����Gr,DGr,flag2=3,���� Gz,DGz


switch flag1 % ʵ�����鲿
    case 1 % ʵ��
        switch flag2 % G DG,Gr DGr,Gz DGz
            case 1 % G DG
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('G','DG');
                set(hl,'FontSize',10);
                h1=title({['Green ����ʵ��,΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
            case 2 % Gr DGr
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('Gr','DGr');
                set(hl,'FontSize',10);
                h1=title({['Green ����r�ݶ�ʵ��,΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
            case 3 % Gz DGz
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('Gz','DGz');
                set(hl,'FontSize',10);
                h1=title({['Green ����z�ݶ�ʵ��,΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
        end
    case 2 % �鲿
        switch flag2 % G DG,Gr DGr,Gz DGz
            case 1 % G DG
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('G','DG');
                set(hl,'FontSize',10);
                h1=title({['Green �����鲿,΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
            case 2 % Gr DGr
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('Gr','DGr');
                set(hl,'FontSize',10);
                h1=title({['Green ����r�ݶ�,΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
            case 3 % Gz DGz
                [xx,RIG,Optional,timeD]=GreenSD(internal,h,r,z,flag1,flag2);
                figure;
                plot(xx,RIG(:,1),'-r',xx,RIG(:,2),'--k');
                hl=legend('Gz','DGz');
                set(hl,'FontSize',10);
                h1=title({['Green ����z�ݶ��鲿,΢�ַ�ʱ��: ',num2str(timeD)];['(r,z)=( ',num2str(r),',',num2str(z),'), h=',num2str(h),...
                    ',\omega =',num2str(Optional)]});
                set(h1,'FontSize',10);
        end
end


end

