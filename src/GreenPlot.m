function GreenPlot(internal,h,r,z,flag)
% �������д��2014��8��13�գ�Green������ͼ��˶ʿ����
% internal Ƶ���������
% h ΢�ַ��̼��㲽��
% r,z λ��ֵ
% flag ��ţ�flag=1 ��ʾ���� OG, ODG, flag=2 ��ʾ����OG,ODG,OGr,ODGr,OGz,ODGz

switch flag
    case 1
        [xx,OG,ODG,OGr,ODGr,OGz,ODGz, Optional]=GreenD(internal,h,r,z,flag);
        % OG,ODG,OGr,ODGr,OGz,ODGz �ֱ������ֺ�������һ�׵��������ֺ����ݶȷ��̼���
        % һ�׵���ֵ����һ��Ϊʵ�����ڶ���Ϊ�鲿
        % Optional TBNM������Ż�����
        figure;
        plot(xx,OG(:,1));
        title(['G ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(1,1))]);
        legend('Real G');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,OG(:,2));
        title(['G ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(2,1))]);
        legend('Image G');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,ODG(:,1));
        title(['DG ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(1,1))]);
        legend('Real DG');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,ODG(:,2));
        title(['DG ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(2,1))]);
        legend('Image DG');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,OG(:,1),'--r',xx,OG(:,2),'-k');
        title(['G ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),')']);
        legend('Real G','Image G');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,ODG(:,1),'--r',xx,ODG(:,2),'-k');
        title(['DG ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),')']);
        legend('Real DG','Image DG');
        xlabel('Ƶ��\omega');
    case 2
        [xx,OG,ODG,OGr,ODGr,OGz,ODGz, Optional]=GreenD(internal,h,r,z,flag);
        % OG,ODG,OGr,ODGr,OGz,ODGz �ֱ������ֺ�������һ�׵��������ֺ����ݶȷ��̼���
        % һ�׵���ֵ����һ��Ϊʵ�����ڶ���Ϊ�鲿
        % Optional TBNM������Ż�����
        figure;
        plot(xx,OG(:,1));
        title(['G ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(1,1))]);
        legend('Real G');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,OG(:,2));
        title(['G ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(2,1))]);
        legend('Image G');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,OGr(:,1));
        title(['Gr ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(1,2))]);
        legend('Real Gr');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,OGr(:,2));
        title(['Gr ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(2,2))]);
        legend('Image Gr');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,OGz(:,1));
        title(['Gz ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(1,3))]);
        legend('Real Gz');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,OGz(:,2));
        title(['Gz ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(2,3))]);
        legend('Image Gz');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,ODG(:,1));
        title(['DG ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(1,1))]);
        legend('Real DG');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,ODG(:,2));
        title(['DG ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(2,1))]);
        legend('Image DG');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,ODGr(:,1));
        title(['DGr ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(1,2))]);
        legend('Real DGr');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,ODGr(:,2));
        title(['DGr ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(2,2))]);
        legend('Image DGr');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,ODGz(:,1));
        title(['DGz ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(1,3))]);
        legend('Real DGz');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,ODGz(:,2));
        title(['DGz ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),'),΢���Ż�����\omega=',num2str(Optional(2,3))]);
        legend('Image DGz');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,OG(:,1),'--r',xx,OGr(:,1),'-k',xx,OGz(:,1),'-.b');
        title(['G,Gr,Gz ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),')']);
        legend('Real G','Real Gr','Real Gz');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,OG(:,2),'--r',xx,OGr(:,2),'-k',xx,OGz(:,2),'-.b');
        title(['G,Gr,Gz ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),')']);
        legend('Image G','Image Gr','Image Gz');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,ODG(:,1),'--r',xx,ODGr(:,1),'-k',xx,ODGz(:,1),'-.b');
        title(['DG,DGr,DGz ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),')']);
        legend('Real DG','Real DGr','Real DGz');
        xlabel('Ƶ��\omega');
        figure;
        plot(xx,ODG(:,2),'--r',xx,ODGr(:,2),'-k',xx,ODGz(:,2),'-.b');
        title(['DG,DGr,DGz ����Ƶ��ͼ��(r,z)=(',num2str(r),',',num2str(z),')']);
        legend('Image DG','Image DGr','Image DGz');
        xlabel('Ƶ��\omega');
end
end

