function [xx,RG,IG,time,optional]=GreenDD(internal,r,z,h,flag)
% �������д��2014��8��24�գ�ͨ������Ƶ����ֺ���΢�ַ�����ʽ����Ƶ����ֺ���ֵ
% ֱ�ӵ��ú��� [xx,OG,ODG,OGr,ODGr,OGz,ODGz,Optional]=GreenD(internal,h,r,z,flag)

% r,z  ����λ��
% internal Ƶ������
% h ����
% flag ��ţ�flag=1,���� G,DG, flag=2 ����G,DG,Gr,DGr,Gz,DGz

% RG ���� [G,DG,Gr,DGr,Gz,DGz]ʵ��
% IG ���� [G,DG,Gr,DGr,Gz,DGz]�鲿
% time ���ؼ���ʱ��

switch flag
    case 1
        tic;
        [xx,OG,ODG,~,~,~,~,optional]=GreenD(internal,h,r,z,flag);
        RG=[OG(:,1),ODG(:,1)];
        IG=[OG(:,2),ODG(:,2)];
        time=toc;
    case 2
        tic;
        [xx,OG,ODG,OGr,ODGr,OGz,ODGz,optional]=GreenD(internal,h,r,z,flag);
        RG=[OG(:,1),ODG(:,1),OGr(:,1),ODGr(:,1),OGz(:,1),ODGz(:,1)];
        IG=[OG(:,2),ODG(:,2),OGr(:,2),ODGr(:,2),OGz(:,2),ODGz(:,2)];
        time=toc;
end

end

