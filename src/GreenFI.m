function [RG,IG,time]=GreenFI(r,z,omega,flag)
% �������д��2014��8��24�գ�ͨ��ֱ�Ӽ���Ƶ����ֺ������ַ�����ʽ����Ƶ����ֺ���ֵ
% ֱ�ӵ��ú��� RG=GreenF(r,z,omega,flag)

% r,z  ����λ��
% omega Ƶ��
% flag ��ţ�flag=1,���� G,DG, flag=2 ����G,DG,Gr,DGr,Gz,DGz

% RG ���� [G,DG,Gr,DGr,Gz,DGz]ʵ��
% IG ���� [G,DG,Gr,DGr,Gz,DGz]�鲿
% time ���ؼ���ʱ��

switch flag
    case 1
        tic;
        n=length(omega);
        RG=zeros(n,2);
        IG=zeros(n,2);
        for i=1:n
            RGF=GreenF(r,z,omega(i),flag);
            RG(i,:)=RGF(1,1:2);
            IG(i,:)=RGF(2,1:2);
        end
        time=toc;
    case 2
        tic;
        n=length(omega);
        RG=zeros(n,6);
        IG=zeros(n,6);
        for i=1:n
            RGF=GreenF(r,z,omega(i),flag);
            RG(i,:)=RGF(1,:);
            IG(i,:)=RGF(2,:);
        end
        time=toc;
end

end

