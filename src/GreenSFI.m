function [RG,time]=GreenSFI(r,z,omega,flag1,flag2)
% �������д��2014��8��24�գ�ͨ��ֱ�Ӽ���Ƶ����ֺ������ַ�����ʽ����Ƶ����ֺ���ֵ
% �ֱ���� G DG; Gr DGr; Gz DGz ʵ�����鲿
% ֱ�ӵ��ú��� RG=GreenSF(r,z,omega,flag1,flag2)

% r,z  ����λ��
% omega Ƶ��
% flag1 ��ţ�flag=1,����ʵ����flag1=2,�����鲿
% flag2 ��ţ�flag2=1,���� G,DG, flag2=2 ����Gr,DGr,flag2=3,���� Gz,DGz

% RG ���� [ԭʼֵ������ֵ]
% time ���ؼ���ʱ��

tic;
switch flag1 % ʵ�����鲿
    case 1 % ʵ��
        switch flag2 % G DG, Gr DGr, Gz DGz
            case 1 % G DG
                n=length(omega);
                RG=zeros(n,2);
                for i=1:n
                    RGF=GreenSF(r,z,omega(i),1,1);
                    RG(i,:)=RGF(:);
                end
            case 2 % Gr DGr
                n=length(omega);
                RG=zeros(n,2);
                for i=1:n
                    RGF=GreenSF(r,z,omega(i),1,2);
                    RG(i,:)=RGF(:);
                end
            case 3 % Gz DGz
                n=length(omega);
                RG=zeros(n,2);
                for i=1:n
                    RGF=GreenSF(r,z,omega(i),1,3);
                    RG(i,:)=RGF(:);
                end
        end
    case 2 % �鲿
        switch flag2 % G DG, Gr DGr, Gz DGz
            case 1 % G DG
                n=length(omega);
                RG=zeros(n,2);
                for i=1:n
                    RGF=GreenSF(r,z,omega(i),2,1);
                    RG(i,:)=RGF(:);
                end
            case 2 % Gr DGr
                n=length(omega);
                RG=zeros(n,2);
                for i=1:n
                    RGF=GreenSF(r,z,omega(i),2,2);
                    RG(i,:)=RGF(:);
                end
            case 3 % Gz DGz
                n=length(omega);
                RG=zeros(n,2);
                for i=1:n
                    RGF=GreenSF(r,z,omega(i),2,3);
                    RG(i,:)=RGF(:);
                end
        end
end
time=toc;

end

