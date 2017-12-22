function [xx,OG,ODG,OGr,ODGr,OGz,ODGz, Optional]=GreenD(internal,h,r,z,flag)
% �������д��2014��8��13�գ�˶ʿ���ģ����ֺ���΢�ַ������
% internal Ƶ��������䣬��������t����omega
% h ΢�ַ��̼��㲽��
% r,z λ��ֵ

% xx Ƶ��ֵ
% OG,ODG,OGr,ODGr,OGz,ODGz �ֱ������ֺ�������һ�׵��������ֺ����ݶȷ��̼���
% һ�׵���ֵ����һ��Ϊʵ�����ڶ���Ϊ�鲿
% flag ��ţ�flag=1 ��ʾ���� OG, ODG, flag=2 ��ʾ����OG,ODG,OGr,ODGr,OGz,ODGz
% Optional TBNM������Ż�����

Optional=zeros(2,3);
if r>0 && z<0
    switch flag
        case 1 % ����OG, ODG
            % �����ֵ
            RG=GreenF(r,z,internal(1),flag); % RG=[G,DG,Gr,DGr,Gz,DGz] ���Ϊ����
                                             % ���������ݶȷ���ֵ������ Omega һ��
                                             % ����ֵ����һ��Ϊʵ�����ڶ���Ϊ�鲿
            G=RG(:,1); 
            DG=RG(:,2);
            [x,yR,OptR]=GreenTBNM([G(1),DG(1)], internal, h, r, z, 1, 1);  % ����ʵ��
            [x,yI,OptI]=GreenTBNM([G(2),DG(2)], internal, h, r, z, 1, 2);  % �����鲿
            OG=[yR(:,1),yI(:,1)];
            ODG=[yR(:,2),yI(:,2)];
            Optional(:,1)=[OptR; OptI];
            OGr=0;
            ODGr=0;
            OGz=0;
            ODGz=0;
            xx=x;
        case 2 % ����OG,ODG,OGr,ODGr,OGz,ODGz
            % �����ֵ
            RG=GreenF(r,z,internal(1),flag); % RG=[G,DG,Gr,DGr,Gz,DGz] ���Ϊ����
                                             % ���������ݶȷ���ֵ������ Omega һ��
                                             % ����ֵ����һ��Ϊʵ�����ڶ���Ϊ�鲿
            G=RG(:,1); 
            DG=RG(:,2);
            [x,yR,OptR]=GreenTBNM([G(1),DG(1)], internal, h, r, z, 1, 1);  % G����ʵ��
            [x,yI,OptI]=GreenTBNM([G(2),DG(2)], internal, h, r, z, 1, 2);  % G�����鲿
            OG=[yR(:,1),yI(:,1)];
            ODG=[yR(:,2),yI(:,2)];
            Optional(:,1)=[OptR; OptI];
            Gr=RG(:,3); 
            DGr=RG(:,4);
            [x,yR,OptR]=GreenTBNM([Gr(1),DGr(1)], internal, h, r, z, 2, 1);  % Gr����ʵ��
            [x,yI,OptI]=GreenTBNM([Gr(2),DGr(2)], internal, h, r, z, 2, 2);  % Gr�����鲿
            OGr=[yR(:,1),yI(:,1)];
            ODGr=[yR(:,2),yI(:,2)];
            Optional(:,2)=[OptR; OptI];
            Gz=RG(:,5); 
            DGz=RG(:,6);
            [x,yR,OptR]=GreenTBNM([Gz(1),DGz(1)], internal, h, r, z, 3, 1);  % Gz����ʵ��
            [x,yI,OptI]=GreenTBNM([Gz(2),DGz(2)], internal, h, r, z, 3, 2);  % Gz�����鲿
            OGz=[yR(:,1),yI(:,1)];
            ODGz=[yR(:,2),yI(:,2)];
            Optional(:,3)=[OptR; OptI];
            xx=x;
    end
else
    warning('r,z ��������');
    return;
end

end

