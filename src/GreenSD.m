function [xx,RIG,Optional,time]=GreenSD(internal,h,r,z,flag1,flag2)
% �������д��2014��8��28�գ�˶ʿ���ģ����ֺ���΢�ַ������,�ֱ���� G,DG;Gr,DGr;Gz,DGz;
% internal Ƶ��������䣬��������t����omega
% h ΢�ַ��̼��㲽��
% r,z λ��ֵ

% xx Ƶ��ֵ
%RIG �ֱ������ֺ�������һ�׵��������ֺ����ݶȷ��̼���һ�׵���ֵ,��һ��Ϊԭʼֵ���ڶ���Ϊ����ֵ
% flag1 ��ţ�flag1=1,����ʵ����flag1=2,�����鲿
% flag2 ��ţ�flag2=1,���� G,DG, flag2=2 ����Gr,DGr,flag2=3,���� Gz,DGz
% Optional TBNM������Ż�����
% time ����ʱ��

tic;
if r>0 && z<0
    switch flag1 % ʵ�����鲿
        case 1 % ʵ��
            switch flag2 % ���� G DG,Gr DGr,Gz DGz
                case 1 % G DG
                    % �����ֵ
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[ԭʼֵ������ֵ] ���Ϊ����
                    [x,yR,OptR]=GreenTBNM([RG(1),RG(2)], internal, h, r, z, 1, 1);  % ����ʵ��
                    RIG=yR;  % ���У���һ��Ϊԭʼֵ���ڶ���Ϊ����ֵ
                    Optional=OptR;
                    xx=x;
                case 2 % Gr DGr
                    % �����ֵ
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[ԭʼֵ������ֵ] ���Ϊ����
                    [x,yR,OptR]=GreenTBNM([RG(1),RG(2)], internal, h, r, z, 2, 1);  % ����ʵ��
                    RIG=yR;  % ���У���һ��Ϊԭʼֵ���ڶ���Ϊ����ֵ
                    Optional=OptR;
                    xx=x;
                case 3 % Gz DGz
                    % �����ֵ
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[ԭʼֵ������ֵ] ���Ϊ����
                    [x,yR,OptR]=GreenTBNM([RG(1),RG(2)], internal, h, r, z, 3, 1);  % ����ʵ��
                    RIG=yR;  % ���У���һ��Ϊԭʼֵ���ڶ���Ϊ����ֵ
                    Optional=OptR;
                    xx=x;
            end
        case 2  % �鲿
            switch flag2 % ���� G DG,Gr DGr,Gz DGz
                case 1 % G DG
                    % �����ֵ
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[ԭʼֵ������ֵ] ���Ϊ����
                    [x,yR,OptR]=GreenTBNM([RG(1),RG(2)], internal, h, r, z, 1, 2);  % �����鲿
                    RIG=yR;  % ���У���һ��Ϊԭʼֵ���ڶ���Ϊ����ֵ
                    Optional=OptR;
                    xx=x;
                case 2 % Gr DGr
                    % �����ֵ
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[ԭʼֵ������ֵ] ���Ϊ����
                    [x,yR,OptR]=GreenTBNM([RG(1),RG(2)], internal, h, r, z, 2, 2);  % ����ʵ��
                    RIG=yR;  % ���У���һ��Ϊԭʼֵ���ڶ���Ϊ����ֵ
                    Optional=OptR;
                    xx=x;
                case 3 % Gz DGz
                    % �����ֵ
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[ԭʼֵ������ֵ] ���Ϊ����
                    [x,yR,OptR]=GreenTBNM([RG(1),RG(2)], internal, h, r, z, 3, 2);  % ����ʵ��
                    RIG=yR;  % ���У���һ��Ϊԭʼֵ���ڶ���Ϊ����ֵ
                    Optional=OptR;
                    xx=x;
            end
    end
else
    warning('r,z ��������');
    return;
end
time=toc;

end

