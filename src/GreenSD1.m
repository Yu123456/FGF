function [xx,RIG,Optional,time]=GreenSD1(internal,h,r,z,flag1,flag2)
% �������д��2014��11��19�գ�˶ʿ���ģ����ֺ���΢�ַ������,�ֱ���� G,DG;Gr,DGr;Gz,DGz;
% ��Ҫʵ�ֱ䲽������ omega>=1 ʱ�� h'=10*h
% ע�⣬�˳�����Ҫ��֤ omega>1
% internal Ƶ��������䣬��������t����omega
% h ΢�ַ��̼��㲽������ʼ����
% r,z λ��ֵ

% xx Ƶ��ֵ
%RIG �ֱ������ֺ�������һ�׵��������ֺ����ݶȷ��̼���һ�׵���ֵ,��һ��Ϊԭʼֵ���ڶ���Ϊ����ֵ
% flag1 ��ţ�flag1=1,����ʵ����flag1=2,�����鲿
% flag2 ��ţ�flag2=1,���� G,DG, flag2=2 ����Gr,DGr,flag2=3,���� Gz,DGz
% Optional TBNM������Ż�����
% time ����ʱ��

if internal(2)<=1
    error('Ƶ�ʱ������ 1 ���ܽ��б䲽����');
    return;
end

tic;

% ������ֿ�
h1=10*h;
hx=internal(1):h:1;
hn=length(hx);
internal1=[internal(1),hx(hn)];
internal2=[hx(hn), internal(2)];

if r>0 && z<0
    switch flag1 % ʵ�����鲿
        case 1 % ʵ��
            switch flag2 % ���� G DG,Gr DGr,Gz DGz
                case 1 % G DG
                    % �����ֵ
                    RG=GreenSF(r,z,internal(1),flag1,flag2); % RG=[ԭʼֵ������ֵ] ���Ϊ����
                    [x1,yR1,OptR1]=GreenTBNM([RG(1),RG(2)], internal1, h, r, z, 1, 1);  % ����ʵ��
                    [x2,yR2,OptR2]=GreenTBNM([yR1(hn,1),yR1(hn,2)], internal2, h1, r, z, 1, 1);  % ����ʵ��
                    RIG=[yR1(1:hn-1,:);yR2];  % ���У���һ��Ϊԭʼֵ���ڶ���Ϊ����ֵ
                    Optional=[OptR1,OptR2];
                    xx=[x1(1:hn-1);x2];
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



