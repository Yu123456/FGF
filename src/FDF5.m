function Rev = FDF5(X,Y)
% ˶ʿ����������5�Ĳ��Գ��򣬲����Լ�д��˶ʿ���� 4.2.5 ��

FXY=0;
DxF=FXY;
DyF=FXY;
DxxF=FXY;
DxyF=FXY;
DyyF=FXY;
% ָ�����ּ�ָ���˻�
ey=exp(-Y)*ei(Y);
facn=1;  % n �Ľ׳˼���
X2=(X/2)^2; % X ����2��ƽ��
Xn=X2;
Ym=Y^2;
Y1=1/Y;
% n=1 ʱ
facm=1;  % ���� FXY, DxF
facm1=1;  % ���� DyF
summ=1/Y+1/Ym;   % ���� FXY,DxF
summ1=1/Ym+2/(Ym*Y);  % ���� DyF
flag=-1;
sumn=-Xn*(summ-ey);  % ���� FXY
sumn1=-Xn*(summ1+1/Y-ey); % ���� DyF
sumn2=-Xn*2*(summ-ey); % ���� DxF

for n=2:15
    flag=-flag;  %  ���� ��-1��^n
    facn=facn*n;
    Xn=Xn*X2;  % ���� (X/2)^(2n)
    %�ڲ���ͣ�ֻ������� m=2n-1,m=2n
    n2=2*n-2;
    n1=2*n-1;
    n3=2*n;
    facm=facm*n2;
    facm1=facm1*n1;
    Ym=Ym*Y;
    summ=summ+facm/Ym;
    summ1=summ1+facm1/Ym*Y1;  % �˴����һ�� Y,�ȵ��������ͺ�ͳһ����һ�� Y
    facm=facm*n1;
    facm1=facm1*n3;
    Ym=Ym*Y;
    summ=summ+facm/Ym;
    summ1=summ1+facm1/Ym*Y1;
    facn2=facn^2;
    sumn=sumn+flag*Xn*(summ-ey)/facn2;
    sumn1=sumn1+flag*Xn*(summ1+Y1-ey)/facn2;
    sumn2=sumn2+flag*n3*Xn*(summ-ey)/facn2;
end

FXY=-2*ey+2*sumn;
DxF=sumn2*2/X;
DyF=2*(ey-Y1)-2*sumn1;

Rev=[FXY, DxF, DyF];
end

