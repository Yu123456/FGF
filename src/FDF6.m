function Rev = FDF6(X,Y)
% ˶ʿ����������5�Ĳ��Գ��򣬲����Լ�д��˶ʿ���� 4.2.5 ��
% ���� FXY, DxF, DyF, DxxF, DxyF, DyyF


% ָ�����ּ�ָ���˻�
ey=exp(-Y)*ei(Y);
facn=1;  % n �Ľ׳˼���
X2=(X/2)^2; % X ����2��ƽ��
Xn=X2;
Ym=Y^2;
Y1=1/Y;
Y2=Y1/Y;
Y21=Y2+Y1-ey;
% n=1 ʱ
facm=1;  % ���� FXY, DxF, DxxF
facm1=1;  % ���� DyF, DxyF
facm2=8; % ���� DyyF
summ=1/Y+1/Ym;   % ���� FXY,DxF,DxxF
summ1=1/Ym+2/Ym*Y1;  % ���� DyF, DxyF
summ2=2/Ym*Y1+6/Ym*Y2; % ���� DyyF
flag5=-1;
sumn=-Xn*(summ-ey);  % ���� FXY
sumn1=-Xn*(summ1+1/Y-ey); % ���� DyF
sumn2=-Xn*2*(summ-ey); % ���� DxF
sumn3=-Xn*2*(summ1+Y1-ey); % ���� DxyF
sumn4=-(summ-ey); % ���� DxxF
sumn5=-Xn*(summ2+Y2+Y1-ey); % ���� DyyF

for n=2:15
    flag5=-flag5;  %  ���� ��-1��^n
    facn=facn*n;
    Xn=Xn*X2;  % ���� (X/2)^(2n)
    %�ڲ���ͣ�ֻ������� m=2n-1,m=2n
    n2=2*n-2;
    n1=2*n-1;
    n3=2*n;
    facm=facm*n2;
    facm1=facm1*n1;
    facm2=facm2*n3;
    Ym=Ym*Y;
    summ=summ+facm/Ym;
    summ1=summ1+facm1/Ym*Y1; 
    summ2=summ2+facm2/Ym*Y2;
    facm=facm*n1;
    facm1=facm1*n3;
    facm2=facm2*(n3+1);
    Ym=Ym*Y;
    summ=summ+facm/Ym;
    summ1=summ1+facm1/Ym*Y1;
    summ2=summ2+facm2/Ym*Y2;
    facn2=facn^2;
    sumn=sumn+flag5*Xn*(summ-ey)/facn2;
    sumn1=sumn1+flag5*Xn*(summ1+Y1-ey)/facn2;
    sumn2=sumn2+flag5*n3*Xn*(summ-ey)/facn2;
    sumn3=sumn3+flag5*n3*Xn*(summ1+Y1-ey)/facn2;
    sumn4=sumn4+flag5*n*n1*Xn*(summ-ey)/facn2;
    sumn5=sumn5+flag5*Xn*(summ2+Y21)/facn2;
end

FXY=-2*ey+2*sumn;
DxF=sumn2*2/X;
DyF=2*(ey-Y1)-2*sumn1;
DxyF=-sumn3*2/X;
DxxF=sumn4/X2;
DyyF=2*Y21+2*sumn5;

Rev=[FXY, DxF, DyF, DxxF, DxyF, DyyF];
end

