function Rev = FDF5(X,Y)
% 硕士论文中区域5的测试程序，参照自己写的硕士论文 4.2.5 节

FXY=0;
DxF=FXY;
DyF=FXY;
DxxF=FXY;
DxyF=FXY;
DyyF=FXY;
% 指数积分及指数乘积
ey=exp(-Y)*ei(Y);
facn=1;  % n 的阶乘计算
X2=(X/2)^2; % X 除以2的平方
Xn=X2;
Ym=Y^2;
Y1=1/Y;
% n=1 时
facm=1;  % 计算 FXY, DxF
facm1=1;  % 计算 DyF
summ=1/Y+1/Ym;   % 计算 FXY,DxF
summ1=1/Ym+2/(Ym*Y);  % 计算 DyF
flag=-1;
sumn=-Xn*(summ-ey);  % 计算 FXY
sumn1=-Xn*(summ1+1/Y-ey); % 计算 DyF
sumn2=-Xn*2*(summ-ey); % 计算 DxF

for n=2:15
    flag=-flag;  %  计算 （-1）^n
    facn=facn*n;
    Xn=Xn*X2;  % 计算 (X/2)^(2n)
    %内层求和，只需加两项 m=2n-1,m=2n
    n2=2*n-2;
    n1=2*n-1;
    n3=2*n;
    facm=facm*n2;
    facm1=facm1*n1;
    Ym=Ym*Y;
    summ=summ+facm/Ym;
    summ1=summ1+facm1/Ym*Y1;  % 此处提出一个 Y,等到最后求完和后统一除以一个 Y
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

