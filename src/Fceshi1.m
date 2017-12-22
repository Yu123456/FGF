function Re_v=Fceshi1(X,Y)
% 本程序编写于2014年8月27日，测试F(X,Y)计算 0<X<8,Y>20


if 0<X && X<8 && Y>=20  % 区域 Sigma4
    ey=exp(-Y);
    In1=(1/Y-2/Y^3)*ey-2*(1/Y^2-1/Y^3);  % 递推初值, I_1,后续代表 I_{n-1}
    In2=(-24/Y^5+4/Y^3-1/Y)*ey+(24/Y^5-24/Y^4+8/Y^3);  % 递推初值，I_2,后续代表 I_{n-2}
    R=sqrt(X^2+Y^2);
    YR=Y/R;
    YR2=YR^2;
    pm=1;   % -1 的n次幂
    % 计算 FXY
    su1=0;
    % 计算 DxF
    su2=0;
    % 计算 n=1
    pm=-pm;
    YRn=YR*YR2;
    factor=1/2; % 计算 FXY 类阶乘
    su1=su1+pm*factor*YRn*In1;
    factor1=3/2; % 计算 DxF 类阶乘
    su2=su2+pm*factor1*YRn*In1;
    % 计算 n=2
    pm=-pm;
    YRn=YR*YR2;
    factor=factor*3/2;
    su1=su1+pm*factor*YRn*In2;
    factor1=factor1*5/2;
    su2=su2+pm*factor1*YRn*In2;
    Y2=Y^2;
    for i=3:10
        pm=-pm;
        n2=2*i;
        % 类阶乘
        factor=factor*(1-1/n2);
        factor1=factor1*(1+1/n2);
        YRn=YRn*YR2;
        % 计算递推式
        In=-pm*ey/Y+n2*(n2-1)*In1/Y2+2*n2*(i-1)*In2/Y2;
        su1=su1+pm*factor*YRn*In;
        su2=su2+pm*factor1*YRn*In;
        % 递推式指标增加
        In2=In1;
        In1=In;
    end
    H0=StruveH0(X);
    H1=StruveH1(X);
    by0=bessely(0,X);
    by1=bessely(1,X);
    FXY=-pi*ey*(H0+by0)-2*(1-ey)/R-2*su1;
    DxF=-pi*ey*(2/pi-H1-by1)+2*X*(1-ey)/R^3+2*X*su2*YR2/Y2;
else
    str='X,Y 不在此区间';
    warning(str);
    return;
end

Re_v=[FXY, DxF];

end

