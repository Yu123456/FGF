function RG=GreenF(r,z,omega,flag)
% 本程序编写于2014年8月13日，硕士论文，参照手稿P87-89
% r,z 位置输入
% omega 频率（即初始频率值）
% flag 标号，flag=1,计算 G,DG, flag=2 计算G,DG,Gr,DGr,Gz,DGz
% RG=[G,DG,Gr,DGr,Gz,DGz] 输出为格林函数及其梯度方程值，及对 Omega 一阶导数值，第一行为
% 实部，第二行为虚部

G=[0;0];
DG=G;
Gr=G;
DGr=G;
Gz=G;
DGz=G;
if r>0 && z<0
    switch flag
        case 1  % 计算 G,DG
            k0=omega^2;
            X=k0*r;
            Y=-k0*z;
            Re_v=FDF(X,Y,flag);  % Re_v=[FXY, DxF, DyF, DxxF, DxyF, DyyF]
            FXY=Re_v(1);
            DxF=Re_v(2);
            DyF=Re_v(3);
            ekz=exp(-Y);
            bj0=besselj(0,X);
            bj1=besselj(1,X);
            G(1)=2/sqrt(r^2+z^2)+k0*FXY; % G 的实部
            G(2)=-2*pi*k0*ekz*bj0;  % G 的虚部
            DG(1)=2*omega*FXY+2*omega*k0*(r*DxF-z*DyF);  % DG 的实部
            DG(2)=-4*pi*omega*ekz*((1-Y)*bj0-X*bj1); % DG 的虚部
        case 2  % 计算G,DG,Gr,DGr,Gz,DGz
            k0=omega^2;
            k2=k0^2;
            X=k0*r;
            Y=-k0*z;
            rz=sqrt(r^2+z^2);
            Re_v=FDF(X,Y,flag);  % Re_v=[FXY, DxF, DyF, DxxF, DxyF, DyyF]
            FXY=Re_v(1);
            DxF=Re_v(2);
            DyF=Re_v(3);
            DxxF=Re_v(4);
            DxyF=Re_v(5);
            DyyF=Re_v(6);
            ekz=exp(-Y);
            bj0=besselj(0,X);
            bj1=besselj(1,X);
            G(1)=2/sqrt(r^2+z^2)+k0*FXY; % G 的实部
            G(2)=-2*pi*k0*ekz*bj0;  % G 的虚部
            DG(1)=2*omega*FXY+2*omega*k0*(r*DxF-z*DyF);  % DG 的实部
            DG(2)=-4*pi*omega*ekz*((1-Y)*bj0-X*bj1); % DG 的虚部
            Gr(1)=-2*r/rz^3+k2*DxF; % Gr 的实部
            Gr(2)=2*pi*k2*ekz*bj1;  % Gr 的虚部
            DGr(1)=4*k0*omega*DxF+2*k2*omega*(r*DxxF-z*DxyF); % DGr 的实部
            DGr(2)=4*pi*k0*omega*ekz*(X*bj0+(1-Y)*bj1);  % DGr 的虚部
            Gz(1)=-2*z/rz^3-k2*DyF; % Gz 的实部
            Gz(2)=-2*pi*k2*ekz*bj0; % Gz 的虚部
            DGz(1)=-4*k0*omega*DyF-2*k2*omega*(r*DxyF-z*DyyF); % DGz 的实部
            DGz(2)=-4*pi*k0*omega*ekz*((2-Y)*bj0-X*bj1); % DGz 的虚部
    end
else
    warning('r,z输入有误');
    return;
end

RG=[G,DG,Gr,DGr,Gz,DGz];

end

