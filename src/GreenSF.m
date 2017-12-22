function RG=GreenSF(r,z,omega,flag1,flag2)
% 本程序编写于2014年8月28日，硕士论文，分别计算G,DG,Gr,DGr,Gz,DGz
% r,z 位置输入
% omega 频率（即初始频率值）
% flag1 标号，flag=1,计算实部，flag1=2,计算虚部
% flag2 标号，flag2=1,计算 G,DG, flag2=2 计算Gr,DGr,flag2=3,计算 Gz,DGz
% RG=[原始值，导数] 输出为格林函数及其梯度方程值，及对 Omega 一阶导数值，第一行为
% 实部，第二行为虚部


if r>0 && z<0
    switch flag1 % 实部或虚部
        case 1  % 实部
            switch flag2  % G DG Gr DGr Gz DGz
                case 1 % 计算 G DG
                    k0=omega^2;
                    X=k0*r;
                    Y=-k0*z;
                    Re_v=FDF(X,Y,1);  % Re_v=[FXY, DxF, DyF, DxxF, DxyF, DyyF]
                    FXY=Re_v(1);
                    DxF=Re_v(2);
                    DyF=Re_v(3);
                    GG=2/sqrt(r^2+z^2)+k0*FXY; % G 的实部
                    DGG=2*omega*FXY+2*omega*k0*(r*DxF-z*DyF);  % DG 的实部
                case 2 % 计算 Gr DGr
                    k0=omega^2;
                    k2=k0^2;
                    X=k0*r;
                    Y=-k0*z;
                    rz=sqrt(r^2+z^2);
                    Re_v=FDF(X,Y,2);  % Re_v=[FXY, DxF, DyF, DxxF, DxyF, DyyF]
                    DxF=Re_v(2);
                    DxxF=Re_v(4);
                    DxyF=Re_v(5);
                    GG=-2*r/rz^3+k2*DxF; % Gr 的实部
                    DGG=4*k0*omega*DxF+2*k2*omega*(r*DxxF-z*DxyF); % DGr 的实部
                case 3 % 计算 Gz DGz
                    k0=omega^2;
                    k2=k0^2;
                    X=k0*r;
                    Y=-k0*z;
                    rz=sqrt(r^2+z^2);
                    Re_v=FDF(X,Y,2);  % Re_v=[FXY, DxF, DyF, DxxF, DxyF, DyyF]
                    DyF=Re_v(3);
                    DxyF=Re_v(5);
                    DyyF=Re_v(6);
                    GG=-2*z/rz^3-k2*DyF; % Gz 的实部
                    DGG=-4*k0*omega*DyF-2*k2*omega*(r*DxyF-z*DyyF); % DGz 的实部
            end
        case 2  % 计算虚部
            switch flag2  % G DG Gr DGr Gz DGz
                case 1 % G DG
                    k0=omega^2;
                    X=k0*r;
                    Y=-k0*z;
                    ekz=exp(-Y);
                    bj0=besselj(0,X);
                    bj1=besselj(1,X);
                    GG=-2*pi*k0*ekz*bj0;  % G 的虚部
                    DGG=-4*pi*omega*ekz*((1-Y)*bj0-X*bj1); % DG 的虚部
                case 2 % Gr DGr
                    k0=omega^2;
                    k2=k0^2;
                    X=k0*r;
                    Y=-k0*z;
                    ekz=exp(-Y);
                    bj0=besselj(0,X);
                    bj1=besselj(1,X);
                    GG=2*pi*k2*ekz*bj1;  % Gr 的虚部
                    DGG=4*pi*k0*omega*ekz*(X*bj0+(1-Y)*bj1);  % DGr 的虚部
                case 3 % Gz DGz
                    k0=omega^2;
                    k2=k0^2;
                    X=k0*r;
                    Y=-k0*z;
                    ekz=exp(-Y);
                    bj0=besselj(0,X);
                    bj1=besselj(1,X);
                    GG=-2*pi*k2*ekz*bj0; % Gz 的虚部
                    DGG=-4*pi*k0*omega*ekz*((2-Y)*bj0-X*bj1); % DGz 的虚部
            end
    end
else
    warning('r,z输入有误');
    return;
end

RG=[GG,DGG];

end

