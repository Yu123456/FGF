function RG=GreenSF(r,z,omega,flag1,flag2)
% �������д��2014��8��28�գ�˶ʿ���ģ��ֱ����G,DG,Gr,DGr,Gz,DGz
% r,z λ������
% omega Ƶ�ʣ�����ʼƵ��ֵ��
% flag1 ��ţ�flag=1,����ʵ����flag1=2,�����鲿
% flag2 ��ţ�flag2=1,���� G,DG, flag2=2 ����Gr,DGr,flag2=3,���� Gz,DGz
% RG=[ԭʼֵ������] ���Ϊ���ֺ��������ݶȷ���ֵ������ Omega һ�׵���ֵ����һ��Ϊ
% ʵ�����ڶ���Ϊ�鲿


if r>0 && z<0
    switch flag1 % ʵ�����鲿
        case 1  % ʵ��
            switch flag2  % G DG Gr DGr Gz DGz
                case 1 % ���� G DG
                    k0=omega^2;
                    X=k0*r;
                    Y=-k0*z;
                    Re_v=FDF(X,Y,1);  % Re_v=[FXY, DxF, DyF, DxxF, DxyF, DyyF]
                    FXY=Re_v(1);
                    DxF=Re_v(2);
                    DyF=Re_v(3);
                    GG=2/sqrt(r^2+z^2)+k0*FXY; % G ��ʵ��
                    DGG=2*omega*FXY+2*omega*k0*(r*DxF-z*DyF);  % DG ��ʵ��
                case 2 % ���� Gr DGr
                    k0=omega^2;
                    k2=k0^2;
                    X=k0*r;
                    Y=-k0*z;
                    rz=sqrt(r^2+z^2);
                    Re_v=FDF(X,Y,2);  % Re_v=[FXY, DxF, DyF, DxxF, DxyF, DyyF]
                    DxF=Re_v(2);
                    DxxF=Re_v(4);
                    DxyF=Re_v(5);
                    GG=-2*r/rz^3+k2*DxF; % Gr ��ʵ��
                    DGG=4*k0*omega*DxF+2*k2*omega*(r*DxxF-z*DxyF); % DGr ��ʵ��
                case 3 % ���� Gz DGz
                    k0=omega^2;
                    k2=k0^2;
                    X=k0*r;
                    Y=-k0*z;
                    rz=sqrt(r^2+z^2);
                    Re_v=FDF(X,Y,2);  % Re_v=[FXY, DxF, DyF, DxxF, DxyF, DyyF]
                    DyF=Re_v(3);
                    DxyF=Re_v(5);
                    DyyF=Re_v(6);
                    GG=-2*z/rz^3-k2*DyF; % Gz ��ʵ��
                    DGG=-4*k0*omega*DyF-2*k2*omega*(r*DxyF-z*DyyF); % DGz ��ʵ��
            end
        case 2  % �����鲿
            switch flag2  % G DG Gr DGr Gz DGz
                case 1 % G DG
                    k0=omega^2;
                    X=k0*r;
                    Y=-k0*z;
                    ekz=exp(-Y);
                    bj0=besselj(0,X);
                    bj1=besselj(1,X);
                    GG=-2*pi*k0*ekz*bj0;  % G ���鲿
                    DGG=-4*pi*omega*ekz*((1-Y)*bj0-X*bj1); % DG ���鲿
                case 2 % Gr DGr
                    k0=omega^2;
                    k2=k0^2;
                    X=k0*r;
                    Y=-k0*z;
                    ekz=exp(-Y);
                    bj0=besselj(0,X);
                    bj1=besselj(1,X);
                    GG=2*pi*k2*ekz*bj1;  % Gr ���鲿
                    DGG=4*pi*k0*omega*ekz*(X*bj0+(1-Y)*bj1);  % DGr ���鲿
                case 3 % Gz DGz
                    k0=omega^2;
                    k2=k0^2;
                    X=k0*r;
                    Y=-k0*z;
                    ekz=exp(-Y);
                    bj0=besselj(0,X);
                    bj1=besselj(1,X);
                    GG=-2*pi*k2*ekz*bj0; % Gz ���鲿
                    DGG=-4*pi*k0*omega*ekz*((2-Y)*bj0-X*bj1); % DGz ���鲿
            end
    end
else
    warning('r,z��������');
    return;
end

RG=[GG,DGG];

end

