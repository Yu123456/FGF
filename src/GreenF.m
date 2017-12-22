function RG=GreenF(r,z,omega,flag)
% �������д��2014��8��13�գ�˶ʿ���ģ������ָ�P87-89
% r,z λ������
% omega Ƶ�ʣ�����ʼƵ��ֵ��
% flag ��ţ�flag=1,���� G,DG, flag=2 ����G,DG,Gr,DGr,Gz,DGz
% RG=[G,DG,Gr,DGr,Gz,DGz] ���Ϊ���ֺ��������ݶȷ���ֵ������ Omega һ�׵���ֵ����һ��Ϊ
% ʵ�����ڶ���Ϊ�鲿

G=[0;0];
DG=G;
Gr=G;
DGr=G;
Gz=G;
DGz=G;
if r>0 && z<0
    switch flag
        case 1  % ���� G,DG
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
            G(1)=2/sqrt(r^2+z^2)+k0*FXY; % G ��ʵ��
            G(2)=-2*pi*k0*ekz*bj0;  % G ���鲿
            DG(1)=2*omega*FXY+2*omega*k0*(r*DxF-z*DyF);  % DG ��ʵ��
            DG(2)=-4*pi*omega*ekz*((1-Y)*bj0-X*bj1); % DG ���鲿
        case 2  % ����G,DG,Gr,DGr,Gz,DGz
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
            G(1)=2/sqrt(r^2+z^2)+k0*FXY; % G ��ʵ��
            G(2)=-2*pi*k0*ekz*bj0;  % G ���鲿
            DG(1)=2*omega*FXY+2*omega*k0*(r*DxF-z*DyF);  % DG ��ʵ��
            DG(2)=-4*pi*omega*ekz*((1-Y)*bj0-X*bj1); % DG ���鲿
            Gr(1)=-2*r/rz^3+k2*DxF; % Gr ��ʵ��
            Gr(2)=2*pi*k2*ekz*bj1;  % Gr ���鲿
            DGr(1)=4*k0*omega*DxF+2*k2*omega*(r*DxxF-z*DxyF); % DGr ��ʵ��
            DGr(2)=4*pi*k0*omega*ekz*(X*bj0+(1-Y)*bj1);  % DGr ���鲿
            Gz(1)=-2*z/rz^3-k2*DyF; % Gz ��ʵ��
            Gz(2)=-2*pi*k2*ekz*bj0; % Gz ���鲿
            DGz(1)=-4*k0*omega*DyF-2*k2*omega*(r*DxyF-z*DyyF); % DGz ��ʵ��
            DGz(2)=-4*pi*k0*omega*ekz*((2-Y)*bj0-X*bj1); % DGz ���鲿
    end
else
    warning('r,z��������');
    return;
end

RG=[G,DG,Gr,DGr,Gz,DGz];

end

