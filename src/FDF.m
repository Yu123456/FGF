function Re_v=FDF(X,Y,flag)
% 本程序编写于2014年8月8日，硕士毕业论文计算 F(X,Y)
% 参照手稿 p85-98,p111-113
% X,Y 输入
% Re_v 返回 FXY, DxF, DxxF, DyF, DyyF, DxyF，分别表示函数F(X,Y)函数值，对x一阶偏导数，
% 对x二阶偏导数，对y一阶偏导数，对y二阶偏导数，对x，y混合偏导数
% flag 标号   flag=1 计算 FXY，DxF,DyF,对应G, flag=2 计算FXY,DxF,DyF,DxxF,DxyF,
% DyyF,对应Gr,Gz，当然也包含G

FXY=0;
DxF=FXY;
DyF=FXY;
DxxF=FXY;
DxyF=FXY;
DyyF=FXY;


if X<=0 || Y<=0 
    warning('X,Y 输入有误！');
    return;
end
if flag==1 || flag==2 
    switch flag
        case 1 %计算FXY，DxF,DyF
            if (X>0 && X<=8 && Y>0 && Y<=18)  % Sigma1，对应 D1
                f0=fn(X,Y,0);
                % 计算 FXY,DyF
                su1=f0;
                % 计算 DxF
                su2=f0;
                factor=1;
                for i=1:40
                    % 阶乘
                    factor=factor*i;
                    fv=fn(X,Y,i)/factor;
%                     disp([i,vpa(fv,10)]);
                    su1=su1+(1+i)*fv;
                    su2=su2+fv;
                end
                R=sqrt(X^2+Y^2);
                ey=exp(-Y);
                piey=pi*ey*bessely(0,X);
                FXY=-piey-2*R/X^2+2*ey*R*su1/X^2;
                DxF=pi*ey*bessely(1,X)+2*Y/(X*R)-2*ey*R*su2/X;
                DyF=(2*Y^2)/(X^2*R)+piey-2*R*ey*su1/X^2;   
%             elseif X>0 && X<=8 && Y>18 && Y<20 % Sigma2   划归 D5
%                 Y0=18;
%                 R0=sqrt(X^2+Y0^2);
%                 R=sqrt(X^2+Y^2);
%                 f0=fn(X,Y0,0);
%                 % 计算 FXY,DyF
%                 su1=f0;
%                 GF1=Gauss16_1(R0,R,X);  % 计算高斯积分
%                 % 计算 DxF
%                 su2=f0;
%                 GF2=Gauss16_2(R0,R,X);  % 计算高斯积分
%                 factor=1;
%                 for i=1:40
%                     % 阶乘
%                     factor=factor*i;
%                     fv=fn(X,Y0,i)/factor;
%                     su1=su1+(1+i)*fv;
%                     su2=su2+fv;
%                 end
%                 by0=bessely(0,X);  % bessel Y0(X)
%                 by1=bessely(1,X);  % bessel Y1(X)
%                 ey=exp(-Y);
%                 ey18=exp(18-Y);
%                 FXY=-pi*ey*by0-2*R0*ey18/X^2-2*ey*GF1+2*R0*ey*su1/X^2;
%                 DxF=pi*ey*by1+2*Y/(X*R)-2*ey*GF2/X-2*R0*ey*su2/X;
%                 DyF=-2/R+pi*ey*by0+2*R0*ey18/X^2+2*ey*GF1-2*R0*ey*su1/X^2;
            elseif X>8 && Y<=20 && X/Y<2 % 区域 Sigma6，对应 D2
                %disp('Sigma6');
                ey=exp(-Y);
                In1=(1/Y-2/Y^3)*ey-2*(1/Y^2-1/Y^3);  % 递推初值, I_1,后续代表 I_{n-1}
                In2=(1-ey)/Y;  % 递推初值，I_0,后续代表 I_{n-2}
                R=sqrt(X^2+Y^2);
                YR=Y/R;
                YR2=YR^2;
                pm=1;   % -1 的n次幂
                % 计算 FXY,DyF
                su1=0;
                % 计算 DxF
                su2=0;
                % 计算 n=1
                pm=-pm;
                YRn=YR*YR2;
                factor=1/2; % 计算 FXY,DyF类阶乘
                su1=su1+pm*factor*YRn*In1;
                factor1=3/2; % 计算 DxF 类阶乘
                su2=su2+pm*factor1*YRn*In1;
                Y2=Y^2;
                for i=2:20
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
                DyF=-2*ey/R+pi*ey*(H0+by0)+2*su1;
           elseif X>8 && X/Y>=2 && Y<=20  % 区域 Sigma5,对应 D3
                ey=exp(-Y);
                I2n=1-ey;  % 递推初值
                X2=X^2;
                Xn=1/X;   % X 负幂次初值
                % 计算 FXY,DyF
                su1=0;
                % 计算 DxF
                su2=0;
                factor=1; % 阶乘
                factor1=1; % 奇数阶乘除以2的初值
                pm=1;   % -1 的n次幂
                % 项数截断
                if X/Y>4
                    N=7;
                else
                    N=15;
                end
                for i=1:N
                    pm=-pm;
                    % 阶乘
                    factor=factor*i;
                    factor1=factor1*(2*i-1)/2;
                    Xn=Xn/X2;
                    I2n=I2nF(I2n,Y,i);
                    su1=su1+pm*factor1*I2n*Xn/factor;
                    su2=su2+pm*(2*i+1)*factor1*I2n*Xn/(X*factor);
                end
                H0=StruveH0(X);
                H1=StruveH1(X);
                by0=bessely(0,X);
                by1=bessely(1,X);
                R=sqrt(X^2+Y^2);
                FXY=-pi*ey*(H0+by0)-2*(1-ey)/X-2*su1;
                DxF=-pi*ey*(2/pi-H1-by1)+2*(1-ey)/X2+2*su2;
                DyF=-2/R+pi*ey*(H0+by0)+2*(1-ey)/X+2*su1;
            elseif X>0 && X<8 && Y>=18 % 对应 D5
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
                flag5=-1;
                sumn=-Xn*(summ-ey);  % 计算 FXY
                sumn1=-Xn*(summ1+1/Y-ey); % 计算 DyF
                sumn2=-Xn*2*(summ-ey); % 计算 DxF

                for n=2:15
                    flag5=-flag5;  %  计算 （-1）^n
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
                    summ1=summ1+facm1/Ym*Y1;  
                    facm=facm*n1;
                    facm1=facm1*n3;
                    Ym=Ym*Y;
                    summ=summ+facm/Ym;
                    summ1=summ1+facm1/Ym*Y1;
                    facn2=facn^2;
                    sumn=sumn+flag5*Xn*(summ-ey)/facn2;
                    sumn1=sumn1+flag5*Xn*(summ1+Y1-ey)/facn2;
                    sumn2=sumn2+flag5*n3*Xn*(summ-ey)/facn2;
                end

                FXY=-2*ey+2*sumn;
                DxF=sumn2*2/X;
                DyF=2*(ey-Y1)-2*sumn1;
            else  % Sigma3,对应 D4
                R=sqrt(X^2+Y^2);
                stheta=Y/R;
                % 计算 FXY,DyF
                Pm1=Pmn(stheta,0,0);
                factor1=1;
                su1=factor1*Pm1/R;
                % 计算 DxF
                Pm2=Pmn(stheta,1,-1);
                factor2=2;
                su2=factor2*Pm2/R^2;
                for i=1:5
                    factor1=factor1*i;
                    Pm1=Pmn(stheta,i,0);
                    su1=su1+factor1*Pm1/R^(i+1);
                    factor2=factor2*(i+2);
                    Pm2=Pmn(stheta,i+1,-1);
                    su2=su2+factor2*Pm2/R^(i+2);
                end
                H0=StruveH0(X);
                H1=StruveH1(X);
                by0=bessely(0,X);
                by1=bessely(1,X);
                ey=exp(-Y);
                FXY=-pi*ey*(H0+by0)-2*su1;
                DxF=-pi*ey*(2/pi-H1-by1)+2*su2;
                DyF=-2/R+pi*ey*(H0+by0)+2*su1;
            end
        case 2 %计算FXY,DxF,DyF,DxxF,DxyF,DyyF
             if X>0 && X<=8 && Y>0 && Y<=18  % 对应 D1
                f0=fn(X,Y,0);
                % 计算 FXY,DyF,DyyF
                su1=f0;
                % 计算 DxF,DxyF
                su2=f0;
                % 计算 DxxF
                su3=0;
                factor=1;
                for i=1:40
                    % 阶乘
                    factor=factor*i;
                    fv=fn(X,Y,i)/factor;
%                     disp([i,vpa(fv,10)]);
                    su1=su1+(1+i)*fv;
                    su2=su2+fv;
                    su3=su3+i*fv;
                end
                R=sqrt(X^2+Y^2);
                ey=exp(-Y);
                by0=bessely(0,X);
                by1=bessely(1,X);
                by2=bessely(2,X);
                piey=pi*ey*by0;
                FXY=-piey-2*R/X^2+2*ey*R*su1/X^2;
                DxF=pi*ey*by1+2*Y/(X*R)-2*ey*R*su2/X;
                DyF=(2*Y^2)/(X^2*R)+piey-2*R*ey*su1/X^2;   
                DxxF=pi*ey*(by0-by2)/2-2*Y*(2*X^2+Y^2)/(X^2*R^3)+2*Y^2/(X^2*R)-2*ey*R*su3/X^2;
                DxyF=-pi*ey*by1+2*X/R^3-2*Y/(X*R)+2*R*ey*su2/X;
                DyyF=2*Y/R^3-2*Y^2/(X^2*R)-pi*ey*by0+2*R*ey*su1/X^2;
%             elseif X>0 && X<=8 && Y>18 && Y<=20  % 该区域划归 D5
%                 disp('running vD');
%                 Y0=18;
%                 R0=sqrt(X^2+Y0^2);
%                 R=sqrt(X^2+Y^2);
%                 f0=fn(X,Y0,0);
%                 % 计算 FXY,DyF,DyyF
%                 su1=f0;
%                 GF1=Gauss16_1(R0,R,X);  % 计算高斯积分
%                 % 计算 DxF,DxyF
%                 su2=f0;
%                 GF2=Gauss16_2(R0,R,X);  % 计算高斯积分
%                 % 计算 DxxF
%                 su3=0;
%                 factor=1;
%                 for i=1:40
%                     % 阶乘
%                     factor=factor*i;
%                     fv=fn(X,Y0,i)/factor;
%                     su1=su1+(1+i)*fv;
%                     su2=su2+fv;
%                     su3=su3+i*fv;
%                 end
%                 by0=bessely(0,X);  % bessel Y0(X)
%                 by1=bessely(1,X);  % bessel Y1(X)
%                 by2=bessely(2,X);  % bessel Y2(X)
%                 ey=exp(-Y);
%                 ey18=exp(18-Y);
%                 FXY=-pi*ey*by0-2*R0*ey18/X^2-2*ey*GF1+2*R0*ey*su1/X^2;
%                 DxF=pi*ey*by1+2*Y/(X*R)-2*ey*GF2/X-2*R0*ey*su2/X;
%                 DyF=-2/R+pi*ey*by0+2*R0*ey18/X^2+2*ey*GF1-2*R0*ey*su1/X^2;
%                 DxxF=pi*ey*(by0-by2)/2-2*Y*(2*X^2+Y^2)/(X^2*R^3)+2*ey*GF2/X^2-...
%                     2/R+2*R0*ey18/X^2+2*ey*GF1-2*R0*ey*su3/X^2;
%                 DxyF=-pi*ey*by1+2*X/R^3-2*Y/(X*R)+2*ey*GF2/X+2*R0*ey*su2/X;
%                 DyyF=2*(Y+R^2)/R^3-pi*ey*by0-2*R0*ey18/X^2-2*ey*GF1+2*R0*ey*su1/X^2;
            elseif X>8 && Y<=20 && X/Y<2 % 区域 Sigma6,对应 D2
                ey=exp(-Y);
                In1=(1/Y-2/Y^3)*ey-2*(1/Y^2-1/Y^3);  % 递推初值, I_1,后续代表 I_{n-1}
                In2=(1-ey)/Y;  % 递推初值，I_0,后续代表 I_{n-2}
                R=sqrt(X^2+Y^2);
                YR=Y/R;
                YR2=YR^2;
                Y2=Y^2;
                X2R2=X^2/R^2;
                R3=R^3;
                pm=1;   % -1 的n次幂
                % 计算 FXY,DyF,DyyF
                su1=0;
                % 计算 DxF,DxyF
                su2=0;
                % 计算 DxxF
                su3=0;   % 只需要在 su2 基础上乘以 (1-X^2(2*n+3)/R^2)
                % 计算 n=1
                pm=-pm;
                YRn=YR*YR2;
                factor=1/2; % 计算 FXY,DyF类阶乘
                su1=su1+pm*factor*YRn*In1;
                factor1=3/2; % 计算 DxF 类阶乘
                su23=pm*factor1*YRn*In1;  % su2,su3 共同项
                su2=su2+su23;
                su3=su3+su23*(1-5*X2R2);
                for i=2:16
                    pm=-pm;
                    n2=2*i;
                    % 类阶乘
                    factor=factor*(1-1/n2);
                    factor1=factor1*(1+1/n2);
                    YRn=YRn*YR2;
                    % 计算递推式
                    In=-pm*ey/Y+n2*(n2-1)*In1/Y2+2*n2*(i-1)*In2/Y2;
                    su1=su1+pm*factor*YRn*In;
                    su23=pm*factor1*YRn*In;  % su2,su3 共同项
                    su2=su2+su23;
                    su3=su3+su23*(1-X2R2*(n2+3));
                    % 递推式指标增加
                    In2=In1;
                    In1=In;
                end
                H0=StruveH0(X);
                H1=StruveH1(X);
                by0=bessely(0,X);
                by1=bessely(1,X);
                by2=bessely(2,X);
                FXY=-pi*ey*(H0+by0)-2*(1-ey)/R-2*su1;
                DxF=-pi*ey*(2/pi-H1-by1)+2*X*(1-ey)/R3+2*X*su2*YR2/Y2;
                DyF=-2*ey/R+pi*ey*(H0+by0)+2*su1;
                DxxF=pi*ey*(2*H0+by0-2*H1/X-by2)/2+(2/R3-6*X^2/R^5)*(1-ey)+2*su3*YR2/Y2;
                DxyF=(2*X/R3+2-pi*(H1+by1))*ey-2*X*su2*YR2/Y2;
                DyyF=2*Y/R3+2*ey/R-pi*ey*(H0+by0)-2*su1;
            elseif X>8 && X/Y>=2 && Y<=20  % 区域 Sigma5，对应 D3
                ey=exp(-Y);
                I2n=1-ey;  % 递推初值
                X2=X^2;
                Xn=1/X;   % X 负幂次初值
                % 计算 FXY,DyF,DyyF
                su1=0;
                % 计算 DxF,DxyF
                su2=0;
                % 计算 DxxF
                su3=0;
                factor=1; % 阶乘
                factor1=1; % 奇数阶乘除以2的初值
                pm=1;   % -1 的n次幂
                % 项数截断
                if X/Y>4
                    N=7;
                else
                    N=15;
                end
                for i=1:N
                    pm=-pm;
                    i2=2*i;
                    % 阶乘
                    factor=factor*i;
                    factor1=factor1*(i2-1)/2;
                    Xn=Xn/X2;
                    I2n=I2nF(I2n,Y,i);
                    su1=su1+pm*factor1*I2n*Xn/factor;
                    su2=su2+pm*(i2+1)*factor1*I2n*Xn/(X*factor);
                    su3=su3+pm*(i2+1)*(i2+2)*factor1*I2n*Xn/(factor);
                end
                H0=StruveH0(X);
                H1=StruveH1(X);
                by0=bessely(0,X);
                by1=bessely(1,X);
                by2=bessely(2,X);
                R=sqrt(X^2+Y^2);
                FXY=-pi*ey*(H0+by0)-2*(1-ey)/X-2*su1;
                DxF=-pi*ey*(2/pi-H1-by1)+2*(1-ey)/X2+2*su2;
                DyF=-2/R+pi*ey*(H0+by0)+2*(1-ey)/X+2*su1;
                DxxF=pi*ey*(2*H0+by0-2*H1/X-by2)/2-4*(1-ey)/X^3-2*su3/X2;
                DxyF=2*X/R^3+pi*ey*(2/pi-H1-by1)-2*(1-ey)/X2-2*su2;
                DyyF=2*(Y+R^2)/R^3-2*(1-ey)/X-pi*ey*(H0+by0)-2*su1;
             elseif X<8 && Y>=18  % 对应 D5
                 % 指数积分及指数乘积
                ey=exp(-Y)*ei(Y);
                facn=1;  % n 的阶乘计算
                X2=(X/2)^2; % X 除以2的平方
                Xn=X2;
                Ym=Y^2;
                Y1=1/Y;
                Y2=Y1/Y;
                Y21=Y2+Y1-ey;
                % n=1 时
                facm=1;  % 计算 FXY, DxF, DxxF
                facm1=1;  % 计算 DyF, DxyF
                facm2=8; % 计算 DyyF
                summ=1/Y+1/Ym;   % 计算 FXY,DxF,DxxF
                summ1=1/Ym+2/Ym*Y1;  % 计算 DyF, DxyF
                summ2=2/Ym*Y1+6/Ym*Y2; % 计算 DyyF
                flag5=-1;
                sumn=-Xn*(summ-ey);  % 计算 FXY
                sumn1=-Xn*(summ1+1/Y-ey); % 计算 DyF
                sumn2=-Xn*2*(summ-ey); % 计算 DxF
                sumn3=-Xn*2*(summ1+Y1-ey); % 计算 DxyF
                sumn4=-(summ-ey); % 计算 DxxF
                sumn5=-Xn*(summ2+Y2+Y1-ey); % 计算 DyyF

                for n=2:15
                    flag5=-flag5;  %  计算 （-1）^n
                    facn=facn*n;
                    Xn=Xn*X2;  % 计算 (X/2)^(2n)
                    %内层求和，只需加两项 m=2n-1,m=2n
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
             else  % 对应 D4
                % disp('running 4')
                R=sqrt(X^2+Y^2);
                stheta=Y/R;
                % 计算 FXY,DyF,DyyF
                Pm1=Pmn(stheta,0,0);
                factor1=1;
                su1=factor1*Pm1/R;
                % 计算 DxF,DxyF
                Pm2=Pmn(stheta,1,-1);
                factor2=2;
                su2=factor2*Pm2/R^2;
                % 计算 DxxF
                Pm3=Pmn(stheta,2,-2);
                Pm4=Pmn(stheta,2,0);
                factor3=factorial(4);
                su3=(factor3*Pm3-factor2*Pm4)/R^3;
                for i=1:5
                    factor1=factor1*i;
                    Pm1=Pmn(stheta,i,0);
                    su1=su1+factor1*Pm1/R^(i+1);
                    factor2=factor2*(i+2);
                    Pm2=Pmn(stheta,i+1,-1);
                    su2=su2+factor2*Pm2/R^(i+2);
                    factor3=factor3*(i+4);
                    Pm3=Pmn(stheta,i+2,-2);
                    Pm4=Pmn(stheta,i+2,0);
                    su3=su3+(factor3*Pm3-factor2*Pm4)/R^(i+3);
                end
                H0=StruveH0(X);
                H1=StruveH1(X);
                by0=bessely(0,X);
                by1=bessely(1,X);
                by2=bessely(2,X);
                ey=exp(-Y);
                FXY=-pi*ey*(H0+by0)-2*su1;
                DxF=-pi*ey*(2/pi-H1-by1)+2*su2;
                DyF=-2/R+pi*ey*(H0+by0)+2*su1;
                DxxF=pi*ey*(2*H0+by0-2*H1/X-by2)/2-su3;
                DxyF=2*X/R^3+2*ey-pi*ey*(H1+by1)-2*su2;
                DyyF=2*Y/R^3+2/R-pi*ey*(H0+by0)-2*su1;
            end
    end
else
    warning('flag 输入错误（flag=1,2,3）！');
    return;
end
Re_v=[FXY, DxF, DyF, DxxF, DxyF, DyyF];

end

