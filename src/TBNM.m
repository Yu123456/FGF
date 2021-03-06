function [x, y]=TBNM(initialValue, internal, Omega, h)
% 本程序写于2014年8月6日，实现二阶线性微分方程求解
% 参照： Samuel N.Jator, S.Swindell, R.French. Trigonometrically fitted block
% Numerov type method for y''=f(x,y,y')[J]. Numer Algor. 2013,62:13-26.
% 手稿：2014年8月5日
% 微分方程类型 y''+g0(x)y'+g1(x)y=g2(x)

% 测试论文中 example 5.7 即 BESSELJ(1/2,x)

% initialValue: 初值 [y(x0), y'(x0)]
% internal: 求解区间 [x0, x1]
% Omega：基底三角函数参数
% h:步长

% y中第一列为y值，第二列为一阶导数值y'

x=internal(1):h:internal(2);
maxn=length(x);
% 转换成奇数个值
if mod(maxn,2)==0
    maxn=maxn-1;
end
x=x(1:maxn)';
% y中第一列为y值，第二列为一阶导数值y'
y=zeros(maxn,2);
y(1,:)=initialValue;

% 构造系数矩阵A
% 计算beta
u=Omega*h;
csc1=csc(u);
sec1=sec(u/2);
csc2=csc(u/2)^2;
u2=u^2;
u4=u^4;
u6=u^6;
u8=u^8;
if u<0.01
     beta0=1/12+u2/240+u4/6048+u6/172800+u8/5322240;
     beta1=5/6-u2/120-u4/3024-u6/86400-u8/2661120;
     beta2=beta0;
     beta00=-7/24-7*u2/480-71*u4/60480-53*u6/483840-23*u8/2128896;
     beta10=-1/4+u2/144+u4/4320+u6/134400+u8/4354560;
     beta20=1/24+11*u2/1440+19*u4/20160+247*u6/2419200+1013*u8/95800320;
     beta01=1/8+17*u2/1440+67*u4/60480+29*u6/268800+1031*u8/95800320;
     beta11=5/12-u2/240-u4/6048-u6/172800-u8/5322240;
     beta21=-1/24-11*u2/1440-19*u4/20160-247*u6/2419200-1013*u8/95800320;
     beta02=1/24-u2/288-47*u4/60480-233*u6/2419200-199*u8/19160064;
     beta12=13/12-11*u2/720-17*u4/30240-23*u6/1209600-29*u8/47900160;
     beta22=3/8+3*u2/160+3*u4/2240+31*u6/268800+13*u8/1182720;
else
     beta0=(-1+u2*csc2/4)/u2;
     beta1=(2+u2*(1-csc2/2))/u2;
     beta2=beta0;
     beta00=(-8-u*csc2*(u-2*sec1*sin(3*u/2)))/(8*u2);
     beta10=csc2*(2+(u2-2)*cos(u)-2*u*sin(u))/(4*u2);
     beta20=(-u*csc2+4*csc1)/(8*u);
     beta01=(-8+u*(u*csc2+4*csc1))/(8*u2);
     beta11=(2*u2+4-u2*csc2)/(4*u2);
     beta21=(u*csc2-4*csc1)/(8*u);
     beta02=(-8+u*(3*u*csc2-4*csc1))/(8*u2);
     beta12=(csc2*(2-(2+3*u2)*cos(u)+2*u*sin(u)))/(4*u2);
     beta22=(3*u*csc2-4*sin(u)*csc2+4*csc1)/(8*u);
end
A1=[beta1, beta1, beta2, beta2; beta10, beta10, beta20, beta20; beta11, beta11,...
    beta21, beta21; beta12, beta12, beta22, beta22];
A3=[-2, 0, 1, 0; -1, 0, 0, 0; -1, h, 0, 0; -1, 0, 0, h];
Y1=[beta0, beta1, beta1; beta00, beta10, beta20; beta01, beta11, beta21;...
    beta02, beta12, beta22];
Y3=[-1, 0; -1, -h; -1, 0; -1, 0];
Y4=[beta0, beta00, beta01, beta02]';
h2=h^2;

for i=1:2:maxn-1
    A2=[g1(x(i+1)), 0, 0, 0; 0, g0(x(i+1)), 0, 0; 0, 0, g1(x(i+2)), 0;...
        0, 0, 0, g0(x(i+2))];
    Y2=[g2(x(i)), g2(x(i+1)), g2(x(i+2))]';
    A=h2*A1*A2+A3;
    Y=h2*Y1*Y2+Y3*y(i,:)'-h2*(y(i,1)*g1(x(i))+y(i,2)*g0(x(i)))*Y4;
    X=A\Y;
    y(i+1:i+2,:)=[X(1), X(2); X(3), X(4)];
end
end

% 微分方程系数函数
function y=g0(x)

y=1./x;

end

function y=g1(x)

y=(x.^2-0.25)./x.^2;

end

function y=g2(x)

y=0;

end
