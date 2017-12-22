function y=TwoRK4(initialValue, internal, h)
% 本程序编写于2014年8月6日，Runge-Kutta 4阶法直接求解二阶微分方程
% 参照： E. 卡姆克. 常微分手册[M]. p193-194.
% 微分方程类型 y''=f(x,y,y')

% 测试 BESSELJ(1/2,x)

% initialValue: 初值 [y(x0), y'(x0)]
% internal: 求解区间 [x0, x1]
% h:步长

% 为了匹配 TBNM 函数做对比测试，此处点数与 TBNM 中一致
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

for i=1:maxn-1
    y(i+1,:)=RK4(x(i),y(i,1),y(i,2),h);
end

end

function Re_y=RK4(x,y,dy,h)
% dy 代表y的一阶导数
% h 步长

Re_y=zeros(1,2);
k1=h*f(x,y,dy);
k2=h*f(x+h/2, y+h*dy/2+h*k1/8, dy+k1/2);
k3=h*f(x+h/2, y+h*dy/2+h*k1/8, dy+k2/2);
k4=h*f(x+h, y+h*dy+h*k3/2, dy+k3);
A=dy+(k1+k2+k3)/6;
Re_y(1)=y+h*A;
Re_y(2)=A+(k2+k3+k4)/6;

end

% 函数 f(x,y,y')
function value=f(x,y,dy)

value=-dy/x-(x^2-0.25)*y/x^2;

end