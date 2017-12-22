function y=TwoRK4(initialValue, internal, h)
% �������д��2014��8��6�գ�Runge-Kutta 4�׷�ֱ��������΢�ַ���
% ���գ� E. ��ķ��. ��΢���ֲ�[M]. p193-194.
% ΢�ַ������� y''=f(x,y,y')

% ���� BESSELJ(1/2,x)

% initialValue: ��ֵ [y(x0), y'(x0)]
% internal: ������� [x0, x1]
% h:����

% Ϊ��ƥ�� TBNM �������ԱȲ��ԣ��˴������� TBNM ��һ��
x=internal(1):h:internal(2);
maxn=length(x);
% ת����������ֵ
if mod(maxn,2)==0
    maxn=maxn-1;
end
x=x(1:maxn)';
% y�е�һ��Ϊyֵ���ڶ���Ϊһ�׵���ֵy'
y=zeros(maxn,2);
y(1,:)=initialValue;

for i=1:maxn-1
    y(i+1,:)=RK4(x(i),y(i,1),y(i,2),h);
end

end

function Re_y=RK4(x,y,dy,h)
% dy ����y��һ�׵���
% h ����

Re_y=zeros(1,2);
k1=h*f(x,y,dy);
k2=h*f(x+h/2, y+h*dy/2+h*k1/8, dy+k1/2);
k3=h*f(x+h/2, y+h*dy/2+h*k1/8, dy+k2/2);
k4=h*f(x+h, y+h*dy+h*k3/2, dy+k3);
A=dy+(k1+k2+k3)/6;
Re_y(1)=y+h*A;
Re_y(2)=A+(k2+k3+k4)/6;

end

% ���� f(x,y,y')
function value=f(x,y,dy)

value=-dy/x-(x^2-0.25)*y/x^2;

end