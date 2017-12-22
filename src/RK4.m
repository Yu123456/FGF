function [t,y]= RK4(int,y0,h)
% Input: interval [a,b], initial value y0, step size h
% Output: time steps t, solution y
% int is interval [a,b]
% 4阶Runge-Kutta法

t(1)=int(1);
y(1)=y0;
n=round((int(2)-int(1))/h);
for i=1:n
   t(i+1)=t(i)+h;
   y(i+1)=RK4step(t(i),y(i),h);
end
plot(t,y);
end

function ry=RK4step(t,y,h)
% Runge-Kutta Method
% Input: current time t, current value y, stepsize h
% Output: approximate solution value at time t+h
s1=ydot(t,y);
s2=ydot(t+h/2,y+h/2*s1);
s3=ydot(t+h/2,y+h/2*s2);
s4=ydot(t+h,y+h*s3);
ry=y+h/6*(s1+2*s2+2*s3+s4);
end

function z=ydot(t,y)
% 函数变量
z=t*y+t^3;
end

