function [t y]=euler(int,y0,h)

t(1)=int(1);
y(1)=y0;
n=round((int(2)-int(1))/h);
for i=1:n
    t(i+1)=t(i)+h;
    y(i+1)=eulerstep(t(i),y(i),h);
end
plot(t,y);
end

function y=eulerstep(t,y,h)

y=y+h*ydot(t,y);

end

function z=ydot(t,y)
% 在此处输入方程形式

z=y;

end