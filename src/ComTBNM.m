function ComTBNM(initialValue, internal, Omega, h)
% 本程序写于2014年8月6日，对TBNM函数进行对比测试
% 测试论文中 example 5.7 即 BESSELJ(1/2,x)

% initialValue: 初值 [y(x0), y'(x0)]
% internal: 求解区间 [x0, x1]
% Omega：基底三角函数参数
% h:步长

[x, y]=TBNM(initialValue, internal, Omega, h);
yRK4=TwoRK4(initialValue, internal, h);
% 真值
y_real=besselj(1/2,x);

% TBNM 误差
error1=y(:,1)-y_real;
% RK4 误差
error2=yRK4(:,1)-y_real;

% 作图
figure;
plot(x,y(:,1),'-r',x,yRK4(:,1),':b',x,y_real,'--k');
title(['TBNM--RK4--真值对比图',' \omega=',num2str(Omega),',h=',num2str(h)]);
legend('TBNM','RK4','真解');

figure;
plot(x,error1,'-r',x,error2,':k');
title(['TBNM--RK4--误差图',' \omega=',num2str(Omega),',h=',num2str(h)]);
legend('TBNM','RK4');
end

