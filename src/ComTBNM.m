function ComTBNM(initialValue, internal, Omega, h)
% ������д��2014��8��6�գ���TBNM�������жԱȲ���
% ���������� example 5.7 �� BESSELJ(1/2,x)

% initialValue: ��ֵ [y(x0), y'(x0)]
% internal: ������� [x0, x1]
% Omega���������Ǻ�������
% h:����

[x, y]=TBNM(initialValue, internal, Omega, h);
yRK4=TwoRK4(initialValue, internal, h);
% ��ֵ
y_real=besselj(1/2,x);

% TBNM ���
error1=y(:,1)-y_real;
% RK4 ���
error2=yRK4(:,1)-y_real;

% ��ͼ
figure;
plot(x,y(:,1),'-r',x,yRK4(:,1),':b',x,y_real,'--k');
title(['TBNM--RK4--��ֵ�Ա�ͼ',' \omega=',num2str(Omega),',h=',num2str(h)]);
legend('TBNM','RK4','���');

figure;
plot(x,error1,'-r',x,error2,':k');
title(['TBNM--RK4--���ͼ',' \omega=',num2str(Omega),',h=',num2str(h)]);
legend('TBNM','RK4');
end

