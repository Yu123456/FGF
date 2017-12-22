function ceshi()

x=0:0.01:pi;
y1=sin(x);
y2=cos(x);
figure;
subplot(1,2,1);
plot(x,y1,'--r');
legend('sin');
title('sin');
xlabel('sin(x)');
subplot(1,2,2);
plot(x,y2,'-.k');
legend('cos');
title('cos');
xlabel('cos(x)');
hs=suptitle('∂‘±»Õº');
set(hs,'FontSize',16);

end

