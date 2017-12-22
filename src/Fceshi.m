function Fceshi(n)
% �������д��2014��8��20�գ�����LegendreP�ƽ�(t^2+u^2)^(-1/2)
% n ȡ��LegendreP�Ľ���

t=0.1:0.2:2;
m=length(t);
u=0:0.05:1;
mu=length(u);
for i=1:m
    an=Fan(n,t(i));
    y=zeros(n+1,mu);
    for j=1:n+1
        leg=legendre(j-1,u);
        y(j,:)=leg(1,:);
    end
    py=zeros(1,mu);
    for k=1:mu
        py(k)=an*y(:,k);
    end
    % ��ֵ
    real_y=(t(i)^2+u.^2).^(-1/2);
    error=py-real_y;
    figure;
    plot(u,error,'*-');
    title(['t=',num2str(t(i)),',Legendre ����Ϊ ',num2str(n)]);
end

end

% ���� a_n(t)
function value=Fan(n,t)

value=zeros(1,n+1);
for i=0:n
    value(i+1)=(i+1/2)*Gauss32(-1,1,i,t);
end

end