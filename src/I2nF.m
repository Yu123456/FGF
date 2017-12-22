function value=I2nF(I0,Y,n)
% 本程序编写于2014年8月17日，硕士论文，手稿P102中式（102-2），递推计算
% I0 上一次递推值
% Y 
% n 指标

n2=2*n;
Y2=Y^(n2-1);
value=Y2*Y-n2*Y2+n2*(n2-1)*I0;

end

