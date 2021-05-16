function  [c,t]= poly_coef(fcn,var)
% 输入：
%      fcn 为待确定系数值多项式的符号表达式
%      var 为多项式的主变量，如果是单变量可省略之。
% 输出：
%      c 为系数向量，所对应的变量幂次有高到低。
%      t 与系数向量相对应项的变量向量，即 fcn=c*conj(t)'
syms n
if nargin<2
    var=findsym(fcn);
    if length(var)>1
        error('Specify var please!')
    end
end
[~,t0]=coeffs(fcn,var);
m=log2(subs(t0(1),2));
f=fcn+symsum(var^n,n,0,m);
[c1,t]=coeffs(f,var);
c=c1-ones(1,length(c1));