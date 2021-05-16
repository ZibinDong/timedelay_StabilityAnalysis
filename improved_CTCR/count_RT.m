function [RT] = count_RT(s_p,t_p, params_a)
%COUNT_RT 此处显示有关此函数的摘要
%   此处显示详细说明
syms s;
a = 0;
b = 0;
for j = 1:length(params_a)
    a = a + (subs(diff(params_a(j), s),s, s_p)*exp(-j*t_p*s_p));
end
for j = 1:length(params_a)
    b = b + (j*subs(params_a(j),s,s_p)*exp(-j*t_p*s_p));
end
RT = imag(a/b);
if RT>=0
    RT = 1;
else
    RT = -1;
end

end

