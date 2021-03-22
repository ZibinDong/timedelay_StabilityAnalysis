function [Tck, NSs] = count_T(b_charac, n)
%COUNT_T 此处显示有关此函数的摘要
%   此处显示详细说明
syms T;
NS = 0;
NSs = [];
Tck = [];
now = 0;
for Tt = -1:0.0001:1
    now = now+0.0001;
    comp = now/(2);
    disp(comp);
    T_params = subs(b_charac, T, Tt);
    T_params = sym2poly(T_params);
    routh_list = Routh(T_params);
    routh_list = routh_list(:,1);
    temp_NS = 0;
    for i = 1:length(routh_list)-1
        if routh_list(i)*routh_list(i+1) < 0
            temp_NS = temp_NS+1;
        end
    end
    if temp_NS ~= NS
        Tck = [Tck, Tt];
        NSs = [NSs, temp_NS];
        NS = temp_NS;
    end
end
Tck = Tck(2:end);
NSs = [NSs, NSs(1)];
NSs = NSs(2:end);
if ismember(0, Tck)
    index = find(Tck==0);
    if index~=1 && abs(NSs(index-1)-NSs(index)) == n
        Tck(index) = [];
        NSs = [NSs, NSs(index)];
        NSs(index) = [];
    end
end


end

