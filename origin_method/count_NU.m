
function [NU_t] = count_NU(t, wck_o,tkl_o, RTs_o, NU_0)
%COUNT_NU 此处显示有关此函数的摘要
%   此处显示详细说明
NU_t = NU_0;
U = 0;
for k = 1:length(wck_o)
    if t<tkl_o(k)
        U = 0;
    elseif t>tkl_o(k) && wck_o(k) == 0
        U = 1;
    elseif t>tkl_o(k) && wck_o(k) ~= 0
        U = 2;
    end
    NU_t = NU_t + (ceil((t-tkl_o(k))/(2*pi/wck_o(k)))*U*RTs_o(k));
end

