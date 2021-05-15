clc;clear
% ����������ʼ��
A = input('A = ');
B = input('B = ');
n = size(A, 1);
syms s t x T f(s,t) g(s,T);
f = exp(-s*t);
g = (1-T*s)/(1+T*s);
I = eye(n);

% �����������̣�����aϵ��
charac = det(s*I - A - B*x);
params_a = coeffs(charac, x);
b_charac = subs(charac, x, g);
b_charac = collect(b_charac, s);
b_charac = numden(b_charac);
b_charac = collect(b_charac, s);

% ����Tck
disp("counting Tck");
[Tck,NSs] = count_T(b_charac, n);
m = length(Tck);

% ���wck
disp("counting wck");
wck = rand(1,m);
for k = 1:m
    roots_vector = roots(sym2poly(subs(b_charac, T, Tck(k))));
    for j = 1:length(roots_vector)
        if abs(real(roots_vector(j))) < 0.001
            wck(k) = abs(imag(roots_vector(j)));
            break
        end
    end
end

% ����1���ڵ���Tck������tkl
disp("counting tkl");
tkl = [];
for k = 1:m
    l = 0;
    while (2/wck(k))*(atan(wck(k)*Tck(k))+l*pi) <= 0
        l= l+1;
    end
    new_t = (2/wck(k))*(atan(wck(k)*Tck(k))+l*pi);
    tkl = [tkl, new_t];
end

tkl_o = tkl;
wck_o = wck;
Tck_o = Tck;

for k = 1:m
    l = 1;
    while (tkl(k) + 2*pi/wck(k)*l) < 1
        tkl = [tkl, tkl(k) + 2*pi/wck(k)*l];
        wck = [wck, wck(k)];
        Tck = [Tck, Tck(k)];
        l = l+1;
    end
end
m = length(Tck);

% �ϲ��ɱ�񲢰���tkl��������
sheet = [tkl', wck', Tck'];
sheet = sortrows(sheet,1);
tkl = sheet(:,1);
wck = sheet(:,2);
Tck = sheet(:,3);

% ����RT
RTs = [];
RTs_o = [];
for k = 1:m
    RTs = [RTs; count_RT(wck(k)*1i, tkl(k), params_a)];
end
for k = 1:length(Tck_o)
    RTs_o = [RTs_o; count_RT(wck_o(k)*1i, tkl_o(k), params_a)];
end
sheet = [tkl, RTs, wck, Tck];

% ����NU(0)
NU_0 = 0;
t_charac = subs(charac, x, f);
t_charac = subs(t_charac, t, 0);
roots_0 = roots(sym2poly(t_charac));
for j = 1:length(roots_0)
    if real(roots_0(j)) >= 0
        NU_0 = NU_0+1;
    end
end

% ����NU_t
NU_t = [];
for k = 1:length(tkl)
    NU_t = [NU_t; count_NU(tkl(k)+0.0001, wck_o, tkl_o,RTs_o,NU_0)];
end

disp("sheet=t  RT  w  T");
disp(sheet);
disp("NU(t+delta)");
disp(NU_t);

