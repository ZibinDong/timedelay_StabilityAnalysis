clc;clear
% collect input data
% A = [-1 13.5 -1;-3 -1 -2;-2 -1 -4]
% B = [-5.9 7.1 -70.3;2 -1 5;2 0 6]
A = input('A = ');
B = input('B = ');
time_delay_limit = 1;
n = size(A, 1);
syms s t x T f(s,t) g(s,T);
f = exp(-s*t);
g = (1-T*s)/(1+T*s);
I = eye(n);

% count charac_eq
charac = det(s*I - A - B*x);
params_a = coeffs(charac, x);
b_charac = subs(charac, x, g);
b_charac = collect(b_charac, s);
b_charac = numden(b_charac);
b_charac = collect(b_charac, s);

% count Tck
disp("counting Tck");
b_coef = poly_coef(b_charac, s);
Routh_table = Routh(b_coef);
R1 = numden(Routh_table(end-1,1));
R21 = Routh_table(end-2,1);
R22 = Routh_table(end-2,2);
T_roots = vpa(solve(R1));
Tck = [];
for i = 1:length(T_roots)
    if(abs(imag(T_roots(i))) < 0.0001)
        Tck = [Tck, T_roots(i)];
    end
end
digits(5);
Tck = vpa(Tck);
m = length(Tck);

% count wck
disp("counting wck");
wck = zeros(1,m);
for k = 1:m
    R21_tmp = subs(R21, T, Tck(k));
    R22_tmp = subs(R22, T, Tck(k));
    if R21_tmp*R22_tmp > 0
        wck(k) = sqrt(R22_tmp/R21_tmp);
    else
        disp('error Tck & wck');
        exit();
    end
end

% count min tkl
disp("counting tkl");
tkl = zeros(1,m);
for k = 1:m
    l = 0;
    while (2/wck(k))*(atan(wck(k)*Tck(k))+l*pi) <= 0
        l= l+1;
    end
    while (2/wck(k))*(atan(wck(k)*Tck(k))+(l-1)*pi) >= 0
        l = l-1;
    end
    new_t = (2/wck(k))*(atan(wck(k)*Tck(k))+l*pi);
    tkl(k) = new_t;
end
% count RT
RTs = zeros(1,m);
for k = 1:m
    RTs(k) = count_RT(wck(k)*1i, tkl(k), params_a);
end

% build map sheet [tkl, RTk, wck, Tck]
map_sheet = [tkl', RTs', wck', Tck'];

expand_m = m;
result_table = map_sheet;
% expand tkl 
for k = 1:m
    l = 1;
    while tkl(k)+2*pi/wck(k)*l < time_delay_limit
        result_table = [result_table; result_table(k, :)];
        result_table(end, 1) = tkl(k)+2*pi/wck(k)*l;
        l = l+1;
        expand_m = expand_m + 1;
    end
end
result_table = sortrows(result_table, 1);


% count NU(0)
NU_0 = 0;
t_charac = subs(charac, x, f);
t_charac = subs(t_charac, t, 0);
roots_0 = roots(sym2poly(t_charac));
for j = 1:length(roots_0)
    if real(roots_0(j)) >= 0
        NU_0 = NU_0+1;
    end
end

% count NU_t
NU_t = zeros(expand_m, 1);
for k = 1:expand_m
    NU_t(k) = count_NU(result_table(k,1)+0.0001, map_sheet, NU_0);
end

disp("sheet=t  RT  w  T");
disp(double(result_table));
disp("NU(t+delta)");
disp([NU_0;NU_t]);

