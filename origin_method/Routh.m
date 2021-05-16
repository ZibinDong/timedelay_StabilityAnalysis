function [routh_list] = Routh(chara_equ)
%UNTITLED �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
n = length(chara_equ);
chara_equ = reshape(chara_equ, 1 ,n);
if mod(n,2) == 0
    n1 = n/2;
else
    n1 = (n+1)/2;
    chara_equ = [chara_equ, 0];
end
routh = reshape(chara_equ, 2, n1);
routh_list = zeros(n, n1);
routh_list(1:2, :) = routh;
i = 3;
while 1
    if routh_list(i-1,1)==0 && sum(routh_list(i-1,2:n1)~= 0)
        chara_equ = conv(chara_equ, [1 3]);
        n = length(chara_equ);
        if mod(n,2)==0
            n1 = n/2;
        else
            n1 = (n+1)/2;
            chara_equ = [chara_equ, 0];
        end
        routh = reshape(chara_equ,2,n1);
        routh_list=zeros(n,n1);
        routh_list(1:2, :) = routh;
        i = 3;
    end
    ai = routh_list(i-2,1)/routh_list(i-1,1);
    for j=1:n1-1
        routh_list(i,j)=routh_list(i-2,j+1)-ai*routh_list(i-1,j+1);
    end
    if sum(routh_list(i,:)) == 0
        k=0;
        l=1;
        F = zeros(1,n1);
        while n-i-k>=0
            F(1)=n-i+l-k;
            k = k+2;
            l = l+1;
        end
        routh_list(i,:)=routh_list(i-1,:).*F(1,:);
    end
    i = i+1;
    if i>n
        break
    end
end

end

