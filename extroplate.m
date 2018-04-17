function [epsilon,t]=extroplate(epsilon,t,L)
%This function is used to get Wilson chain parameters if the chain lenth is very long and the Wilson chain parameters
%at the tail of chain are much smaller than machine precision, where the normal numerics won't get the correct results.

wp=1e-30;                       %Wrong precision. The inputing vectors goes wrong below this presion therefore need to be corrected.
a=max(find(epsilon>wp*1e5));
b=max(find(epsilon>wp));
if L>b
x=a:b;
y=log(epsilon(a:b));
p=polyfit(x',y,1);

l=b+1:L;
epsilon(l)=exp(p(1).*l+p(2));

a=max(find(t>wp*1e5));
b=max(find(t>wp));
x=a:b;
y=log(t(a:b));
p=polyfit(x',y,1);

l=b+1:L;
t(l)=exp(p(1).*l+p(2));

end
end
