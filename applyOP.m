function mpsnew=applyOP(OP,mps)
L=length(mps);
assert(L==length(OP));
mpsnew=cell(L,1);
for j=1:L
    mpsnew{j}=contracttensors(mps{j},3,3,OP{j},2,2);
end
end