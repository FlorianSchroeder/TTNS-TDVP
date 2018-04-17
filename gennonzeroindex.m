function para=gennonzeroindex(mps,Vmat,para,j)
if para.parity~='n'
    if nargin==4
            para.Anzi{j}=find(mps{j});
            [ll,rr,dd] = size(mps{j});
            assert(length(para.Anzi{j})==ll*rr*dd/2);
            para.Vnzi{j}=find(Vmat{j});
            assert(length(para.Vnzi{j})==para.dk(j)*para.d_opt(j)/2);


            A=permute(mps{j},[3,1,2]);
            para.Anzilr{j}=find(A);

            A=permute(mps{j},[2,3,1]);
            para.Anzirl{j}=find(A);
    else
        for j=1:para.L
            para=gennonzeroindex(mps,Vmat,para,j);
        end
    end
end
end
