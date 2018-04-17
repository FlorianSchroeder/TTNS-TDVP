function spin=calspin(mps,Vmat,para,results)
%   DEPRECATED: Use getObservable()
%
%Calculate the spin expectation value
% Has to be modified if site 1 changes dimension!!
% Modified:
%   FS 24/05/2014:  excluded MLSpinBoson. Use with try-catch block

N=para.L;
assert(N==length(mps) && N==length(Vmat));
assert(~strcmp(para.model,'MLSpinBoson'),'not possible for MLSBM');

ndset=cell(1,N);

for j=1:N
    ndset{1,j}=eye(size(Vmat{j},1));
end
sx=ndset;sy=sx;sz=sy;

%debug:
sx{1,1}

[sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
sx{para.spinposition}=sigmaX;
sy{para.spinposition}=sigmaY;
sz{para.spinposition}=sigmaZ;
if strcmp(para.model,'2SpinPhononModel')
    sx{para.spinposition}=kron(sigmaZ,eye(2));  %measures excitation of site1
    sz{para.spinposition}=kron(eye(2),sigmaZ);  %measures excitation of site2
    %try to find a good way to measure this!
    sy{para.spinposition}=kron(sigmaY,eye(2));  %measures only sy of site1
end
spin.sx=expectationvalue(sx,mps,Vmat,mps,Vmat);
spin.sy=expectationvalue(sy,mps,Vmat,mps,Vmat);
spin.sz=expectationvalue(sz,mps,Vmat,mps,Vmat);
spin.sx=real(spin.sx);
spin.sy=real(spin.sy);
spin.sz=real(spin.sz);

end
