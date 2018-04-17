function vNE = vonNeumannEntropy(S)
% Calculates for eigenvalues s_i : vNE = -  sumOver_i(s_i^2 ln(s_i^2))

dim=size(S,1);			% = dk
s=diag(S);
s=s./sqrt(s'*s);			% normalize eigenvalues (probabilities)

%maxvNE=log(dim)/log(2.0);

entropy=0;
for x=1:dim
    if(s(x)>1e-16)
        entropy=entropy-s(x).^2.*(log(s(x).^2));
    end
end

vNE=entropy;
