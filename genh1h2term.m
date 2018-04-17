function [op,para]=genh1h2term(para,op)
% Generate the hamiltonian terms

op.h1term = cell(1,para.L); % One body terms of the Hamiltonian

% two bodies terms of the Hamiltonian,
% op.h2term(:,1,j) is the operators of site j when it coupled to site j+1, and
% op.h2term(:,2,j) is the operater when it coupled to j-1.
% The factors of the Hamiltonian terms are included only in the op.h2term(:,1,j)

op.h2term = cell(para.M, 2,para.L);			% default: nearest neighbor interaction

for s=1:para.L
    op=genh1h2term_onesite(para,op,s);
end

end
