function op=update_sitej_h1h2(para,op,sitej)
% 1. Updates Hamiltonian terms with spatial shift para.shift(j)
% 2. Takes new terms into h1j, h2j and rescales
% Commented by Florian Schroeder 13/01/2014

op=genh1h2term_onesite(para,op,sitej);
op=gen_sitej_h1h2(op,para,sitej);

end