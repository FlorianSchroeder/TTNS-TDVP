function op=gen_sitej_op(op,para,sitej,leftge)
% Take: 	op.Opstorage(:,1,j) 	-> op.Opleft	*rescaling: x := Lambda^(j-2) * x
%			op.Opstorage(:,2,j+1)	-> op.Opright
%			op.Hlrstorage{j}		-> op.Hleft     *rescaling: x := Lambda*(x - leftge);
%			op.Hlrstorage{j+1}      -> op.Hright	*rescaling: x := Lambda^(j-2) * x
% Rescale if para.rescaling == 1
%
% leftge = 0 in first loop?
% Commented by Florian Schroeder 30/01/2014
%   Changed
%       - FS 20/10/2014: introduced para.sweepto for operator selection

switch para.sweepto
    case 'r'
        op.Opleft = op.Opstorage(:,1,sitej);				% empty in first l -> r sweep.	From updateop()
        op.Opright = op.Opstorage(:,2,sitej+1);             % interaction term op of next site. In eff basis of current r_j
        op.Hleft = op.Hlrstorage{sitej};					% =0 for j = 1 in first sweep.	From updateop()
        op.Hright = op.Hlrstorage{sitej + 1};				% H_r which is non-interacting with site j, In eff basis of current r_j

        % Rescale if wanted
        if sitej >= 3 && isfield(para,'rescaling') && para.rescaling
            Lambda = para.chain{1}.Lambda;
            %old: op.Hleft =Lambda.^(sitej-2).*(op.Hlrstorage{sitej}-leftge(sitej).*eye(size(op.Hlrstorage{sitej})));

            op.Hleft = Lambda.*(op.Hleft - leftge(sitej).*eye(size(op.Hleft)));
            op.Hright = Lambda.^(sitej-2).*op.Hright;
            for m=1:para.M
                op.Opleft{m} = op.Opleft{m}.*Lambda^(sitej-2);
            end
        end

        % Get rescaled bare h1 and h2 terms of current site in on-site basis. (of h2term, only {m,1,j} are rescaled)
        op = gen_sitej_h1h2(op,para,sitej);
    case 'l'
        % Does exactly same as 'r' sweep
        op.Opleft = op.Opstorage(:,1,sitej);				% empty in first l -> r sweep.	From updateop()
        op.Opright = op.Opstorage(:,2,sitej+1);             % interaction term op of next site. In eff basis of current r_j
        op.Hleft = op.Hlrstorage{sitej};					% =0 for j = 1 in first sweep.	From updateop()
        op.Hright = op.Hlrstorage{sitej + 1};				% H_r which is non-interacting with site j, In eff basis of current r_j
        op = gen_sitej_h1h2(op,para,sitej);
end

end
