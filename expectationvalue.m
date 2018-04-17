function [e] = expectationvalue(hset,mps1,V1,mps0,V0)
% Calculates the overlap <mps1|hset|mps0>
%	hset:	MxN cell
%		expectation values taken over each row in M seperately and summed in the end.
%		N does not have to equal para.L
%
% can be used for multi-site expectation values!
% If hset is empty then it's eqvilent to calculate the overlap.
%
% Modified:
%		FS 02/06/2014:	- changed em=1 to em = [] as starting value.

[M, N] = size(hset);

% expectation value
e = 0;
for m = 1:M
    em = [];					%	allows starting somewhere in the middle of the chain.
    for j = N:-1:1
        em = updateCright(em, mps0{j},V0{j}, hset{m,j}, mps1{j},V1{j});
    end
    e = e + em;
end

% % norm
% n0=1;
% n1 = 1;
%
% for j = N:-1:1
%     d = size(mps0{j}, 3);
%     X = eye(d);
%    n0 = updateCright(n0, mps0{j}, [],X, mps0{j},[]);
% end
%
% for j = N:-1:1
%     d = size(mps1{j}, 3);
%     X = eye(d);
%    n1 = updateCright(n1, mps1{j}, [],X, mps1{j},[]);
% end
%
% e = e/sqrt(n1)/sqrt(n0);


end