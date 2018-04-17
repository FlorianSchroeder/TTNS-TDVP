function [Cstorage] = initCstorage(mpsB, mpoX, mpsA, N)
    Cstorage = cell(1, N + 1);
    Cstorage{1} = 1;
    Cstorage{N + 1} = 1;
    for i = N:-1:2
        if isempty(mpoX), X = [];
        else X = mpoX{i};
        end
        Cstorage{i} = updateCright(Cstorage{i + 1}, mpsB{i}, X, mpsA{i});
    end