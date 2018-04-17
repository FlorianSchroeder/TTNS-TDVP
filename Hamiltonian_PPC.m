function [H0, H1] = Hamiltonian_PPC(para)
% Creates the Hamiltonian for a PPC. Uses Cyclic symmetry
%   H entries in units of eV
%
% Modified:
%   FS 20140609:    - Can apply diagonal static disorder using para.MLSB_disDiag and para.MLSB_staticDisorder
%
%
% zero disorder Hamiltonian:
switch para.MLSB_system
    case 'RsMolischianumB850'
        H0 = RsMolischianumB850();
    case 'RsMolischianumB800B850'
        H0 = RsMolischianumB800B850();
    otherwise
end

% apply disorder:
if isfield(para,'MLSB_staticDisorder') && para.MLSB_staticDisorder == 1
    % apply diagonal disorder
    H0 = H0 + diag(cmToeV(para.MLSB_disDiag));
end

if para.MLSB_tOff == 0
    H1  = diag(para.MLSB_t(1:size(H0)));
        % produces a Diagonal matrix, defining the coupling between each system site and the first Wilson chain site.
        % only define relative to para.t(1)
        % para.MLSB_t should have some periodicity, n = #states;
        % also works for MLSB_t = vector
elseif para.MLSB_tOff >= 1
        % para.MLSB_tOff defines the on which off diagonal the coupling should be applied
        % could be combined with above code.
    dim = size(H0,1);
    H1  = diag(para.MLSB_t(1:dim));
    H1  = [zeros(dim,para.MLSB_tOff), H1];
    b1  = H1(:,dim+1:end);            % coupling to close ring
    H1  = H1(1:dim,1:dim)+padarray(b1',[dim-size(b1',1) 0],'post');
% comment following line if jaynes-cummings type coupling is desired!
% Symmetrize H1:
    H1  = H1+H1' -diag(diag(H1));            % right-triangular to symmetric H0
end
end

function H0 = RsMolischianumB850()
    %% from Tretiak 2000 for Rs. molischianum
    H0para.n_unit = 2;                 % number of states/sites per unit(alpha+beta bchl a)
    H0para.n_symmetry = 8;             % cyclic symmetry in pigment complex
    H0para.H_dim = H0para.n_unit*H0para.n_symmetry;  % number of excited states to consider
    H0para.n_interaction_range = 3;    % # interacting unit cells

    %% Define matrix elements in eV, right-diagonal Form
    H0para.matrix_unit_cell = [1.59           cmToeV(258);...   % Hamiltonian of unit cell. Dim: n_unit x n_unit
                               0              1.61];            % Diagonal similar to Hu 1997: 1.6191eV

    % in cm^-1
    H0para.matrix_cell_coupling = [-67    22      0       0;... % Hamiltonian with cell interactions. Dimensions: n_unit x n_unit*(n_interaction_range-1)
                                   210    -40     17      0];
    H0para.matrix_cell_coupling = cmToeV(H0para.matrix_cell_coupling);

    %% Construct Ring Hamiltonian:
    H0 = RingHamiltonian(H0para);
%    [V,D] = eig(H0);              % get eigenstates in V and eigenvalues in D
%     plot(diag(D)) % to plot the eigenvalues. See 20140520 presentation.
%     figure(2)
%     plot(V.*V)

end

function H0 = RsMolischianumB800B850()
%% TO complete

end

function H0 = RingHamiltonian(H0para)
    % creates hamiltonian with pbc
    % H0para: struct with fields:
    %   n_unit, n_symmetry, H_dim, n_interaction_range, matrix_unit_cell, matrix_cell_coupling

    H = zeros(H0para.H_dim,H0para.H_dim+H0para.n_unit*(H0para.n_interaction_range-1));   % make longer to accomodate interaction terms
    for i=1:H0para.n_unit:H0para.n_symmetry*H0para.n_unit
        H(i:i+H0para.n_unit-1, i:i+H0para.n_unit*H0para.n_interaction_range-1) = ...
            [H0para.matrix_unit_cell H0para.matrix_cell_coupling];
    end
    b = H(:,H0para.H_dim+1:end);            % coupling to close ring
    H0 = H(1:H0para.H_dim,1:H0para.H_dim)+padarray(b',[H0para.H_dim-size(b',1) 0],'post');
    H0 = H0+H0' -diag(diag(H0));            % right-triangular to symmetric H0
end