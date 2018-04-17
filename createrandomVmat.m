function [Vmat] = createrandomVmat(para)
% Creates random V-matrix used for optimal boson basis.
% Dimensions per site: dk x d_opt
% Modified
%	FS 16/07/2015: - use prod(para.dk(:,I)) to allow multiple chains at OBB level
%	FS 21/08/2015: - Can also create Vtens for Multi-Chain code

Vmat = cell(1, para.L);

for i = 1:para.L
    if para.useVmat == 1 && prod(i~= para.spinposition)				% ready for array in spinposition
        if para.parity~='n'
			assert(para.nChains == 1, 'only single chains allowed with this parity');
            halfdkdim=para.dk(i)/2;
            halfddim=para.d_opt(i)/2;
            Vmat{i} = zeros(para.dk(i), para.d_opt(i));
            Vmat{i}(1:halfdkdim,1:halfddim)=randn(halfdkdim,halfddim)./sqrt(halfdkdim*halfddim);
            Vmat{i}(halfdkdim+1:end,halfddim+1:end)=randn(halfdkdim,halfddim)./sqrt(halfdkdim*halfddim);
		else
			if para.useVtens == 0
				Vmat{i} = randn(prod(para.dk(:,i)), para.d_opt(i));
				Vmat{i} = Vmat{i}./sqrt(numel(Vmat{i}));
			else
				NC = para.nChains;
				Vmat{i} = cell(1,NC+1);
				for mc = 1:NC
					Vmat{i}{mc} = randn(para.dk(mc,i),para.d_opt(mc,i));
					Vmat{i}{mc} = Vmat{i}{mc}./sqrt(numel(Vmat{i}{mc}));
				end
				Vmat{i}{end} = randn(para.d_opt(:,i)');
				Vmat{i}{end} = Vmat{i}{end}./sqrt(numel(Vmat{i}{end}));
			end
        end
    else
        Vmat{i} = eye(para.dk(1,i), para.d_opt(1,i));                   % why still para.d_opt??
    end
end
end
