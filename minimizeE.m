function [mps,Vmat,para,results,op] = minimizeE(op,para)
% Ground state calculation
% Modified:
%	FS 30/05/2014: - First verion of a energy flowdiagram added. See Cheng 2012 for more information.
%   FS 06/06/2014: - 'para.loop == 20 && size(results.d_opt,2) < 2' added to adjustdopt.m criterium, as wrongly chosen
%                    initial D, d_opt or dk can prevent convergence!
%   FS 20/10/2014: - using para.sweepto everywhere.
%
% Commented ba Florian Schroeder 03/02/2014

% randn('state', 0)   %% why this line? this always generates the same random Vmat and MPS
rng(0,'v5normal');     %% should be replaced by: rng(0);
L = para.L;
M = size(op.h2term,1);

if para.resume == 1 && para.savedexist==1
    [Vmat,mps,loop,para,results,op]=loadsaved(para);
else
    loop=1;
    Vmat    = createrandomVmat(para);
    mps     = createrandommps(para);
    %Preassign space for the results structure.
    results = initresults(para);
end

para = gennonzeroindex(mps,Vmat,para);
para

%% Override if preparing artificial vacuum Ground State
if (strcmp(para.model,'SpinBoson') || strcmp(para.model,'SpinBoson2folded')|| strcmp(para.model,'SpinBoson2C')) && strcmp(para.SpinBoson.GroundStateMode,'artificial')
	if strcmp(para.SpinBoson.InitialState, 'sz')
		%% prepare +Sz eigenstate
		mps{1} = reshape([1,zeros(1,numel(mps{1})-1)],[1,para.D(1),para.d_opt(1)]);
		Vmat{1} = eye(para.dk(1));
	elseif strcmp(para.SpinBoson.InitialState, '-sz')
		mps{1} = reshape([  zeros(1,numel(mps{1})/2),...
						  1,zeros(1,numel(mps{1})/2-1)],[1,para.D(1),para.d_opt(1)]);
		Vmat{1} = eye(para.dk(1));
	elseif strcmp(para.SpinBoson.InitialState, 'sx')
		mps{1} = reshape([1/sqrt(2),zeros(1,numel(mps{1})/2-1),...
						  1/sqrt(2),zeros(1,numel(mps{1})/2-1)],[1,para.D(1),para.d_opt(1)]);
		Vmat{1} = eye(para.dk(1));
	elseif strcmp(para.SpinBoson.InitialState, '-sx')
		mps{1} = reshape([-1/sqrt(2),zeros(1,numel(mps{1})/2-1),...
						   1/sqrt(2),zeros(1,numel(mps{1})/2-1)],[1,para.D(1),para.d_opt(1)]);
		Vmat{1} = eye(para.dk(1));
	elseif strcmp(para.SpinBoson.InitialState, 'sy')
		mps{1} = reshape([ 1/sqrt(2),zeros(1,numel(mps{1})/2-1),...
						  1i/sqrt(2),zeros(1,numel(mps{1})/2-1)],[1,para.D(1),para.d_opt(1)]);
		Vmat{1} = eye(para.dk(1));
	elseif strcmp(para.SpinBoson.InitialState, '-sy')
		mps{1} = reshape([-1/sqrt(2),zeros(1,numel(mps{1})/2-1),...
						  1i/sqrt(2),zeros(1,numel(mps{1})/2-1)],[1,para.D(1),para.d_opt(1)]);
		Vmat{1} = eye(para.dk(1));
	else
		error('VMPS:minimizeE:DefineInitialState','InitialState=none is not implemented yet');
	end
	for j = 2:para.L
		if j == para.L
			Dr = 1;
		else
			Dr = para.D(j);
		end

		mps{j} = reshape([1, zeros(1,numel(mps{j})-1)],para.D(j-1),Dr,para.d_opt(j));
		if para.foldedChain
			locDim      = sqrt(para.dk(j));
			occEstimate = kron(locDim:-1:1,ones(1,locDim))+kron(ones(1,locDim),locDim:-1:1);		% kronecker sum to estimate lowest states
			[~,order]   = sort(occEstimate);														% find order of lowest states
		else
			order = para.dk(j):-1:(para.dk(j)-para.d_opt(j)+1);										% inverted order
		end
		if para.nChains == 1
% 			Vmat{j}	= sparse(order(1:para.d_opt(j)),1:para.d_opt(j),1,para.dk(j),para.d_opt(j));	% inverted order
			Vmat{j}	= sparse(1:para.d_opt(j),1:para.d_opt(j),1,para.dk(j),para.d_opt(j));			% 1:n order
		else
			Vmat{j} = sparse(1,1,1, prod(para.dk(:,j)),para.d_opt(j));								% 1:n order
% 			error('VMPS:minimizeE:artificial','not yet done for Multi-Chain Vmat');
		end
	end
	[mps,Vmat,para] = prepare(mps,Vmat,para);
	[op] = initstorage(mps, Vmat, op,para);
	para.trustsite = para.L;		% needed for TDVP
	return;
	% exit the function here since no optimization is necessary!
end


[mps,Vmat,para] = prepare(mps,Vmat,para);
% storage-initialization sweep l <- r
[op] = initstorage(mps, Vmat, op,para);

%The relative change of d_opt and D, 1 means 100%.
para.d_opt_change=1;
para.D_change=1;

while loop<=para.loopmax;
    fprintf('\nloop = %d\n',loop);
    para.loop=loop;
	para.shifted = 0;													% to log if bosons were shifted
    results.Evalues = [];
    results.flowdiag{loop} = [];
    % ********* cycle 1: j ? j + 1 (from 1 to L - 1)*********
    for j = 1:L
        fprintf('%g-', j); para.sweepto = 'r'; para.sitej = j;
        op=gen_sitej_op(op,para,j,results.leftge);  					% take Site h1 & h2 Operators apply rescaling to Hleft, Hright, Opleft ...???
        [Amat,Vmat,results,para,op]   = optimizesite(mps,Vmat,op,para,results,j);
		if j~=L
            [mps{j}, U, para,results] = prepare_onesite(Amat,para,j,results);
            mps{j+1} = contracttensors(U,2,2,mps{j+1},3,1);
        else
            mps{j}=Amat;
		end
		% optimisation finished -> update effective Hamiltonian operators
        op = updateop(op,mps,Vmat,j,para);                                % calls updateHleft, updateCleft to update for next sweep
        % sweep finished

		% get Energy eigenvalues for analysis
        Hleft     = op.Hlrstorage{j+1};
		eigvalues = sort(eig((Hleft'+Hleft)/2));
        % log the flowdiagram
        if para.logging && j~=L                                                         % Last site has different length in EV
            if j~=1            % check for different dimensions to apply vertcat
                dimCat = min(size(results.flowdiag{loop},2),size(eigvalues,1));
                results.flowdiag{loop} = [results.flowdiag{loop}(:,1:dimCat); eigvalues(1:dimCat)'];
            else
                results.flowdiag{loop} = [results.flowdiag{loop}; eigvalues'];
            end
        end
		results.leftge(j+1)=eigvalues(1);
        results.Evalues = [results.Evalues, results.E];		% results.E from optimizesite()
    end
	if (para.logging && para.shifted) 		% log the shift if shifted (trustsite > 0)
		results.shift{para.loop} = para.shift;
	end
	if para.logging
		results.Vmat_svLog{para.loop} = cellfun(@(x) x(1,1), results.Vmat_sv(~cellfun('isempty',results.Vmat_sv)));		% log the highest SV of Vmat
		results.EvaluesLog{para.loop} = results.Evalues;
	end
	if para.useFloShift3
		para.FloShift3minMaxSV = min(cellfun(@(x) x(1,1), results.Vmat_sv(~cellfun('isempty',results.Vmat_sv))));
	end
    fprintf('\nE = %.10g\t', results.E);							% print last energy Eigenvalue calculated
    results.Eerror(loop)=std(results.Evalues)/abs(mean(results.Evalues));
    fprintf('E_error =  %.16g\t',results.Eerror(loop));
    %Calculate the relative change of Vmat von Neumann entropy as one criteria for convergence check
    % TODO: Need to exclude Vmat_vNE == 0 as this gives either inf or NaN
    vNEdiff=(results.Vmat_vNE(2:end)-results.lastVmat_vNE(2:end))./results.Vmat_vNE(2:end);
    vNEdiff(results.Vmat_vNE(2:end)==0) = 0;        % if Vmat_vNE == 0 set diff to 0 to fix inf & NaN
    vNEdiff=abs(vNEdiff);
    fprintf('\nvNEdiff = \n'); disp(mat2str(vNEdiff));
    results.vNEdiff=vNEdiff;
    fprintf('para.shift = \n'); disp(mat2str(para.shift));
    results.lastVmat_vNE=results.Vmat_vNE;
    para=trustsite(para,results);
    fprintf('precise_sites k <= %d\t',para.precisesite);
    fprintf('trustable_sites k <= %d\n',para.trustsite(loop));

    %% start sweep l <- r
    para.sweepto = 'l';

    if (para.dimlock==0 && para.trustsite(end)>3 && mod(loop,2)==0) || results.Eerror(end)<para.precision || para.loop == 20 && size(results.d_opt,2) < 2 %&& para.trustsite(end)/para.precisesite>0.6 && sqrt(var(para.trustsite(end-2:end)))/mean(para.trustsite(end-2:end))<0.1 %The trust site has not been improved during the last 2 sweeps
        if para.loop == 20 && size(results.d_opt,2) < 2
            para.trustsite(end) = 5;        % arbitrary setting to cause an optimization!
        end
        %%Expand or Truncate D and d_opt
        para.adjust = 1;
%		 if results.Eerror(end)<para.precision						% PERHAPS GOOD STATEMENT!! TEST THIS! ADD 1 sweep after that before exit!
%			para.trustsite(end) = para.L;								% this ensures last adjustment of dimensions before finishing minimizeE.m
%		 end
        [op,para,results,mps,Vmat] = adjustdopt(op,para,results,mps,Vmat);      % optimise d_opt and dk, calls genh1h2term -> all effective H are invalid now
        [mps,Vmat,para, results]   = rightnormA(mps,Vmat,para,results);         % changed to also update results
        %para.trustsite(loop)=0;
        fprintf('d_opt = ');
        disp(mat2str(para.d_opt));
        fprintf('para.D = ');
        disp(mat2str(para.D));
		dispif('para.dk = ', para.useDkExpand);
        dispif(mat2str(para.dk),para.useDkExpand);
     else
        para.adjust=0;
        [mps,Vmat] = rightnormA(mps,Vmat,para,results);			%	sweeps r->l to right-normalize A matrices.
    end

    [op] = initstorage(mps, Vmat, op,para);							%	recreate operators for next optimization sweep

    fprintf('d_opt_change=%g\t',para.d_opt_change);
    fprintf('D_change=%g\t',para.D_change);
    if abs(para.d_opt_change)<para.minDimChange && abs(para.D_change)<para.minDimChange
      para.dimlock=1;
    end
    results.maxVmatsv=max(cellfun(@(x) x(end), results.Vmat_sv(max(para.spinposition+1):para.trustsite(end))));fprintf('maxVmatsv=%g\t',results.maxVmatsv);
    results.maxAmatsv=max(cellfun(@(x) x(end), results.Amat_sv(2:para.trustsite(end))));fprintf('maxAmatsv=%g\n',results.maxAmatsv);
    %if para.trustsite(loop)>=para.L-2 && max(vNEdiff)<1e-3 && results.maxVmatsv<para.svmintol && results.maxAmatsv<para.svmintol
    if (results.Eerror(end)<para.precision && para.dimlock == 1) || (para.loop > 10 && prod(results.Eerror(end-10:end)<para.precision)) %&& para.trustsite(end)>para.L-5
		break;
    end

    %if para.trustsite(end)>para.L-10 && results.maxVmatsv<para.svmaxtol && results.maxAmatsv<para.svmaxtol
%	    break;
%    end
    %Save only every 10 sweeps to save time. (compared to save every step)
    if mod(loop,10)==0
    save(para.filename,'para','Vmat','mps','results','op');
	end

	fprintf('\nOccupation:\n');
	fprintf('%s',mat2str(getObservable({'occupation'},mps,Vmat,para),4));
	fprintf('\n');

    loop=loop+1;
end
end
