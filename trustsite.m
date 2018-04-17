function para=trustsite(para,results)
%Determine until which site in the Wilson chain the results are reliable
%The cretieria used here is the order of magnitude of the smallest singular value ssv should be
%smaller than the order of magnitude of the coupling/onsite energy.


% for k=2:para.L
%     sv=results.Vmat_sv{k}(end-1:end);%Two smallest singular values
%     ssv_V(k)=sv(2)^2/sv(1);
% end
% ssv_V=ssv_V(2:end);
%
% for k=1:para.L-1
%     sv=results.Amat_sv{k}(end-1:end);%Two smallest singular values
%     ssv_A(k)=sv(2)^2/sv(1);
% end
%
% %ssv_V=cellfun(@min,results.Vmat_sv(2:end));
% %ssv_A=cellfun(@min,results.Amat_sv(:))';
% ssv=max(ssv_V,ssv_A);
%
% sv_dig=round(log10(ssv));
% t_dig=round(log10(para.t)');
%
% site=min(find(sv_dig>t_dig))-1;

%%%%%%%%%%%%%
 site1=min(find(results.vNEdiff>0.01));
  if isempty(site1)
     site1=para.L;
 end

 if results.Eerror(end)<1e-4
	site2=max(find(para.chain{1}.epsilon>results.Eerror(end)));
 else
	 site2=1;
 end
 site=min(site1,site2);
 site=site-1;
 para.trustsite(para.loop)=site;
 para.precisesite=site2;
end
