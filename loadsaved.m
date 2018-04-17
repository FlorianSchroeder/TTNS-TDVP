function [Vmat,mps,loop,para,results,op]=loadsaved(para)
    r=load(para.filename);
    Vmat=r.Vmat;
    mps=r.mps;
    loop=r.para.loop+1;
    results=r.results;
    para.dk=r.para.dk;
    para.d_opt=r.para.d_opt;
    para.D=r.para.D;
   % para.loopmax=r.para.loopmax;
    para.shift=r.para.shift;
    para.hasexpanded=r.para.hasexpanded;
    para.adjust=r.para.adjust;%Initialize as 0, no need the edit.
    para.trustsite=r.para.trustsite;
  % para.trustsiteconverged=r.para.trustsiteconverged; %Initialize as 0, no need the edit.
    op=r.op;
end
