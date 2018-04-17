function [bpx,bmx,nx,bpy,bmy,ny]=paritykron(bp,bosonparity)
dim=size(bp,1);
assert(mod(dim,2)==0);

oo=zeros(dim/2);
oe=bp(1:dim/2,1+dim/2:end);
eo=bp(1+dim/2:end,1:dim/2);
ee=zeros(dim/2);
id=eye(dim/2);

switch bosonparity
    case 'xy'
        ze=zeros(dim*dim/4);
        bpx=[kron(oo,id) ze ze kron(oe,id);  %This needs some time to get.
             ze kron(oo,id) kron(eo,id) ze;
             ze kron(oe,id) kron(oo,id) ze;
             kron(eo,id) ze ze kron(ee,id)];
        bpy=[kron(id,ee) ze kron(id,eo) ze;
             ze kron(id,oo) ze kron(id,oe);
             kron(id,oe) ze kron(id,oo) ze;
             ze kron(id,eo) ze kron(id,ee)];
    case 'x'
        idfull=eye(dim);
        ze=zeros(dim*dim/2);
        bpx=[kron(oo,idfull) kron(oe,idfull);
             kron(eo,idfull) kron(ee,idfull)];
        bpy=[kron(id,bp) ze;
             ze kron(id,bp)];
    case 'y'
        idfull=eye(dim);
        ze=zeros(dim*dim/2);
        bpx=[kron(bp,id) ze;
             ze kron(bp,id);];
        bpy=[kron(idfull,oo) kron(idfull,oe);
            kron(idfull,eo) kron(idfull,ee)];
end


 bmx=bpx';
 bmy=bpy';
 nx=bpx*bmx;
 ny=bpy*bmy;
end