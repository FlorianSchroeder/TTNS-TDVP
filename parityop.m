function pa_op=parityop(para)
switch para.model
    case 'SpinDoubleBoson'
        pa_op=cell(1,para.L);
        [sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
        switch para.bosonparity
            case 'xy'
                pa_op{1}=sigmaZ;
            case 'x'
                pa_op{1}=sigmaY;
            case 'y'
                pa_op{1}=sigmaX;
        end
        for j=2:para.L
            dim=para.dk(j);
            [bp,bm,n]=bosonop(sqrt(para.dk(j)),para.shift(j),para.parity);
            if para.parity=='n'
                idm=eye(size(n));
                bpx=kron(bp,idm);bmx=bpx';nx=kron(n,idm);
                bpy=kron(idm,bp);bmy=bpy';ny=kron(idm,n);
            else
                [bpx,bmx,nx,bpy,bmy,ny]=paritykron(bp,para.bosonparity);
            end
            switch para.bosonparity
            case 'xy'
                n=nx+ny;
            case 'x'
                n=nx;
            case 'y'
                n=ny;
            end
            paop=spdiags((-1).^diag(n),0,dim,dim);
            pa_op{j}=real(paop);
        end
    case 'SpinBoson'
        pa_op=cell(1,para.L);
        [sigmaX,sigmaY,sigmaZ]=spinop(para.spinbase);
        pa_op{1}=sigmaX;
        for j=2:para.L
            dim=para.dk(j);
            paop=sparse(dim,dim);
            if para.parity~='n'
                for k=1:dim/2
                    paop(k,k)=-1;
                    paop(k+dim/2,k+dim/2)=1;
                end
            else
                paop(end:end)=1;
                for k=para.dk(j)-1:-1:1
                    paop(k,k)=paop(k+1,k+1)*-1;
                end
            end
            pa_op{j}=paop;
        end
end
end