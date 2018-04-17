function [sigmaX,sigmaY,sigmaZ]=spinop(base)
switch base
    case 'X' %\sigma_x eigenstates as bases
        sigmaX=[-1 0;0 1]; %In the parity basis
        sigmaY=[0 1i;-1i 0];
        sigmaZ=[0 -1;-1 0];
    case 'Y'
        sigmaX=[0 1i;-1i 0]; %In the parity basis, the bases are simgaY eigenstates
        sigmaY=[-1 0;0 1];
        sigmaZ=[0 1;1 0];
    case 'Z'
        sigmaX=[0 1;1 0]; %In the non parity basis
        sigmaY=[0 -1i;1i 0];
        sigmaZ=[1 0;0 -1];
end

end