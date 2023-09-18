function [flagg,proj] = Proj_B1S2(z,t)
% Gradient Projection Method
% P3_Omega attempts to solve the projection problem:
% min (||x-z||_2)^2
%  x
% subject to: x∈Omega_3:={||x||_1≤t and ||x||_2=1}
% Input:
%  z: point to be projected,column vector;
%  t: normalized parameter of 1-norm,constant;
% Output:
%  flagg: classifier,the value of falgg correspond to three cases in theorem 4.3(Liu2020);
%  proj: projected point;
N=size(z,1);
T=t^2;
vmax=max(abs(z));
I1=sum(abs(z)>=vmax);
flagg=0;
eps=1e-5;%allowable error
if((I1<T) || abs(I1-T)<=eps)
    if(norm(z,1)>t*norm(z,2))
        [x1, lambda, ~, ~]=FindRoot_QASB(z,t);
        w=x1; 
        root=lambda;
        fprintf("Proj_B1S2 root is %f\n",root);
        flagg=1;
    else
        w=z./norm(z,2);
        flagg=2;
    end
else
    d=z;
    sf=sign(d);
    e_I1=zeros(N,1);
    e_I1(find(abs(d)>=vmax))=1;
    e_1=[1;zeros(N-1,1)];
    beta=((I1-T)/(I1-1))^(1/2);
    alpha=(t-beta)/I1;
    wz=alpha*e_I1+beta*e_1;
    w=wz.*sf;
    flagg=3;
end
proj=w;
end
