function [flagg,proj] = Proj_S1S2(z,t)
% Gradient Projection Method
% P2_Omega attempts to solve the projection problem:
% min (||x-z||_2)^2
%  x
% subject to: x∈Omega_2:={||x||_1=t and ||x||_2=1}
% Input:
%  z: point to be projected,column vector;
%  t: normalized parameter of 1-norm,constant;
% Output:
%  flagg: C,the value of falgg correspond to three cases in theorem 4.2(Liu2020);
%  proj: projected point;
N=size(z,1);
T=t^2;
vmax=max(abs(z));
I1=sum(abs(z)>=vmax);
flagg=0;
eps=1e-5;%allowable error
d=z;
sf=sign(d);
e_I1=zeros(N,1);
e_I1(find(abs(d)>=vmax))=1;
if(I1>T)
    e_1=[1;zeros(N-1,1)];
    beta=((I1-T)/(I1-1))^(1/2);
    alpha=(t-beta)/I1;
    wz=alpha*e_I1+beta*e_1;
    w=wz.*sf;
    flagg=1;
end
if(abs(I1-T)<=eps)
    wz=e_I1*I1^(-1/2);
    w=wz.*sf;
    flagg=2;
end
if(I1<T)
    [x1, lambda,~,~]=FindRoot_QASB(z,t);
    w=x1;
    root=lambda;
    fprintf("Proj_S1S2 root is %f\n",root);
    flagg=3;
end
proj=w;
end