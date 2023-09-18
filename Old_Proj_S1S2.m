function [flagg,proj] = Old_Proj_S1S2(z,t)
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
if(I1>T)
    %pass;
    d=z;
    sf=sign(d);
    cen=(t/I1).*ones(I1,1);
    cbound=zeros(I1,1);
    cbound(1)=t;
    h=cbound-cen;

    fenzi=1-norm(cen,2)^2;
    fenmu=norm(h,2)^2;
    ss=sqrt(fenzi/fenmu);
    cqiu=cen+ss*h;
    shz=cqiu;
    wz=zeros(N,1);
    k=1;
    for ii=1:N
        if(abs(d(ii))>=vmax)
            wz(ii)=shz(k);
            k=k+1;
        end
    end
    w=wz.*sf;
    flagg=1;
end
if(abs(I1-T)<=eps)
    d=z;
    sf=sign(d);
    for j=1:N
        if(abs(d(j))>=vmax)
            e(j)=1/sqrt(I1);
        else
            e(j)=0;
        end
    end
    w=e'.*sf;
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