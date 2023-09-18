function [flagg,proj] = Old_Proj_B1B2(z,t)
% Gradient Projection Method
% P1_Omega attempts to solve the projection problem:
% min (||x-z||_2)^2
%  x
% subject to: x∈Omega_1:={||x||_1≤t and ||x||_2≤1}
% Input:
%  z: point to be projected,column vector;
%  t: normalized parameter of 1-norm,constant;
% Output:
%  flagg: classifier,the value of falgg correspond to four cases in theorem 4.1(Liu2020);
%  proj: projected point;
N=size(z,1);
T=t^2;
vmax=max(abs(z));
flagg=0;
eps=1e-5;%allowable error
if((norm(z,2)<1 || abs(norm(z,2)-1)<=eps) && (norm(z,1)<t || abs(norm(z,1)-t)<=eps))
    w=z;
    flagg=1;
elseif(norm(z,2)>1 && (norm(z,1)<t*norm(z,2) || abs(norm(z,1)-t*norm(z,2))<=eps))
    w=z./norm(z,2);
    flagg=2;
elseif(norm(z,1)>t && norm(z,1)>t*norm(z,2))
    x=eplb(z,N,t,vmax/2);%
    if(norm(x,2)<1 || abs(norm(x,2)-1)<=eps)
        w=x;
        flagg=3;
    else
        [x1, lambda,~,~]=FindRoot_QASB(z,t);%
        w=x1;
        root=lambda;
        fprintf("Proj_B1B2 root is %f\n",root);
        flagg=4;
    end
else
    error('The input does not meet the requirements');
end
proj=w;
end