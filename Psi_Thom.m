function [Psi_lambda,Psi1_lambda,dd,ss,qq,finishedd]=Psi_Thom(y,lambda,t)
% This function is aim to solve interrelated variable of auxiliary function Psi with algorithm 2(Thom2015)
% Input:
%  y: point to be projected,column vector;
%  lambda: the variable where function Psi should be evaluated;
%  t: the rate of target 1-norm and target 2-norm;
% Output:
%  Psi_lambda: the value of function Psi in lambda;
%  Psi1_lambda: derivative of function Psi in lambda;
%  dd ss qq: three numbers needed to exactly compute the root of function Psi;
%  finishedd: indicater whether the correct interval has been found;
n=size(y,1);
finished=0;
s=0;
q=0;
d=0;
aj=0;t_aj=-lambda;ak=100000;t_ak=1000000;
for i=1:n
    temp=abs(y(i))-lambda;
    %if(temp>=0)
    if(temp>0)
        s=s+abs(y(i));
        q=q+y(i)*y(i);
        d=d+1;
        if(temp<t_ak)
            ak=abs(y(i));
            t_ak=temp;
        end
    else
        if(temp>t_aj)
            aj=abs(y(i));
            t_aj=temp;
        end
    end
end
if ((((s - d * aj) - t * sqrt(q - 2 * aj*s + d * aj*aj)) >= 0) && (((s - d * ak) - t * sqrt(q - 2 * ak*s + d * ak*ak)) < 0))
    finished = 1;
end
s_alp = s - d * lambda;
q_alp = q - 2 * lambda*s + d * lambda*lambda;
Psi_lam = s_alp/sqrt(q_alp)-t;
Psi1_lam = ((s_alp*s_alp)/q_alp-d)/sqrt(q_alp);
Psi_lambda=Psi_lam;
Psi1_lambda=Psi1_lam;
dd=d;
ss=s;
qq=q;
finishedd=finished;
end
    
    
        