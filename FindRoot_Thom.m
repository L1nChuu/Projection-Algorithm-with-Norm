function [proj,ym1,ym2,root,lambda,psi_fun,psi1_fun,bi_steps,nw_steps]=FindRoot_Thom(y,l1,l2)
% This function is aim to solve the root of auxiliary function Psi and the projection of y with algorithm 3(Thom2015)
% Input:
%  y: point to be projected,column vector;
%  l1: target 1-norm;
%  l2: target 2-norm;
% Output:
%  proj: projected point;
%  ym1: the maximum component of projected point;
%  ym2: the second largest component of projected point;
%  root: the root of function Psi;
%  lambda: the independent variable of function Psi in every iteration,colunm vector;
%  psi_fun: a vector composed of values of function Psi on lambda;
%  psi1_fun: a vector composed of values of gradient of function Psi on lambda;
%  bi_steps: iteration steps with Bisection solver;
%  nw_steps: iteration steps with Newton solver;  
n=size(y,1);
t=l1/l2;
bi_iter=0;
nw_iter=0;
eps=1e-6;
i=1;
lambda(i)=0;
yyuan=y;
y=abs(y);
% lambda(i)=vpa(lambda(i));
lambda(i)=lambda(i);
[psi_lambda,psi1_lambda,dd,ss,qq,finishedd]=Psi_Thom(y,lambda(i),t);
psi_fun(i)=psi_lambda;
psi1_fun(i)=psi1_lambda;
if(psi_lambda<=0)
    root=1.0 / dd * (ss - l1 * sqrt(dd*qq - ss * ss) / sqrt(dd*l2*l2 - l1 * l1));
    lambda(i+1)=root;
    [psi_lambda,psi1_lambda,dd,ss,qq,finishedd]=Psi_Thom(y,root,t);
    psi_fun(i+1)=psi_lambda;
    psi1_fun(i+1)=psi1_lambda;
    bi_steps=0;
    nw_steps=0;
    ym1=0;
    ym2=0;
else
    ym1=abs(y(1));
    ym2=0;
    for j=1:n
        if(ym1<abs(y(j)))
            ym1=abs(y(j));
        end
    end
    
    for k=1:n
        if(abs(y(k))~=ym1)
            ym2=abs(y(k));
        end
        for j=1:n
            if(abs(y(j))~=ym1&&ym2<abs(y(j)))
                ym2=abs(y(j));
            end
        end
    end
    lk=0;rk=ym2;
    
    i=i+1;
    lambda(i)=(lk+rk)/2;
    bi_iter=bi_iter+1;
    [psi_lambda,psi1_lambda,dd,ss,qq,finishedd]=Psi_Thom(y,lambda(i),t);
    while (finishedd==0)
        if(psi_lambda>=0)
            lk=lambda(i);
            psi_fun(i)=psi_lambda;
            psi1_fun(i)=psi1_lambda;            
        else
            rk=lambda(i);
            psi_fun(i)=psi_lambda;
            psi1_fun(i)=psi1_lambda;
        end
        i=i+1;
        %digits(12);
%         digits(13);
%         lambda(i)=vpa(lambda(i-1))-psi_lambda/psi1_lambda;
        lambda(i)=lambda(i-1)-psi_lambda/psi1_lambda;
        if((lambda(i)>rk)||(lambda(i)<lk))
            lambda(i)=(lk+rk)/2;
            [psi_lambda,psi1_lambda,dd,ss,qq,finishedd]=Psi_Thom(y,lambda(i),t);
            psi_fun(i)=psi_lambda;
            psi1_fun(i)=psi1_lambda;
            bi_iter=bi_iter+1;
        elseif(i>=3 && (abs(lambda(i)-lambda(i-2))<=eps))
            lambda(i)=(lk+rk)/2;
            [psi_lambda,psi1_lambda,dd,ss,qq,finishedd]=Psi_Thom(y,lambda(i),t);
            psi_fun(i)=psi_lambda;
            psi1_fun(i)=psi1_lambda;
            bi_iter=bi_iter+1;
        else
            [psi_lambda,psi1_lambda,dd,ss,qq,finishedd]=Psi_Thom(y,lambda(i),t);
            nw_iter=nw_iter+1;
            psi_fun(i)=psi_lambda;
            psi1_fun(i)=psi1_lambda;            
        end
    end
end
    root=1.0 / dd * (ss - l1 * sqrt(dd*qq - ss * ss) / sqrt(dd*l2*l2 - l1 * l1));
    lambda(i+1)=root;
    [psi_lambda,psi1_lambda,dd,ss,qq,finishedd]=Psi_Thom(y,root,t);
    psi_fun(i+1)=psi_lambda;
    psi1_fun(i+1)=psi1_lambda;
    bi_steps=bi_iter;
    nw_steps=nw_iter;
% Compute result of the projection in-place
rho=0;
for i=1:n
    if (yyuan(i)>=0)
        t=yyuan(i)-root;
        if(t>0)
            yyuan(i)=t;
            rho=rho+t^2;
        else
            yyuan(i)=0;
        end
    end
    if (yyuan(i)<0)
        t=-yyuan(i)-root;
        if (t>0)
            yyuan(i)=-t;
            rho=rho+t^2;
        else
            yyuan(i)=0;
        end
    end
end
for i=1:n
    proj(i)=l2/sqrt(rho)*yyuan(i);
end
proj=proj';
   


    