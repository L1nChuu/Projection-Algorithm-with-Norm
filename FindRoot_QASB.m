function [x, lambda, iter_step, flag]=FindRoot_QASB(v,z)
% Quadratic Approximation Secant Bisection Method
% This function corresponds to the algorithm 1(Liu2020);
% Input:
%  v: point to be projected,column vector;
%  z: normalized parameter of 1-norm,constant;
% Output:
%  x: projected point;
%  lambda: the root of function phi;
%  iter_step: iteration step;
%  flag: classifier;
n=size(v,1);
delta=1e-9;
flag=0;
iter_step=0;
temp1=0;
z2=z*z;
if(z<0)
    return;
end
s=norm(v,1);
q=norm(v,2);
r=max(abs(v));
i=n-1;
if(s*s-z2*q<0)
    lambda=(s-z*sqrt(((i+1)*q-s*s)/(i+1-z*z)))/(i+1);
    flag=0;
else
    V_i=0;  rho_1=0; rho_2=0; rho_Q=0; rho_S=0; 
    s_1=0; s_2=0; s_Q=0; s_S=0; q_1=0; q_2=0; q_Q=0; q_S=0; s=0; q=0;
    lambda_1=0; lambda_2=r;
    for i=0:n-1
        if(abs(v(i+1))>=lambda_2)
            rho_2=rho_2+1;s_2=s_2+abs(v(i+1));
            q_2=q_2+v(i+1)*v(i+1);
            rho_1=rho_1+1;
            s_1=s_1+abs(v(i+1));
            q_1=q_1+v(i+1)*v(i+1);
            s=s+abs(v(i+1));
            q=q+v(i+1)*v(i+1);
        end
        if(abs(v(i+1))>=lambda_1 && abs(v(i+1))<lambda_2)
            x(V_i+1)=abs(v(i+1));s_1=s_1+x(V_i+1);q_1=q_1+x(V_i+1)*x(V_i+1);rho_1=rho_1+1;
            V_i=V_i+1;
        end
    end
    z2=z*z;
    f_lambda_1=(rho_1-z2)*(rho_1*lambda_1-2*s_1)*lambda_1+s_1*s_1-z2*q_1;
    f_lambda_2=(rho_2-z2)*(rho_2*lambda_2-2*s_2)*lambda_2+s_2*s_2-z2*q_2;

    V_i_b=0; V_i_e=V_i-1;    
    
    while(flag==0)
        iter_step=iter_step+1;
        temp1=sqrt((rho_1*q_1-s_1*s_1)/(rho_1-z2));
        lambda_Q=(s_1-z*temp1)/rho_1;
        
        lambda_S=lambda_1-f_lambda_1*(lambda_2-lambda_1)/(f_lambda_2-f_lambda_1);
        if(abs(lambda_Q-lambda_S)<=delta)
            lambda=lambda_Q;flag=1;
            break;
        end
        lambda=(lambda_Q+lambda_S)/2;
        s_Q=0;s_S=0;s=0;q_Q=0;q_S=0;q=0;rho_Q=0;rho_S=0;rho=0;
        i=V_i_b; j=V_i_e;
        while(i<=j)
            while((i<=V_i_e) && (x(i+1)<=lambda))
                if(x(i+1)>lambda_Q)
                    s_Q=s_Q+x(i+1);q_Q=q_Q+x(i+1)*x(i+1);rho_Q=rho_Q+1;
                end
                i=i+1;
            end
            while((j>=V_i_b) && (x(j+1)>lambda))
                if(x(j+1)>lambda_S)
                    s_S=s_S+x(j+1);q_S=q_S+x(j+1)*x(j+1);rho_S=rho_S+1;
                else
                    s=s+x(j+1);q=q+x(j+1)*x(j+1);rho=rho+1;
                end
                j=j-1;
            end
            if(i<j)
                if(x(i+1)>lambda_S)
                    s_S=s_S+x(i+1);q_S=q_S+x(i+1)*x(i+1);rho_S=rho_S+1;
                else
                    s=s+x(i+1);q=q+x(i+1)*x(i+1);rho=rho+1;
                end
                if(x(j+1)>lambda_Q)
                    s_Q=s_Q+x(j+1);q_Q=q_Q+x(j+1)*x(j+1);rho_Q=rho_Q+1;
                end
                temp=x(i+1);x(i+1)=x(j+1);x(j+1)=temp;
                i=i+1;j=j-1;
            end
        end
        s_S=s_S+s_2; q_S=q_S+q_2; rho_S=rho_S+rho_2;
        s=s+s_S; q=q+q_S; rho=rho+rho_S;
        s_Q=s_Q+s; q_Q=q_Q+q; rho_Q=rho_Q+rho;
        f_lambda_S=(rho_S-z2)*(rho_S*lambda_S-2*s_S)*lambda_S+s_S*s_S-z2*q_S;
        f_lambda=(rho-z2)*(rho*lambda-2*s)*lambda+s*s-z2*q;
        f_lambda_Q=(rho_Q-z2)*(rho_Q*lambda_Q-2*s_Q)*lambda_Q+s_Q*s_Q-z2*q_Q;   
        if(abs(f_lambda)<delta)
            flag=1;
            break;
        end
        if(abs(f_lambda_Q)<delta)
            flag=2;
            lambda=lambda_Q;
            break;
        end
        if(f_lambda<0)
            lambda_2=lambda;  s_2=s; q_2=q; rho_2=rho;
            f_lambda_2=f_lambda;            
            
            lambda_1=lambda_Q; s_1=s_Q; q_1=q_Q; rho_1=rho_Q;
            f_lambda_1=f_lambda_Q;
            V_i_e=j;  i=V_i_b;
            while (i <= j)
                while( (i <= V_i_e) && (x(i+1) <= lambda_Q) )
                    i=i+1;
                end
                while( (j>=V_i_b) && (x(j+1) > lambda_Q) )
                    j=j-1;
                end
                if (i<j)                    
                    x(j+1)=x(i+1);
                    i=i+1;   j=j-1;
                end
            end
            V_i_b=i; V_i=V_i_e-V_i_b+1;
        else
            lambda_1=lambda;  s_1=s;  q_1=q; rho_1=rho;
            f_lambda_1=f_lambda;
            
            lambda_2=lambda_S; s_2=s_S; q_2=q_S; rho_2=rho_S;
            f_lambda_2=f_lambda_S;
            
            V_i_b=i;  j=V_i_e;
            while (i <= j)
                while( (i <= V_i_e) && (x(i+1) <= lambda_S) )
                    i=i+1;
                end
                while( (j>=V_i_b) && (x(j+1) > lambda_S) )
                    j=j-1;
                end
                if (i<j)
                    x(i+1)=x(j+1);
                    i=i+1;   j=j-1;
                end
            end
            V_i_e=j; V_i=V_i_e-V_i_b+1;
        end     
    end
end
norm2x=0;
for i=1:n
    if(v(i)>=0)
        if(v(i)>lambda)
            x(i)=v(i)-lambda;
            norm2x=norm2x+x(i)*x(i);
        else
            x(i)=0;
        end
    end
    if(v(i)<0)
        if(-v(i)>lambda)
            x(i)=v(i)+lambda;
            norm2x=norm2x+x(i)*x(i);
        else
            x(i)=0;
        end
    end
end
for i=1:n
    x(i)=x(i)/sqrt(norm2x);
end
x=x';
return;

            
        
                    
            
                
                    
        

    
