function [Fn1]=ComputeF(F,i,j,beta,m,x,rho,nabla_W0)

F(1,beta,i)=F(1,beta,i)+(m/rho(1,j)*(x(1,1,j)-x(1,1,i))*nabla_W0(1,beta,i,j));  
F(2,beta,i)=F(2,beta,i)+(m/rho(1,j)*(x(1,2,j)-x(1,2,i))*nabla_W0(1,beta,i,j)); 

Fn1=F;