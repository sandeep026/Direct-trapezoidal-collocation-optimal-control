clear
close all

import casadi.*

T=2;
a=1;
b=-2.694528;
c=-1.155356;
N=200;
opti=Opti();
X=opti.variable(2,N+1);
U=opti.variable(1,N);

X0=[0;0];

opti.subject_to(X(:,1)-X0==0)

dt = T/N; % length of a control interval

for k=1:N % loop over control intervals
opti.subject_to(X(:,k+1)-X(:,k)-0.5*dt*(ode_fun(X(:,k),U(:,k))+ode_fun(X(:,k+1),U(:,k)))==0)
end

s=0;
for k=1:N
s=s+0.5*dt*(U(:,k)'*U(:,k)+U(:,k)'*U(:,k))
end

obj=s;

opti.minimize(obj)
opti.subject_to(a*X(1,end)+b*X(2,end)-c==0)
opti.solver('ipopt')
sol=opti.solve();

% Compare the numerical and analytical solution
t=linspace(0,2,N+1);
X=sol.value(X);
U=sol.value(U);
[Xa,Ua]=analytical_solution(t);
plot(Xa(1,:),Xa(2,:))
hold on
scatter(X(1,:),X(2,:))
grid
xlabel('$x_1$')
ylabel('$x_2$')
legend('analytical','single shooting','Location','northwest')


 
 
figure
plot(t,Ua)
hold on
stairs(t,[U U(:,end)])
xlabel('time [s]')
ylabel('f')
legend('analytical','single shooting','Location','northwest')

figure
Lag=opti.f+ opti.lam_g'*opti.g;
H=hessian(Lag,opti.x);
spy(H) 
legend('Hessian sparsity') 
figure 
Jac=jacobian(opti.g,opti.x);
spy(Jac)
legend('Jacobian sparsity') 
