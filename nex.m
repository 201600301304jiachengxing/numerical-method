clear,clc;
%% program 1.1
% ax^2+bx+c=0,p=[a,b,c],delta=b^2-4ac
p=[1,100000,1];
a=p(1);b=p(2);c=p(3);
if b.^2-4*a*c<0
    disp('无实根')
elseif b.^2-4*a*c==0
    x1=-b/(2*a);
else
    d=sqrt(b.^2-4*a*c);
    % error
    e=1e-4;
    if (abs(b)-d)/d>e
        x1=(-b+d)/(2*a);
        x2=(-b-d)/(2*a);
        ans=0
    elseif (abs(b)-d)/d<=e
        if b>0
            x1=(-2*c)/(b+d);
            x2=(-b-d)/(2*a);
            ans=1
        elseif b<=0
            x1=(-b+d)/(2*a);
            x2=(-2*c)/(b-d);
            ans=2
        end
    end
end
format long
x1
x2
(-b+d)/(2*a)
(-b-d)/(2*a)

%% program 1.2
%% (a) rn~r0/2^n,rn->0
clear;
a(1)=0.994;
r(1)=1-a(1);
for i=2:11
    a(i)=a(i-1)/2;
    r(i)=(1/2)^(i-1)-a(i);
end
plot([0:10],r,'-o')
format long
a;
r;

%% (b) rn->r
clear;
a(1)=1;
a(2)=0.497;
r(1)=1-a(1);
r(2)=1/2-a(2);
for i=3:11
    a(i)=(3/2)*a(i-1)-(1/2)*a(i-2);
    r(i)=(1/2)^(i-1)-a(i);
end
plot([0:10],r,'-o')
format long
a;
r

%% (c) r->2r 发散
clear;
a(1)=1;
a(2)=0.497;
r(1)=1-a(1);
r(2)=1/2-a(2);
for i=3:11
    a(i)=(5/2)*a(i-1)-a(i-2);
    r(i)=(1/2)^(i-1)-a(i);
end
plot([0:10],r,'-o')
format long
a;
r;

%% program 2.1
clear
newton(pi/2,1e-5);
clear
newton(5*pi,1e-5);
clear
newton(10*pi,1e-5);

%% program 2.2
%% binary,(0,1)
clear;
a=0;b=1;e=1e-4;
i=1;x(i)=(a+b)/2;fa=5*a-exp(a);fb=5*b-exp(b);f(i)=5*x(i)-exp(x(i));
while abs((b-a)/2)>e
    if f(i)*fa<0
        x(i+1)=(x(i)+a)/2;b=x(i);fb=f(i);
    elseif f(i)*fb<0
        x(i+1)=(x(i)+b)/2;a=x(i);fa=f(i);
    elseif f(i)==0
        break;
    end
    i=i+1;f(i)=5*x(i)-exp(x(i));
end
figure;plot([1:i],x(1:i),'-');xlabel('迭代次数');ylabel('x取值');title('迭代序列随次数变换');
figure;plot([1:i],f(1:i),'-');xlabel('迭代次数');ylabel('函数值');title('函数值随次数变化');

%% newton,初始值1
clear;
e=1e-4;
i=1;x(i)=1;f(i)=5*x(i)-exp(x(i));h(i)=5-exp(x(i));
i=2;x(i)=x(i-1)-f(i-1)/h(i-1);f(i)=5*x(i)-exp(x(i));h(i)=5-exp(x(i));
while abs(x(i)-x(i-1))>e
    i=i+1;x(i)=x(i-1)-f(i-1)/h(i-1);f(i)=5*x(i)-exp(x(i));h(i)=5-exp(x(i));
end
figure;plot([1:i],x(1:i),'-');xlabel('迭代次数');ylabel('x取值');title('迭代序列随次数变换');
figure;plot([1:i],f(1:i),'-');xlabel('迭代次数');ylabel('函数值');title('函数值随次数变化');

%% 割线法,初始值1,0.5
clear;
e=1e-4;
i=1;x(i)=1;f(i)=5*x(i)-exp(x(i));
i=2;x(i)=0.5;f(i)=5*x(i)-exp(x(i));
while abs(x(i)-x(i-1))>e
    i=i+1;x(i)=x(i-1)-f(i-1)*(x(i-1)-x(i-2))/(f(i-1)-f(i-2));f(i)=5*x(i)-exp(x(i));
end
figure;plot([1:i],x(1:i),'-');xlabel('迭代次数');ylabel('x取值');title('迭代序列随次数变换');
figure;plot([1:i],f(1:i),'-');xlabel('迭代次数');ylabel('函数值');title('函数值随次数变化');

%% 试位法,(0,1)
clear;
i=1;a=0;b=1;e=1e-4;fa=5*a-exp(a);fb=5*b-exp(b);c=b-fb*(b-a)/(fb-fa);fc=5*c-exp(c);x(i)=c;f(i)=fc;
while abs((b-a)/2)>e
    if fc*fa<0
        b=c;fb=5*b-exp(b);
    elseif fc*fb<0
        a=c;fa=5*a-exp(a);
    elseif fc==0
        break;
    end
    c=b-fb*(b-a)/(fb-fa);fc=5*c-exp(c);i=i+1;x(i)=c;f(i)=fc;
end
figure;plot([1:i],x(1:i),'-');xlabel('迭代次数');ylabel('x取值');title('迭代序列随次数变换');
figure;plot([1:i],f(1:i),'-');xlabel('迭代次数');ylabel('函数值');title('函数值随次数变化');

%% program 3.1
%% LU分解,LUx=b,Ly=b,Ux=y
clear;
A=[4,-1,1;4,-8,1;-2,1,5];
b=[7;-21;15];
n=length(A);
L=[1,0,0;0,1,0;0,0,1];
U=A;
for i=1:n
    L(i:n,i)=U(i:n,i)./U(i,i);
    U(i+1:n,:)=U(i+1:n,:)-(U(i+1:n,i)./U(i,i)).*U(i,:);
end
L0=[L,b];
for i=1:n
    L0(i+1:n,:)=L0(i+1:n,:)-(L0(i+1:n,i)./L0(i,i)).*L0(i,:);
end
U0=[U,L0(:,n+1)];
for i=n:-1:1
    U0(1:i-1,:)=U0(1:i-1,:)-(U0(1:i-1,i)./U0(i,i)).*U0(i,:);
    U0(i,:)=U0(i,:)/U(i,i);
end
x=U0(:,4);

%% jacobi A=D-L-U,Dx=(L+U)x+b
clear;
num=100;
A=[4,-1,1;4,-8,1;-2,1,5];
b=[7;-21;15];
n=length(A);
D=diag(diag(A));
invD=diag(1./diag(D));
LaddU=D-A;
x(:,1)=rands(3,1);
for i=2:num
    x(:,i)=invD*LaddU*x(:,i-1)+invD*b;
end
x
x=x(:,num);
    
%% guass A=D-L-U,(D-L)x=Ux+b
clear;
num=100;
A=[4,-1,1;4,-8,1;-2,1,5];
b=[7;-21;15];
n=length(A);
G=zeros(n,n);
for i=1:n
    G(i,1:i)=A(i,1:i);
end
U=G-A;
invG=inv(G);
x(:,1)=rands(3,1);
for i=2:num
    x(:,i)=invG*U*x(:,i-1)+invG*b;
end
x=x(:,num);

%% program 3.2
%% p118 3.a
clear;
num=1000;
n=50;
A=diag(ones(n,1)*4)+diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
b=ones(n,1)*3;
X=rands(n,1);
for i=1:num
    X(1,1)=(b(1,1)-A(1,2)*X(2,1))/A(1,1);
    for j=2:n-1
        X(j,1)=(b(j,1)-A(j,j-1)*X(j-1,1)-A(j,j+1)*X(j+1,1))/A(j,j);
    end
    X(n,1)=(b(n,1)-A(n,n-1)*X(n-1,1))/A(n,n);
end
a0=A;
b0=b;
x0=X
%% p118 3.b
clear
num=1000;
n=50;
A=diag(ones(n,1)*4)+diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
b=ones(n,1)*1.5+(ones(n,1).*(mod([1:n]',2)==0)-0.5);
X=rands(n,1);
for i=1:num
    X(1,1)=(b(1,1)-A(1,2)*X(2,1))/A(1,1);
    for j=2:n-1
        X(j,1)=(b(j,1)-A(j,j-1)*X(j-1,1)-A(j,j+1)*X(j+1,1))/A(j,j);
    end
    X(n,1)=(b(n,1)-A(n,n-1)*X(n-1,1))/A(n,n);
end
a0=A;
b0=b;
x0=X
%% p118 4
clear
num=2000;
n=50;
A=diag(ones(n,1)*12)+diag(ones(n-1,1)*(-2),1)+diag(ones(n-1,1)*(-2),-1)+diag(ones(n-2,1),2)+diag(ones(n-2,1),-2);
b=ones(n,1)*5;
X=rands(n,1);
for i=1:num
    X(1,1)=(b(1,1)-A(1,2)*X(2,1)-A(1,3)*X(3,1))/A(1,1);
    X(2,1)=(b(2,1)-A(2,1)*X(1,1)-A(2,3)*X(3,1)-A(2,4)*X(4,1))/A(2,2);
    for j=3:n-2
        X(j,1)=(b(j,1)-A(j,j-1)*X(j-1,1)-A(j,j+1)*X(j+1,1)-A(j,j-2)*X(j-2,1)-A(j,j+2)*X(j+2,1))/A(j,j);
    end
    X(n-1,1)=(b(n-1,1)-A(n-1,n-3)*X(n-3,1)-A(n-1,n-2)*X(n-2,1)-A(n-1,n)*X(n,1))/A(n-1,n-1);
    X(n,1)=(b(n,1)-A(n,n-1)*X(n-1,1)-A(n,n-2)*X(n-2,1))/A(n,n);
end
a0=A;
b0=b;
x0=X

%% program 4.1
clear,clc
n=11;m=100;left=-5;right=5;f=@(x)1./(1+x.*x);
[x0,x1,x2,y0,y1,y2,e]=lagrangePlot(left,right,f,n,m);
vpa(e)
dx=(right-left)/m;
dt=(right-left)/(n-1);
xd(1)=left;
yd(1)=f(xd(1));
for i=2:m+1
    xd(i)=xd(i-1)+dx;
    t=ceil((xd(i)-left)/dt);
    yd(i)=y0(t)*(xd(i)-x0(t+1))/(x0(t)-x0(t+1))+y0(t+1)*(xd(i)-x0(t))/(x0(t+1)-x0(t));
end
syms x;
int0=int(f,x,left,right);
comp=sum(yd*(right-left)/(m-1));
e1=int0-comp;

figure
plot(x0,y0,'go',x1,y1,'b-',x2,y2,'r-');
legend('data','function','lagrange10');

figure
plot(x0,y0,'go',x1,y1,'b-',xd,yd,'r-')
legend('data','function','lagrange2');

%% %%fu2
clear,clc
n=11;m=100;left=-5;right=5;f=@(x)1./(1+x.*x);
[x0,x1,x2,y0,y1,y2,e]=lagrangePlot(left,right,f,n,m);

syms x;
h=@(x)prod(ones(n,1)*x-x0');
f=1./(1+x.^2);
g=diff(f,x,11);
gmax=0;
hmax=0;
for i=-5:0.01:5
    y=subs(g,x,i);
    z=abs(h(i));
    if gmax<y
        gmax=y;
    end
    if hmax<z
        hmax=z;
    end
end
(vpa(gmax)*hmax)/prod([1:11])

%% %%fu

function []=newton(k,e)
f0=@(x)1/2+(1/4)*x.^2-x*sin(x)-(1/2)*cos(2*x);h0=@(x)(1/2)*x-sin(x)-x*cos(x)+sin(2*x);
i=1;x(i)=k;f(i)=f0(x(i));h(i)=h0(x(i));i=i+1;x(i)=x(i-1)-f(i-1)/h(i-1);f(i)=f0(x(i));h(i)=h0(x(i));
while abs(x(i)-x(i-1))>e&i<=100
    i=i+1;x(i)=x(i-1)-f(i-1)/h(i-1);f(i)=f0(x(i));h(i)=h0(x(i));
end
figure;plot([1:i],x(1:i),'-');xlabel('迭代次数');ylabel('x取值');title('迭代序列随次数变换');
figure;plot([1:i],f(1:i),'-');xlabel('迭代次数');ylabel('函数值');title('函数值随次数变化');
end







