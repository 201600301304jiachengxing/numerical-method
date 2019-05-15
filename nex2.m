%% 5.1.1
clear,clc;
% d=1/2*g*t^2,logd=2logt+log1/2*g
t1=[0.2,0.4,0.6,0.8,1.0];
d1=[0.1960,0.7835,1.7630,3.1345,4.8975];
x=t1.^2;y=d1;sumx2=(x)*(x)';sumxy=(x)*(y)';a=sumxy/sumx2;
g=2*a;
plot(t1,d1,'o','LineWidth',1);
hold on
t=[0:0.01:1];
plot(t,1/2*g*t.^2,'-','LineWidth',1);
xlabel('time t/s');ylabel('distance d/m');disp(g);
e=sum(abs(a*t1.^2-d1).^2)

%% 5.1.2
t2=[0.2,0.4,0.6,0.8,1.0];
d2=[0.1965,0.7855,1.7675,3.1420,4.9095];
x=t2.^2;
y=d2;
sumx2=(x)*(x)';
sumxy=(x)*(y)';
a=sumxy/sumx2;
g=2*a;
plot(t2,d2,'o','LineWidth',1);
hold on
t=[0:0.01:1];
plot(t,1/2*g*t.^2,'-','LineWidth',1);
xlabel('time t/s');
ylabel('distance d/m');
disp(g);e=sum(abs(a*t2.^2-d2).^2)

%% 5.2.1
clear,clc;
X=[0,1,2,3,4,5,6];
Y=[1,0,0,1,2,2,1];
dx0=-0.6;
dxn=-1.8;
N=length(X)-1;
H=diff(X);
D=diff(Y)./H;
A=H(2:N-1);
B=2*(H(1:N-1)+H(2:N));
C=H(2:N);
U=6*diff(D);
B(1)=B(1)-H(1)/2;
U(1)=U(1)-3*(D(1)-dx0);
B(N-1)=B(N-1)-H(N-1)/2;
U(N-1)=U(N-1)-3*(dxn-D(N));
for k=2:N-1
    temp=A(k-1)/B(k-1);
    B(k)=B(k)-temp*C(k-1);
    U(k)=U(k)-temp*U(k-1);
end
M(N)=U(N-1)/B(N-1);
for k=N-2:-1:1
    M(k+1)=(U(k)-C(k)*M(k+2))/B(k);
end
M(1)=3*(D(1)-dx0)/H(1)-M(2)/2;
M(N+1)=3*(dxn-D(N))/H(N)-M(N)/2;
for k=0:N-1
    S(k+1,1)=(M(k+2)-M(k+1))/(6*H(k+1));
    S(k+1,2)=M(k+1)/2;
    S(k+1,3)=D(k+1)-H(k+1)*(2*M(k+1)+M(k+2))/6;
    S(k+1,4)=Y(k+1);
end
plot(X,Y,'o')
hold on
k=1;
x1=[0:0.01:1];
y1=S(k,1)*(x1-X(k)).^3+S(k,2)*(x1-X(k)).^2+S(k,3)*(x1-X(k))+S(k,4);
k=2;
x2=[1:0.01:2];
y2=S(k,1)*(x2-X(k)).^3+S(k,2)*(x2-X(k)).^2+S(k,3)*(x2-X(k))+S(k,4);
k=3;
x3=[2:0.01:3];
y3=S(k,1)*(x3-X(k)).^3+S(k,2)*(x3-X(k)).^2+S(k,3)*(x3-X(k))+S(k,4);
k=4;
x4=[3:0.01:4];
y4=S(k,1)*(x4-X(k)).^3+S(k,2)*(x4-X(k)).^2+S(k,3)*(x4-X(k))+S(k,4);
k=5;
x5=[4:0.01:5];
y5=S(k,1)*(x5-X(k)).^3+S(k,2)*(x5-X(k)).^2+S(k,3)*(x5-X(k))+S(k,4);
k=6;
x6=[5:0.01:6];
y6=S(k,1)*(x6-X(k)).^3+S(k,2)*(x6-X(k)).^2+S(k,3)*(x6-X(k))+S(k,4);
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,'b-');
hold on

% nature
clear,clc;
X=[0,1,2,3,4,5,6];
Y=[1,0,0,1,2,2,1];
dx0=-0.6;
dxn=-1.8;
N=length(X)-1;
H=diff(X);
D=diff(Y)./H;
A=H(2:N-1);
B=2*(H(1:N-1)+H(2:N));
C=H(2:N);
U=6*diff(D);
B(1)=B(1)-H(1)/2;
U(1)=U(1)-3*(D(1)-dx0);
B(N-1)=B(N-1)-H(N-1)/2;
U(N-1)=U(N-1)-3*(dxn-D(N));
for k=2:N-1
    temp=A(k-1)/B(k-1);
    B(k)=B(k)-temp*C(k-1);
    U(k)=U(k)-temp*U(k-1);
end
M(N)=U(N-1)/B(N-1);
for k=N-2:-1:1
    M(k+1)=(U(k)-C(k)*M(k+2))/B(k);
end
M(1)=0;
M(N+1)=0;
for k=0:N-1
    S(k+1,1)=(M(k+2)-M(k+1))/(6*H(k+1));
    S(k+1,2)=M(k+1)/2;
    S(k+1,3)=D(k+1)-H(k+1)*(2*M(k+1)+M(k+2))/6;
    S(k+1,4)=Y(k+1);
end
k=1;
x1=[0:0.01:1];
y1=S(k,1)*(x1-X(k)).^3+S(k,2)*(x1-X(k)).^2+S(k,3)*(x1-X(k))+S(k,4);
k=2;
x2=[1:0.01:2];
y2=S(k,1)*(x2-X(k)).^3+S(k,2)*(x2-X(k)).^2+S(k,3)*(x2-X(k))+S(k,4);
k=3;
x3=[2:0.01:3];
y3=S(k,1)*(x3-X(k)).^3+S(k,2)*(x3-X(k)).^2+S(k,3)*(x3-X(k))+S(k,4);
k=4;
x4=[3:0.01:4];
y4=S(k,1)*(x4-X(k)).^3+S(k,2)*(x4-X(k)).^2+S(k,3)*(x4-X(k))+S(k,4);
k=5;
x5=[4:0.01:5];
y5=S(k,1)*(x5-X(k)).^3+S(k,2)*(x5-X(k)).^2+S(k,3)*(x5-X(k))+S(k,4);
k=6;
x6=[5:0.01:6];
y6=S(k,1)*(x6-X(k)).^3+S(k,2)*(x6-X(k)).^2+S(k,3)*(x6-X(k))+S(k,4);
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,'b-');
hold on

% out
clear,clc;
X=[0,1,2,3,4,5,6];
Y=[1,0,0,1,2,2,1];
dx0=-0.6;
dxn=-1.8;
N=length(X)-1;
H=diff(X);
D=diff(Y)./H;
A=H(2:N-1);
B=2*(H(1:N-1)+H(2:N));
C=H(2:N);
U=6*diff(D);
B(1)=B(1)-H(1)/2;
U(1)=U(1)-3*(D(1)-dx0);
B(N-1)=B(N-1)-H(N-1)/2;
U(N-1)=U(N-1)-3*(dxn-D(N));
for k=2:N-1
    temp=A(k-1)/B(k-1);
    B(k)=B(k)-temp*C(k-1);
    U(k)=U(k)-temp*U(k-1);
end
M(N)=U(N-1)/B(N-1);
for k=N-2:-1:1
    M(k+1)=(U(k)-C(k)*M(k+2))/B(k);
end
M(1)=M(2)-H(1)*(M(3)-M(2))/H(2);
M(N+1)=M(N)+H(N)*(M(N)-M(N-1))/H(N-1);
for k=0:N-1
    S(k+1,1)=(M(k+2)-M(k+1))/(6*H(k+1));
    S(k+1,2)=M(k+1)/2;
    S(k+1,3)=D(k+1)-H(k+1)*(2*M(k+1)+M(k+2))/6;
    S(k+1,4)=Y(k+1);
end
k=1;
x1=[0:0.01:1];
y1=S(k,1)*(x1-X(k)).^3+S(k,2)*(x1-X(k)).^2+S(k,3)*(x1-X(k))+S(k,4);
k=2;
x2=[1:0.01:2];
y2=S(k,1)*(x2-X(k)).^3+S(k,2)*(x2-X(k)).^2+S(k,3)*(x2-X(k))+S(k,4);
k=3;
x3=[2:0.01:3];
y3=S(k,1)*(x3-X(k)).^3+S(k,2)*(x3-X(k)).^2+S(k,3)*(x3-X(k))+S(k,4);
k=4;
x4=[3:0.01:4];
y4=S(k,1)*(x4-X(k)).^3+S(k,2)*(x4-X(k)).^2+S(k,3)*(x4-X(k))+S(k,4);
k=5;
x5=[4:0.01:5];
y5=S(k,1)*(x5-X(k)).^3+S(k,2)*(x5-X(k)).^2+S(k,3)*(x5-X(k))+S(k,4);
k=6;
x6=[5:0.01:6];
y6=S(k,1)*(x6-X(k)).^3+S(k,2)*(x6-X(k)).^2+S(k,3)*(x6-X(k))+S(k,4);
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,'b-');
hold on

% snn
clear,clc;
X=[0,1,2,3,4,5,6];
Y=[1,0,0,1,2,2,1];
dx0=-0.6;
dxn=-1.8;
N=length(X)-1;
H=diff(X);
D=diff(Y)./H;
A=H(2:N-1);
B=2*(H(1:N-1)+H(2:N));
C=H(2:N);
U=6*diff(D);
B(1)=B(1)-H(1)/2;
U(1)=U(1)-3*(D(1)-dx0);
B(N-1)=B(N-1)-H(N-1)/2;
U(N-1)=U(N-1)-3*(dxn-D(N));
for k=2:N-1
    temp=A(k-1)/B(k-1);
    B(k)=B(k)-temp*C(k-1);
    U(k)=U(k)-temp*U(k-1);
end
M(N)=U(N-1)/B(N-1);
for k=N-2:-1:1
    M(k+1)=(U(k)-C(k)*M(k+2))/B(k);
end
M(1)=M(2);
M(N+1)=M(N);
for k=0:N-1
    S(k+1,1)=(M(k+2)-M(k+1))/(6*H(k+1));
    S(k+1,2)=M(k+1)/2;
    S(k+1,3)=D(k+1)-H(k+1)*(2*M(k+1)+M(k+2))/6;
    S(k+1,4)=Y(k+1);
end
k=1;
x1=[0:0.01:1];
y1=S(k,1)*(x1-X(k)).^3+S(k,2)*(x1-X(k)).^2+S(k,3)*(x1-X(k))+S(k,4);
k=2;
x2=[1:0.01:2];
y2=S(k,1)*(x2-X(k)).^3+S(k,2)*(x2-X(k)).^2+S(k,3)*(x2-X(k))+S(k,4);
k=3;
x3=[2:0.01:3];
y3=S(k,1)*(x3-X(k)).^3+S(k,2)*(x3-X(k)).^2+S(k,3)*(x3-X(k))+S(k,4);
k=4;
x4=[3:0.01:4];
y4=S(k,1)*(x4-X(k)).^3+S(k,2)*(x4-X(k)).^2+S(k,3)*(x4-X(k))+S(k,4);
k=5;
x5=[4:0.01:5];
y5=S(k,1)*(x5-X(k)).^3+S(k,2)*(x5-X(k)).^2+S(k,3)*(x5-X(k))+S(k,4);
k=6;
x6=[5:0.01:6];
y6=S(k,1)*(x6-X(k)).^3+S(k,2)*(x6-X(k)).^2+S(k,3)*(x6-X(k))+S(k,4);
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,'b-');
hold on

% s2
clear,clc;
X=[0,1,2,3,4,5,6];
Y=[1,0,0,1,2,2,1];
dx0=-0.6;
dxn=-1.8;
N=length(X)-1;
H=diff(X);
D=diff(Y)./H;
A=H(2:N-1);
B=2*(H(1:N-1)+H(2:N));
C=H(2:N);
U=6*diff(D);
B(1)=B(1)-H(1)/2;
U(1)=U(1)-3*(D(1)-dx0);
B(N-1)=B(N-1)-H(N-1)/2;
U(N-1)=U(N-1)-3*(dxn-D(N));
for k=2:N-1
    temp=A(k-1)/B(k-1);
    B(k)=B(k)-temp*C(k-1);
    U(k)=U(k)-temp*U(k-1);
end
M(N)=U(N-1)/B(N-1);
for k=N-2:-1:1
    M(k+1)=(U(k)-C(k)*M(k+2))/B(k);
end
M(1)=1;
M(N+1)=-1;
for k=0:N-1
    S(k+1,1)=(M(k+2)-M(k+1))/(6*H(k+1));
    S(k+1,2)=M(k+1)/2;
    S(k+1,3)=D(k+1)-H(k+1)*(2*M(k+1)+M(k+2))/6;
    S(k+1,4)=Y(k+1);
end
k=1;
x1=[0:0.01:1];
y1=S(k,1)*(x1-X(k)).^3+S(k,2)*(x1-X(k)).^2+S(k,3)*(x1-X(k))+S(k,4);
k=2;
x2=[1:0.01:2];
y2=S(k,1)*(x2-X(k)).^3+S(k,2)*(x2-X(k)).^2+S(k,3)*(x2-X(k))+S(k,4);
k=3;
x3=[2:0.01:3];
y3=S(k,1)*(x3-X(k)).^3+S(k,2)*(x3-X(k)).^2+S(k,3)*(x3-X(k))+S(k,4);
k=4;
x4=[3:0.01:4];
y4=S(k,1)*(x4-X(k)).^3+S(k,2)*(x4-X(k)).^2+S(k,3)*(x4-X(k))+S(k,4);
k=5;
x5=[4:0.01:5];
y5=S(k,1)*(x5-X(k)).^3+S(k,2)*(x5-X(k)).^2+S(k,3)*(x5-X(k))+S(k,4);
k=6;
x6=[5:0.01:6];
y6=S(k,1)*(x6-X(k)).^3+S(k,2)*(x6-X(k)).^2+S(k,3)*(x6-X(k))+S(k,4);
plot(x1,y1,x2,y2,x3,y3,x4,y4,x5,y5,x6,y6,'b-');
legend('data','三次紧压','nature','外推','二次常量','二次指定')
xlabel('xdata');
ylabel('ydata');

%% 6.1
clear,clc;
syms x;
f=@(x)tan(cos((sqrt(5)+sin(x))./(1+x.^2)));
g=diff(f,x);
tol=1e-13;
max1=100;
x=(1+sqrt(5))/3;
h=0.1;
H(1)=h;
D(1)=(f(x+h)-f(x-h))/(2*h);
E(1)=0;
R(1)=0;
for n=1:2
    h=h/10;
    H(n+1)=h;
    D(n+1)=(f(x+h)-f(x-h))/(2*h);
    E(n+1)=abs(D(n+1)-D(n));
    R(n+1)=2*E(n+1)-(abs(D(n+1))+abs(D(n)+eps));
end
n=2;
while((E(n)>E(n+1))&&abs(R(n))>tol)&&n<max1
    h=h/10;
    n=n+1;
    H(n+1)=h;
    D(n+1)=(f(x+h)-f(x-h))/(2*h);
    E(n+1)=abs(D(n+1)-D(n));
    R(n+1)=2*E(n+1)-(abs(D(n+1))+abs(D(n)+eps));
end
figure;
plot([1:n+1],D);
xlabel('迭代次数');
ylabel('数值微分');
figure;
plot([1:n+1],E);
xlabel('迭代次数');
ylabel('误差');
figure
x=[0:0.001:2];
y=f(x);
plot(x,y);
hold on;
x=(1+sqrt(5))/3;
y=f(x);
plot(x,y,'o')
xlabel('xdata');
ylabel('ydata');
f1=matlabFunction(g);
figure
x=[0:0.001:2];
y=f1(x);
plot(x,y);
hold on;
x=(1+sqrt(5))/3;
y=f1(x);
plot(x,y,'o')
xlabel('xdata');
ylabel('ydata');

%% 7.1
clear,clc;
syms x;
g0=@(x)sin(x);g1=diff(g0,x);
a=0;b=pi/4;f=matlabFunction(sqrt(g1.^2+1));

M=10;s=0;h=(b-a)/M;
for k=1:M
    x1=a+h*(k-1);x2=a+h*k;s=s+h*(f(x1)+f(x2))/2;
end
format long
disp(s);

M=5;h=(b-a)/(2*M);s1=0;s2=0;
for k=1:M
    x=a+h*(2*k-1);s1=s1+f(x);
end
for k=1:M-1
    x=a+h*(2*k);s2=s2+f(x);
end
format long
s=h*(f(a)+f(b)+4*s1+2*s2)/3;
disp(s);

x=[0:0.01:1];
y=f(x);
plot(x,y);
xlabel('xdata');
ylabel('ydata');

%% 7.2
clear,clc
n=13;
for k=0:n
    x(k+1)=sqrt(k^2+1);
    y(k+1)=k^(1/3);
end
plot(x,y,'-o','LineWidth',5);xlabel('xdata');ylabel('ydata');
s=0;
e=0;
for k=1:n
    s=s+(y(k)+y(k+1))*(x(k+1)-x(k))/2;
    e=e+y(k+1)*(x(k+1)-x(k)).^3/12;
end
disp(s);
disp(e);

y=[0:0.01:2.5];
x=sqrt(y.^6+1);
hold on
plot(x,y,'LineWidth',5)
legend('f(x)','h(x)');

%% 9.1.1
clear,clc;
L=25000;
k=0.00003;
f=@(t,y)k*(L-y)*y;
a=0;
b=60;
ya=250;
M=300;
h=(b-a)/M;
T=zeros(1,M+1);
Y=zeros(1,M+1);
T=a:h:b;
Y(1)=ya;
for j=1:M
    Y(j+1)=Y(j)+h*feval(f,T(j),Y(j));
end
E=[T' Y'];
plot(T,Y,'LineWidth',2);
xlabel('时间t','FontSize',10);
ylabel('感染者数目y','FontSize',10);
title('感染者数目随时间变化','FontSize',20);
ave=vpa(sum(Y)/(M+1));
% y = L/(1+exp(wx+b)),log(1/(y/L)-1)=wx+b
Z=log(1./(Y/L)-1);
X=[ones(M+1,1) T'];
W=pinv(X)*Z';
t0=[0:0.2:60];
y0=L./(1+exp([ones(length(t0),1) t0']*W));
hold on
plot(t0,y0,'LineWidth',2);
legend('欧拉方法','logistic拟合');
syms t;
y=L./(1+exp(W(1)+t*W(2)));
ave1=vpa(int(y,a,b)/(b-a));

%% 9.1.2
clear,clc;
f=@(t,y,h)(y*1.3-0.25*y^2-0.0001*y*h);
a=0;
b=20;
ya=25;
M=100;
h=(b-a)/M;
T=zeros(1,M+1);
Y=zeros(1,M+1);
H=zeros(1,M+1);
T=a:h:b;
H(1)=0;
Y(1)=ya;
for j=1:M 
    Y(j+1)=Y(j)+h*f(T(j),Y(j),H(j));
    H(j+1)=H(j)+h/2*(Y(j)+Y(j+1));
end
plot(T,Y,'LineWidth',2);
hold on

ya=20;
T=zeros(1,M+1);
Y=zeros(1,M+1);
H=zeros(1,M+1);
T=a:h:b;
H(1)=0;
Y(1)=ya;
for j=1:M 
    Y(j+1)=Y(j)+h*f(T(j),Y(j),H(j));
    H(j+1)=H(j)+h/2*(Y(j)+Y(j+1));
end
plot(T,Y,'LineWidth',2);
hold on

ya=30;T=zeros(1,M+1);Y=zeros(1,M+1);H=zeros(1,M+1);
T=a:h:b;H(1)=0;Y(1)=ya;
for j=1:M 
    Y(j+1)=Y(j)+h*f(T(j),Y(j),H(j));
    H(j+1)=H(j)+h/2*(Y(j)+Y(j+1));
end
%plot(T,Y,'LineWidth',2);
legend('y0=25','y0=20','y0=30');

%% 9.2.1
clear,clc;
f=@(t,y)3*(y+t);
a=0;
ya=1;
M=20;
h=0.1;
T=zeros(1,M+1);
Y=zeros(1,M+1);
R=zeros(1,M+1);
T=a:h:a+M*h;
Y(1)=ya;
R(1)=0;
for j=1:M 
    k1=f(T(j),Y(j));
    k2=f(T(j+1),Y(j)+h*k1);
    Y(j+1)=Y(j)+h/2*(k1+k2);
    R(j+1)=Y(j+1)-(4/3*exp(3*T(j+1))-T(j+1)-1/3);
end
figure;
plot(T,Y,'LineWidth',2);
hold on
mean(R)

f=@(t,y)3*(y+t);
a=0;
ya=1;
M=40;
h=0.05;
T=zeros(1,M+1);
Y=zeros(1,M+1);
R=zeros(1,M+1);
T=a:h:a+M*h;
Y(1)=ya;
R(1)=0;
for j=1:M 
    k1=f(T(j),Y(j));
    k2=f(T(j+1),Y(j)+h*k1);
    Y(j+1)=Y(j)+h/2*(k1+k2);
    R(j+1)=Y(j+1)-(4/3*exp(3*T(j+1))-T(j+1)-1/3);
end
plot(T,Y,'LineWidth',2);
hold on
mean(R)

t=[0:0.01:2];
y=4/3*exp(3*t)-t-1/3;
plot(t,y,'LineWidth',2);
legend('h1','h2','h0');

%% 9.3.1
clear,clc;
f=@(t,y)3*(y+t);
a=0;
ya=1;
M=20;
h=0.1;
T=zeros(1,M+1);
Y=zeros(1,M+1);
R=zeros(1,M+1);
T=a:h:a+M*h;
Y(1)=ya;
R(1)=0;
for j=1:M 
    k1=h*f(T(j),Y(j));
    k2=h*f(T(j)+h/2,Y(j)+k1/2);
    k3=h*f(T(j)+h/2,Y(j)+k2/2);
    k4=h*f(T(j)+h,Y(j)+k3);
    Y(j+1)=Y(j)+(k1+2*k2+2*k3+k4)/6;
    R(j+1)=Y(j+1)-(4/3*exp(3*T(j+1))-T(j+1)-1/3);
end
figure;
plot(T,Y,'LineWidth',2);
hold on
mean(R)

f=@(t,y)3*(y+t);
a=0;
ya=1;
M=40;
h=0.05;
T=zeros(1,M+1);
Y=zeros(1,M+1);
R=zeros(1,M+1);
T=a:h:a+M*h;
Y(1)=ya;
R(1)=0;
for j=1:M 
    k1=f(T(j),Y(j));
    k2=f(T(j)+h/2,Y(j)+h*k1/2);
    k3=f(T(j)+h/2,Y(j)+h*k2/2);
    k4=f(T(j)+h,Y(j)+h*k3);
    Y(j+1)=Y(j)+h/6*(k1+2*k2+2*k3+k4);
    R(j+1)=Y(j+1)-(4/3*exp(3*T(j+1))-T(j+1)-1/3);
end
plot(T,Y,'LineWidth',2);
hold on
mean(R)

t=[0:0.01:2];
y=4/3*exp(3*t)-t-1/3;
plot(t,y,'LineWidth',2);
legend('h1','h2','h0');

%% 11.1
clear,clc;
n=100;
e=1e-10;
A=[4,-1,1;-1,3,-2;1,-2,3];
x(:,1)=[1,1,1]';
lamda=0;
z=zeros(3,1);
for i=1:n
    y(:,i)=A*x(:,i);
    x(:,i+1)=y(:,i)/max(y(:,i));
    if sqrt((x(:,i+1)-x(:,i))'*(x(:,i+1)-x(:,i)))<e
        z=x(:,i);
        lamda=max(y(:,i));
        disp(i);
        break;
    end
end
plot(max(y),'LineWidth',2)
xlabel('迭代次数');
ylabel('特征值');

disp(z);
disp(lamda);
%% 11.2
clear,clc;
n=100;
e=1e-10;
A=[4,-1,1;-1,3,-2;1,-2,3];
x(:,1)=[1,1,1]';
lamda(1)=0;
lm=0;
z=zeros(3,1);
for i=2:n
    x(:,i)=A*x(:,i-1);
    x(:,i)=x(:,i)/sqrt(x(:,i)'*x(:,i));
    lamda(i)=(x(:,i)'*A*x(:,i))/(x(:,i)'*x(:,i));
    if (lamda(i)-lamda(i-1))<e
        z=x(:,i);
        lm=lamda(i);
        break;
    end
end
plot(lamda,'LineWidth',2)
xlabel('迭代次数');
ylabel('特征值');
disp(z);
disp(lm);

%% 11.3
clear,clc;
n=100;
e=1e-10;
A=[4,-1,1;-1,3,-2;1,-2,3];
q=5.9;
B=inv(A-q*eye(3));
x(:,1)=[1,1,1]';
lamda=0;
z=zeros(3,1);
for i=1:n
    y(:,i)=B*x(:,i);
    x(:,i+1)=y(:,i)/max(y(:,i));
    if sqrt((x(:,i+1)-x(:,i))'*(x(:,i+1)-x(:,i)))<e
        z=x(:,i);
        lamda=max(y(:,i));
        break;
    end
end
plot(q+1./max(y),'LineWidth',2)
xlabel('迭代次数');
ylabel('特征值');
disp(z);
disp(q+1/lamda);

