datamatrix = zeros(6,1270)
datamatrix(2,1:length(u)) = u;
datamatrix(3,1:length(u)) = u;
b = u( 1:find( u > 100000, 1) );

function Rdot = f(R,t)
 v = 8.56e-7
y = 69.5e-3
pl = 0.9
fs = 10e-15 
Rdot(1) = R(2)
Rdot(2) = ppval(pp,wat(:,1))/pl*1/R(1) + 

% Runge-Kutta Method
function[tt, yy] =runge_kutta(f,t0,y0,h,N)
k=N+1;
tt=zeros(k,1);
yy=zeros(k,1);
tt(1)=t0;yy(1)=y0;
for i=2:k
  tt(i)=tt(i-1)+h;
  m1=f(tt(i-1),yy(i-1));
  m2=f(tt(i-1)+(h/2),yy(i-1)+(h/2)*m1);
  m3=f(tt(i-1)+(h/2),yy(i-1)+(h/2)*m2);
  m4=f(tt(i),yy(i-1)+h*m3);
  yy(i)=yy(i-1)+h*(m1+2*m2+2*m3+m4)/6;
end
end