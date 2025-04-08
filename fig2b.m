clear all;
tic

N=4;  kappa_a=1;E1 = 0.05;k=1;detla1=0;ph=0;ph1=0;
M=50;

% Define cavity field and atomic operators
a = tensor(destroy(N),identity(2));
b = tensor(identity(N),destroy(2));
 for m=1:M
     m
     g=6*m/M-3;
     
     xx(m)=g;
     for m1=1:M
          E2=0.04*m1/M-0.02;
         
     xx1(m1)=E2;
      
   
     
    
        
H =  detla1*a'*a+ 2*detla1*b'*b ...
+g*(a'*a'*b+b'*a*a)+  E1*(a'*exp(i*ph)+a*exp(-i*ph))+  E2*(b'*exp(i*ph1)+b*exp(-i*ph1));

  LH = -i * (spre(H) - spost(H));
          L1=kappa_a/2*(2*spre(a)*spost(a')-spre(a'*a)-spost(a'*a));
       
          L2=kappa_a/2*(2*spre(b)*spost(b')-spre(b'*b)-spost(b'*b));
       
       
L = LH+L1+L2;
% Find steady state
rhoss = steady(L);
   
         gg(m,m1)=trace((a'*a'*a*a)*rhoss)/(trace(a'*a*rhoss))^2;
%          
        nn(m,m1)=trace(a'*a*rhoss);
           
      
     end
  end
 
% g(m)=expect(a'*a'*a*a,rhoss);
%  e(m)=expect(a'*a+b'*b,rhoss);

for n=1:M
    Ex=0.04*n/M-0.02;
    fb(n)=Ex;
    gn=-E1^2/(Ex);
    gy(n)=gn;
end

 figure
mesh(xx1,xx,log(abs(gg)))
hold on
plot(fb,gy)
%  figure
% mesh(xx1,xx,log(abs(nn)))

% figure
% plot(xx1,p0,xx1,(p01),xx1,(p02),xx1,(p03))
%  figure
%  plot(xx1,log(abs(bb)),xx1,log(abs(nn1)))
%  figure
%  plot(xx1,p10,xx1,nn)
%   figure
%  plot(xx1,p01,xx1,nn1)
