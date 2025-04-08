clear all;
tic

N=4;  kappa_a=1;E1 = 0.01;k=1;g=0.01;E2=-0.01;detla1=0;
M=50;

% Define cavity field and atomic operators
a = tensor(destroy(N),identity(2));
b = tensor(identity(N),destroy(2));
 for m=1:M
     m
     ph1=20*m/M-10;
     
     
     xx(m)=ph1;
     for m1=1:M
          
         ph=10*m1/M-5;
     xx1(m1)=ph;
      
   
     
    
        
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
        
           
      
     end
  end
 
% g(m)=expect(a'*a'*a*a,rhoss);
%  e(m)=expect(a'*a+b'*b,rhoss);

for n=1:M
    phh=10*n/M-5;
    phx(n)=phh;
    ph11=2*phh-4*pi;
    phy(n)= ph11;
    phy1(n)= ph11;
end



 figure
mesh(xx1,xx,log10(abs(gg)))
hold on
plot(phx,phy)
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
