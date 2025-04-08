clear all;
tic


N=6;  kappa_a=1;g=6;E1=0.1;k=1;ph=0;detla1=0;ph1=0;
M=80;
% Define cavity field and atomic operators
a = tensor(destroy(N),identity(2));
b = tensor(identity(N),destroy(2));
 for m=1:M
     m
    E2=0.02*m/M-0.01;
   
     xx(m)=E2;
     
      

    
        
H =  detla1*a'*a+ 2*detla1*b'*b ...
+g*(a'*a'*b+b'*a*a)+  E1*(a'*exp(i*ph)+a*exp(-i*ph))+  E2*(b'*exp(i*ph1)+b*exp(-i*ph1));

  LH = -i * (spre(H) - spost(H));
          L1=kappa_a/2*(2*spre(a)*spost(a')-spre(a'*a)-spost(a'*a));
       
          L2=kappa_a/2*(2*spre(b)*spost(b')-spre(b'*b)-spost(b'*b));
       
       
L = LH+L1+L2;
% Find steady state
rhoss = steady(L);
   
         gg(m)=trace((a'*a'*a*a)*rhoss)/(trace(a'*a*rhoss))^2;
%          
         nn(m)=rhoss(3,3);
           
  
   
  end
 



  figure
  plot(xx,log10(abs(gg)),xx,log10(abs(nn)))

%   figure
%  plot(xx,nn)


