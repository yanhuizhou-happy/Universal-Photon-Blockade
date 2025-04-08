clear all;
tic

%N=6;  kappa_a=1;g=0.1;E1=0.01;k=1;ph=0;ph1=0; E2=-E1^2/g; C00=1; %fig2(a)
%N=6;  kappa_a=1;g=1;E1=0.05;k=1;ph=0;ph1=0; E2=-E1^2/g; C00=1;   %fig2(b)
N=6;  kappa_a=1;g=6;E1=0.1;k=1;ph=0;ph1=0; E2=-E1^2/g; C00=1;

M=80;

% Define cavity field and atomic operators
a = tensor(destroy(N),identity(2));
b = tensor(identity(N),destroy(2));
 for m=1:M
     m
     detla1=12*m/M-6;
    
     xx(m)=detla1;
     
      
     
    
        
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
        p1(m)=rhoss(3,3);
           
      nn(m)=trace(a'*a*rhoss);
   
 end
 
for m=1:M
    detla1=12*m/M-6;
  xx1(m)=detla1;
ggg(m)=(((4*(E1^2 + g^2)* k + k^3)^2 + 8* (8 *(2* E1^2 + g^2)^2 + 4 *(3*E1^2 - g^2)*k^2 + ...
        3*k^4)*detla1^2 - 16*(32*E1^2 + 16*g^2 - 9*k^2)*detla1^4 + ...
     256*detla1^6)*(E2^2*g^2*(k^2 + 4*detla1^2) + E1^4 *(k^2 + 16 *detla1^2) + ...
     2*E1^2*E2*g *((k^2 + 8 *detla1^2)* cos(2*ph - ph1) - 2 *k* detla1* sin(2 *ph - ph1))))/(C00^2 *E1^4 *(16 *E2^2 *g^2 + (4 *g^2 + k^2)^2 + ...
     4 *(-16* g^2 + 5 *k^2)* detla1^2 + 64 *detla1^4 - ...
     8* E2 *g *(4* g^2 + k^2 - 8 *detla1^2)* cos(2* ph - ph1) + ...
     48*E2* g* k* detla1* sin(2*ph - ph1))^2);
 
 p11(m)=(4*C00^2*E1^2 *(16* E2^2 *g^2 + (4* g^2 + k^2)^2 + 4 *(-16 *g^2 + 5* k^2)* detla1^2 + 64 *detla1^4 - ...
   8* E2* g* (4 *g^2 + k^2 - 8 *detla1^2)* cos(2* ph - ph1) + 48* E2* g* k *detla1* sin(2 *ph - ph1)))/((4 *(E1^2 + g^2)*k + ...
   k^3)^2 + 8*(8* (2*E1^2 + g^2)^2 + 4*(3* E1^2 - g^2)*k^2 + 3* k^4)* detla1^2 - 16 *(32 *E1^2 + 16* g^2 - 9* k^2) *detla1^4 + 256* detla1^6);
end

%   figure
%   plot(xx,log10(abs(gg)),xx,log10(abs(p1)))

 hold on 
 plot(xx1,log10(abs(nn)),xx1,log10(abs(p11)))

