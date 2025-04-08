clear all;
tic

N=6;   kappa_a=1;
n = 1;E1 = 0.01; E2=-0.000576592; k = 1; ph = 0; ph1 = 0.3;g=0.208414;detla1=n*g;
%n = 10;E1 = 0.01; E2=-0.00576592; k = 1; ph = 0; ph1 = 0.3;g=0.0208414;detla1=n*g;
%n = 5;E1 = 0.01; E2=-0.00991266; k = 1; ph = 0; ph1 = 0.3;g=0.0102441;detla1=n*g;



M=80;

% Define cavity field and atomic operators
a = tensor(destroy(N),identity(2));
b = tensor(identity(N),destroy(2));
 for m=1:M
     m
   detla1=4*m/M-2;
     
     xx(m)=detla1;
     
              
               
             %  E2 = -((E1^2 *(-5*k + 3*k* cos(4*ph-2*ph1)+ sqrt(2)* sqrt(k^2*(-7 + 9* cos(4*ph-2*ph1)))*cos(2*ph-ph1))*csc(2*ph-ph1))/(2*sqrt(2)* k^2));
               
              
         
          %E2 =((E1^2 *(5*k -3*k* cos(4*ph-2*ph1)+ sqrt(2)* sqrt(k^2*(-7 + 9* cos(4*ph-2*ph1)))*cos(2*ph-ph1))*csc(2*ph-ph1))/(2*sqrt(2)* k^2)); 
          
         
      
   
     
    
        
H =  detla1*a'*a+ 2*detla1*b'*b ...
+g*(a'*a'*b+b'*a*a)+  E1*(a'*exp(i*ph)+a*exp(-i*ph))+  E2*(b'*exp(i*ph1)+b*exp(-i*ph1));

  LH = -i * (spre(H) - spost(H));
          L1=kappa_a/2*(2*spre(a)*spost(a')-spre(a'*a)-spost(a'*a));
       
          L2=kappa_a/2*(2*spre(b)*spost(b')-spre(b'*b)-spost(b'*b));
       
       
L = LH+L1+L2;
% Find steady state
rhoss = steady(L);
   
        
                gg(m)=trace((a'*a'*a*a)*rhoss)/(trace(a'*a*rhoss))^2;
               
         
        nn(m)=rhoss(3,3);
        
           
      
end
 
% g(m)=expect(a'*a'*a*a,rhoss);
%  e(m)=expect(a'*a+b'*b,rhoss);



hold on
plot(xx,log10(abs(gg)),xx,log10(abs(nn)))

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
