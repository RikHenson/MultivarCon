function [Ttilde,gof]=ridgeregmethod(Ya,Yb,lambdas); 
% Computation of the estimated linear transformation T (i.e. such that Yb=T*Ya) 
% via a cross-validated version of the ridge regression method.
% 
% INPUT
% Ya:         time series for the ROI1
% Yb:         time series for the ROI2
% lambdas:    set of possible regulariation parameter
% OUTPUT
% Ttilde:     estimated transformation
% gof:        goodness-of-fit
% Alessio Basti 15/07/2019

Ya=Ya';
Yb=Yb';
k=1;
for i=lambdas                    
   H=Ya'*pinv(Ya*Ya'+i*eye(size(Ya,1)))*Ya;                                      
   for j=1:size(Ya,2)
       A(j,j)=1/(1-H(j,j));
   end            
   CrossValid(k)=(norm(A*((eye(size(Ya,2))-H)*Yb'),'fro'))^2;                    
   k=k+1;
end
[B C]=min(CrossValid);
optlambda=lambdas(C);
gof=(1-B/(size(Ya,2)*size(Yb,1)));
Ttilde=Yb*Ya'*pinv(Ya*Ya'+optlambda*eye(size(Ya,1)));

return
