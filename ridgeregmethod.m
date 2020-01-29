function [Ttilde,gof]=ridgeregmethod(X,Y,lambdas); 
% Computation of the estimated linear transformation T (i.e. such that Y=T*X) 
% via a cross-validated version of the ridge regression method.
% 
% INPUT
% X:         time series for the ROI1 (time x voxel matrix)
% Y:         time series for the ROI2
% lambdas:   set of possible regulariation parameters
% OUTPUT
% Ttilde:    estimated transformation
% gof:       goodness-of-fit
% Alessio Basti 15/07/2019

X=X';
Y=Y';
k=1;
for i=lambdas                    
   H=X'*pinv(X*X'+i*eye(size(X,1)))*X;                                      
   for j=1:size(X,2)
       A(j,j)=1/(1-H(j,j));
   end            
   CrossValid(k)=(norm(A*((eye(size(X,2))-H)*Y'),'fro'))^2;                    
   k=k+1;
end
[B C]=min(CrossValid);
optlambda=lambdas(C);
gof=(1-B/(size(X,2)*size(Y,1)));
Ttilde=Y*X'*pinv(X*X'+optlambda*eye(size(X,1)));

return
