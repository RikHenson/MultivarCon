function [Ttilde,Ytilde]=ridgeregmethod(X,Y,lambdas,binYtilde); 
% Computation of the estimated linear transformation T (i.e. such that Y=T*X) 
% via a cross-validated version of the ridge regression method.
% 
% INPUT
% X:         time series for the ROI1 (time x voxel matrix)
% Y:         time series for the ROI2
% lambdas:   set of possible regulariation parameters
% binYtilde: if it is equal to 1, it leads to have an estimate of Y
% OUTPUT
% Ttilde:    estimated transformation
% Ytilde:    estimated time series 
% Alessio Basti 25/06/2020

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
%gof=(1-B/(size(X,2)*size(Y,1)));
optlambda=lambdas(C);
Ttilde=Y*X'*pinv(X*X'+optlambda*eye(size(X,1)));
Ytilde=zeros(size(Y));
if(binYtilde==1)
    for j=1:size(X,2)
        Xapp=X;Xapp(:,j)=[]; 
        Yapp=Y;Yapp(:,j)=[];
        % in case of autocorrelated data, a neighbourhood of j could be removed, and thus a CV-approach
        % different from LOOCV could be used
        Tjtilde=Yapp*Xapp'*pinv(Xapp*Xapp'+optlambda*eye(size(Xapp,1)));
        Ytilde(:,j)=Tjtilde*X(:,j);
    end
end
Ytilde=Ytilde';

% In addition to remove a neighbourhood of a i-th stimulus/time point, one
% could also remove a neighbourhood of j, for every j; in this case the RDM
% between i and j should be computed in the loop. Below an example with
% radius equal to 0
% RDMYtilde=zeros(1,size(X,2)*(size(X,2)-1)/2);
% k=1;
% for i=1:size(X,2)-1
%     for j=i+1:size(X,2)
%         Xapp=X;Xapp(:,[i,j])=[];
%         Yapp=Y;Yapp(:,[i,j])=[];
%         Tijtilde=Yapp*Xapp'/(Xapp*Xapp'+optlambda*eye(size(Xapp,1)));
%         RDMYtilde(k)=1-corr(Tijtilde*X(:,i),Tijtilde*X(:,j));
%         k=k+1;
%     end
% end

return
