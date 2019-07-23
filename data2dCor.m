function [dCor,dCor_u,fc] = data2dCor(Ya,Yb)
% calculates the multivariate distance correlation between two multivariate
% time series (Geerlings et al. 2016, NI)
%
% Inputs are the multivariate time courses for two brain ROIs. Both Ya and
% Yb are cell arrays containing the MV timeseries for different runs. In
% fMRI, these time series are assumed to be pre-processed raw data after
% high-pass filtering and regressing out nuisance variables (e.g. motion
% params, avg CSF, avg WM time series, etc.). the script first zscores the
% data across voxels/electrodes in the two ROIs and for each run
% independently. Then it computes the Euclidean distances
% between all pairs of timepoints. These matrices will be U-centered. This
% operation ensures that the correlation between matrices are not biased by
% the number of voxels in ROIs. the output is the correlation between the
% two matrices. Hamed Nili, based on the code available at
% http://imaging.mrc-cbu.cam.ac.uk/imaging/Geerligs_DistCor
% Furthermore, it computes the correlation between the centered distance
% matrices obtained from the UV time courses (voxel/electrode averaged time
% series). D-centering is used for computing the UV distance correlation.
% Hamed Nili, based on the codes available at 
% http://imaging.mrc-cbu.cam.ac.uk/imaging/Geerligs_DistCor
% the function returns the run-averaged values.


%% compute the dcor and univariate distance matrices for each run
nruns = numel(Ya); % number of runs

% zscore the data (across voxels in each ROI and for each run)
for r=1:nruns
    %Ya_zs = zscore(Ya{r},0,2);
    Ya_zs = Ya{r};Yb_zs = Yb{r};
    %Yb_zs = zscore(Yb{r},0,2);
    % get the voxel-average ts
    ts_a = mean(Ya{r},2);ts_b = mean(Yb{r},2);
       
    % compute the pearson correlation
    [r_uv(r),~] = corr(ts_a,ts_b);
    
    % MV distance correlation (U-centering)
    dcor(r) = dcor_uc(Ya_zs,Yb_zs);
    
    % UV distance correlation (D-centering)
    dcor_u(r) = dcor_dc(ts_a,ts_b);
end

dCor_u=mean(dcor_u);% run-averaged values
dCor=mean(dcor); % run-averaged values
fc = mean(r_uv);

    function [R] = dcor_dc(X,Y)
        % [R] = dcor_dc(X,Y)
        % Computes the double centered distance correlation between X and Y.
        % Rows represent the examples, and columns the variables.
        
        % Based on: http://www.mathworks.com/matlabcentral/fileexchange/49968-dcorr--x--y--
        % and, the R package energy and papers by Szekely and Rizzo (2007; 2013 and 2014).
        
        % Author: Linda Geerligs (lindageerligs@gmail.com), Date: 22-10-2015
        
        a = pdist2(X, X);
        b = pdist2(Y, Y);
        n=size(X,1);
        
        A = Dcenter(a);
        B = Dcenter(b);
        
        dcovXY = sum(sum(A.*B)) ./ (n.^2);
        dvarX = sum(sum(A.*A)) ./ (n.^2);
        dvarY = sum(sum(B.*B)) ./ (n.^2);
        
        R=sqrt(dcovXY / sqrt(dvarX * dvarY));
        
        function A=Dcenter(a)
            A = a - bsxfun(@plus,mean(a),mean(a,2))+mean(a(:));
        end
        
    end

    function [R,T,df] = dcor_uc(X,Y)
        % [R,T,p,df] = dcor_uc(X,Y)
        % Computes the U-centered (bias corrected) distance correlation between X and Y. Also outputs the
        % associated t-value (T). Rows represent the examples, and columns the variables.
        
        % Based on: http://www.mathworks.com/matlabcentral/fileexchange/49968-dcorr--x--y--
        % and, the R package energy and papers by Szekely and Rizzo (2007; 2013 and 2014).
        
        % Author: Linda Geerligs (lindageerligs@gmail.com), Date: 22-10-2015
        
        a = pdist2(X, X);
        b = pdist2(Y, Y);
        n=size(X,1);
        
        A = Ucenter(a,n);
        B = Ucenter(b,n);
        
        dcovXY = sum(sum(A.*B)) ./ (n*(n-3));
        dvarX = sum(sum(A.*A)) ./ (n*(n-3));
        dvarY = sum(sum(B.*B)) ./ (n*(n-3));
        
        R=dcovXY / sqrt(dvarX * dvarY);
        
        df=(n*(n-3))/2 -1;
        T=sqrt(df).*(R./sqrt(1-R.^2));
        
        if R<0
            R=0;
        else
            R=sqrt(R);
        end
        
        function A = Ucenter(a,n)
            m = sum(a,2);
            M = sum(m)/((n - 1) * (n - 2));
            m = m./(n-2);
            A = a - repmat(m,[1 n]);
            A = A - repmat(m,[1 n])';
            A = A+M;
            A(eye(size(A))==1)=0;
        end
        
    end

end

