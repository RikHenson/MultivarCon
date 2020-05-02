function [MVconn,MVconn_null] = computeMVconn(X,Y,opt)
% this is mainly a wrapper function: given the input MV patterns for 2
% regions, it computes MVPD, LPRD, dCor, RC and MIM and also FC and FC_SVD
% the function also simulates the null hypothesis of no functional
% connectivity by randomly shuffling the time points (first component) in
% each run of every subject. This assumes no temporal autocorrelation in
% data and would not be appropriate for temporally smooth data.
% results from permuted data are returned in MVconn_null, which contains
% one null number per subject per iteration. The number of iterations are
% specified in opt.nRandomisation.
% 
% inputs:
%        X: a cell array, one entry per subject. X{s} is the data from
%        subject "s" and contains mutiple cell arrays, one per run. For
%        example X{2}{1} is the data from run1 of subject2. The data is 
%        time x voxels.
%        Y : the same as X for the second region
% 
% outputs:
%        MVconn: contains the following fields:
%               FC, FCSVD, FCCCA, MVPD, LPRD, dCor, RCA (MIM, ImCohSVD, MVLagCoh, LagCohSVD
%               each is a column vector with one number per subject.
%        MVconn_null: contains the same fields as MVconn. Each would be a
%        matrix with size of nSubjects x nRandomisations.
% Hamed Nili
if ~isfield(opt,'nRandomisation')
    opt.nRandomisation = 1;
end
if ~isfield(opt,'zscore')
    opt.zscore = 0;
end

nSub = length(X);

% zscore patterns within each run if opt.zscore is set in the options
if opt.zscore
    for r = 1:numel(X)
        X{r} = zscore(X{r},0,2);
        Y{r} = zscore(Y{r},0,2);
    end
end

% Calculate connectivity on given data
for s=1:nSub
    if ~isfield(opt,'segleng') 
        [mvpd(s,1),lprd(s,1),fc(s,1),fc_svd(s,1),fc_cca(s,1)] = data2mvpd_lprd_fc(X{s},Y{s},opt); 
        [dcor(s,1),dcor_u(s,1)] = data2dCor(X{s},Y{s});
        [rc(s,1),~] = data2rc(X{s},Y{s},'Correlation');
    else
        [mim(s,1),imcoh_svd(s,1),mvlagcoh(s,1),lagcoh_svd(s,1)] = data2lagconn(X{s},Y{s},opt);
    end
end

% Calculate connectivity when X and Y independent random noise (since
% some connectivity measures, eg dCor, not bounded by 0 or -1)
bmvpd       = NaN(nSub,opt.nRandomisation);
blprd       = NaN(nSub,opt.nRandomisation);
bfc         = NaN(nSub,opt.nRandomisation);
bfc_svd     = NaN(nSub,opt.nRandomisation);
bfc_cca     = NaN(nSub,opt.nRandomisation);
bdcor       = NaN(nSub,opt.nRandomisation);
brc         = NaN(nSub,opt.nRandomisation);
bmim        = NaN(nSub,opt.nRandomisation);
bimcoh_svd  = NaN(nSub,opt.nRandomisation);
bmvlagcoh   = NaN(nSub,opt.nRandomisation);
blagcoh_svd = NaN(nSub,opt.nRandomisation);

if opt.nRandomisation == 1 %if only one, then don't bother with parfor like below
    if nSub < 20
        warning('May not be sufficient subjects/randomisations to estimate null properly')
    end
    iter = 1;
    for s=1:nSub % Ensure reasonably accurate estimate
        bX = {}; bY = {};
        for r=1:length(X{s})
            bX{r} = X{s}{r}(randperm(size(X{s}{r},1)),:);
            bY{r} = Y{s}{r}(randperm(size(Y{s}{r},1)),:);
        end
        if ~isfield(opt,'segleng')
            [bmvpd(s,iter),blprd(s,iter),bfc(s,iter),bfc_svd(s,iter),bfc_cca(s,iter)] = data2mvpd_lprd_fc(bX,bY,opt);
            [bdcor(s,iter),bdcor_u(s,iter)] = data2dCor(bX,bY);
            [brc(s,iter),~] = data2rc(bX,bY,'Correlation');
        else
            [bmim(s,iter),bimcoh_svd(s,iter),bmvlagcoh(s,iter),blagcoh_svd(s,iter)] = data2lagconn(bX,bY,opt);
        end
    end
elseif opt.nRandomisation > 1
    for s=1:nSub % Ensure reasonably accurate estimate
        fprintf('null subject %d from %d \n',s,nSub)
        parfor iter = 1:opt.nRandomisation
            bX = {}; bY = {};
            for r=1:length(X{s})
                bX{r} = X{s}{r}(randperm(size(X{s}{r},1)),:);
                bY{r} = Y{s}{r}(randperm(size(Y{s}{r},1)),:);
            end
            if ~isfield(opt,'segleng')
                [tmp_bmvpd{iter},tmp_blprd{iter},tmp_bfc{iter},tmp_bfc_svd{iter},tmp_bfc_cca{iter}] = data2mvpd_lprd_fc(bX,bY,opt);
                [tmp_bdcor{iter},~] = data2dCor(bX,bY);
                [tmp_brc{iter},~] = data2rc(bX,bY,'Correlation');
            else
                [tmp_bmim{iter},tmp_bimcoh_svd{iter},tmp_bmvlagcoh{iter},tmp_blagcoh_svd{iter}] = data2lagconn(bX,bY,opt);
            end
        end
        if ~isfield(opt,'segleng')
            bmvpd(s,:) = cat(2,tmp_bmvpd{:});
            blprd(s,:) = cat(2,tmp_blprd{:});
            bfc(s,:) = cat(2,tmp_bfc{:});
            bfc_svd(s,:) = cat(2,tmp_bfc_svd{:});
            bfc_cca(s,:) = cat(2,tmp_bfc_cca{:});
            bdcor(s,:) = cat(2,tmp_bdcor{:});
            brc(s,:) = cat(2,tmp_brc{:});
        else
            bmim(s,:) = cat(2,tmp_bmim{:});
            bimcoh_svd(s,:) = cat(2,tmp_bimcoh_svd{:});
            bmvlagcoh(s,:) = cat(2,tmp_bmvlagcoh{:});
            blagcoh_svd(s,:) = cat(2,tmp_blagcoh_svd{:});
        end
    end
end
fprintf('\n')

if ~isfield(opt,'segleng')
    MVconn.FC = fc;
    MVconn.FCSVD = fc_svd;
    MVconn.FCCCA = fc_cca;
    MVconn.MVPD = mvpd;
    MVconn.LPRD = lprd;
    MVconn.dCor = dcor;
    MVconn.RCA = rc;
else
    MVconn.MIM = mim;
    MVconn.ImCohSVD = imcoh_svd;
    MVconn.MVLagCoh = mvlagcoh;
    MVconn.LagCohSVD = lagcoh_svd;
end

if ~isfield(opt,'segleng')
    MVconn_null.FC = mean(bfc,2);
    MVconn_null.FCSVD = mean(bfc_svd,2);
    MVconn_null.FCCCA = mean(bfc_cca,2);
    MVconn_null.MVPD = mean(bmvpd,2);
    MVconn_null.LPRD = mean(blprd,2);
    MVconn_null.dCor = mean(bdcor,2);
    MVconn_null.RCA = mean(brc,2);
else
    MVconn_null.MIM = mean(bmim,2);
    MVconn_null.ImCohSVD = mean(bimcoh_svd,2);
    MVconn_null.MVLagCoh = mean(bmvlagcoh,2);
    MVconn_null.LagCohSVD = mean(blagcoh_svd,2);
end

return

