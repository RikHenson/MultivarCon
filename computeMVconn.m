function [MVconn,MVconn_null] = computeMVconn(X,Y,opt)

if ~isfield(opt,'nRandomisation')
    opt.nRandomisation = 1;
end


% Calculate connectivity on data given
for s=1:length(X)
    if ~isfield(opt,'segleng') 
        [mvpd(s,1),gof(s,1),fc(s,1),fc_pc(s,1)] = data2mvpd_gof_fc(X{s},Y{s},opt); 
        [dcor(s,1),dcor_u(s,1)] = data2dCor(X{s},Y{s});
        [rc(s,1),~] = data2rc(X{s},Y{s},'Correlation');
    else
        [mim(s,1),imcoh(s,1),imcoh_pc(s,1)] = data2mim(X{s},Y{s},opt);
    end
end

% Calculate connectivity when X and Y independent random noise (since
% some connectivity measures, eg Dcor, not bounded by 0 or -1)
if length(X)<20
    warning('Insufficient subjects (<20) to estimate baseline error')
else
    for iter = 1:opt.nRandomisation
        if opt.nRandomisation > 1
            fprintf('iteration %d from %d \n',iter,opt.nRandomisation)
        end
        for s=1:length(X) % Ensure reasonably accurate estimate
            bX = {}; bY = {};
            for r=1:length(X{s})
                bX{r} = X{s}{r}(randperm(size(X{s}{r},1)),:);
                bY{r} = Y{s}{r}(randperm(size(Y{s}{r},1)),:);
            end
            if ~isfield(opt,'segleng') 
                [bmvpd(s,iter),bgof(s,iter),bfc(s,iter),bfc_pc(s,iter)] = data2mvpd_gof_fc(bX,bY,opt);
                [bdcor(s,iter),bdcor_u(s,iter)] = data2dCor(bX,bY);
                [brc(s,iter),~] = data2rc(bX,bY,'Correlation');
            else
               [bmim(s,iter),bimcoh(s,iter),bimcoh_pc(s,iter)] = data2mim(bX,bY,opt);
            end
        end
    end
end

if ~isfield(opt,'segleng')
    MVconn.FC = fc;
    MVconn.FCPC = fc_pc;
    MVconn.MVPD = mvpd;
    MVconn.GOF = gof;
    MVconn.dCor = dcor;
    MVconn.RCA = rc;
else
    MVconn.MIM = mim;
    MVconn.ImCoh = imcoh;
    MVconn.ImCohPC = imcoh_pc;
end

if ~isfield(opt,'segleng')
    MVconn_null.FC = mean(bfc,2);
    MVconn_null.FCPC = mean(bfc_pc,2);
    MVconn_null.MVPD = mean(bmvpd,2);
    MVconn_null.GOF = mean(bgof,2);
    MVconn_null.dCor = mean(bdcor,2);
    MVconn_null.RCA = mean(brc,2);
else
    MVconn_null.MIM = mean(bmim,2);
    MVconn_null.ImCoh = mean(bimcoh,2);
    MVconn_null.ImCohPC = mean(bimcoh_pc,2);
end

return
