function [MVconn,MVconn_null] = computeMVconn(Ya,Yb,opt)

if ~isfield(opt,'nRandomisation')
    opt.nRandomisation = 1;
end


% Calculate connectivity on data given
for s=1:length(Ya)
    [mvpd(s,1),uvpd(s,1),fc(s,1)] = data2mvpd(Ya{s},Yb{s},opt); 
    [dcor(s,1),dcor_u(s,1)] = data2dCor(Ya{s},Yb{s});
    [rc(s,1),~] = data2rc(Ya{s},Yb{s},'correlation');
end

% Calculate connectivity when Ya and Yb independent random noise (since
% some connectivity measures, eg Dcor, not bounded by 0 or -1)
if length(Ya)<20
    warning('Insufficient subjects (<20) to estimate baseline error')
else
    for iter = 1:opt.nRandomisation
        for s=1:length(Ya) % Ensure reasonably accurate estimate
            bYa = {}; bYb = {};
            for r=1:length(Ya{s})
                bYa{r} = Ya{s}{r}(randperm(size(Ya{s}{r},1)),:);
                bYb{r} = Yb{s}{r}(randperm(size(Yb{s}{r},1)),:);
            end
            [bmvpd(s,iter),buvpd(s,iter),bfc(s,iter)] = data2mvpd(bYa,bYb,opt);
            [bdcor(s,iter),bdcor_u(s,iter)] = data2dCor(bYa,bYb);
            [brc(s,iter),~] = data2rc(bYa,bYb,'correlation');
        end
    end
end


MVconn.FC = fc;MVconn.MVPD = mvpd;
MVconn.UVPD = uvpd;
MVconn.dCor = dcor;MVconn.dCor_univar = dcor_u;
MVconn.RCA = rc;
MVconn_null.FC = bfc;MVconn_null.MVPD = bmvpd;
MVconn_null.UVPD = buvpd;
MVconn_null.dCor = bdcor;MVconn_null.dCor_univar = bdcor_u;
MVconn_null.RCA = brc;









