function [y,cy,fnum] = plot_2reg_2vox(S);


try cm = S.cm; catch cm = colormap('jet'); end
Nc = size(cm,1);

try y = S.y;
    N = size(y,1);
catch
    try C = S.C; 
        N = Nc;
        y = mvnrnd(zeros(N,4),C);
    catch error('Pass either data in S.y or data covar in S.C'); 
    end
end

try cc = S.cc; catch error('Pass contrast across voxels in S.cc'); end

try smoothness = S.smoothness; catch smoothness = 0; end
try categories = S.categories; catch categories = 0; end
try timepoints = S.timepoints; catch timepoints = 0; end
try cat_bound = S.cat_bound; catch cat_bound = 0; end
try SNR = S.SNR; catch SNR = Inf; end
try fnum = S.fnum; catch fnum = 0; end

if fnum
    fnum = figure(fnum); clf
else
    fnum = figure;
end

%CV = S.C;
CV = corrcoef(y);
figure(fnum),imagesc(CV); colormap('gray'); caxis([-1 1]), %axis off,
set(gca,'XTick',[]); set(gca,'YTick',[]);
line([0 5]',[2.5 2.5]','Color',[0 0 0]); line([2.5 2.5]',[0 5]','Color',[0 0 0])
colorbar

if smoothness
    for r = 1:size(y,2)
        y(:,r) = smooth(y(:,r),smoothness);
    end
end

my = max(abs(y(:)));

figure(fnum.Number+1), clf
for r=1:2
    subplot(1,2,r); hold on
    axis square
    axis([-my my -my my])
    set(gca,'XTick',0); set(gca,'YTick',0);
    rr = (r-1)*2 + 1;
    
    mker = '';
    for p=1:N
        mker(p) = 'o';
        if categories == 1
            if r==2
                if y(p,rr+1) - y(p,rr) > cat_bound*my, mker(p) = '*'; 
                elseif y(p,rr+1) - y(p,rr) < -cat_bound*my, mker(p) = '+'; end
            else
                if y(p,rr+1) + y(p,rr) > cat_bound*my, mker(p) = '+'; 
                elseif y(p,rr+1) + y(p,rr) < -cat_bound*my, mker(p) = '*'; end
            end
        elseif categories == -1;
            rp = rand(1);
            if rp<1/2; mker(p) = '+'; else mker(p) = '*'; end
        end

        if timepoints
            plot(y(p,rr),y(p,rr+1),mker(p),'Color',cm(ceil(p*Nc/N),:))
%        plot(y(p,rr),y(p,rr+1),mker(p),'Color',repmat((1-p/N),1,3))
        else
            plot(y(p,rr),y(p,rr+1),mker(p))
        end
        
    end
end

cy = y*cc';

[R,p]=corrcoef(cy)

my = max(abs(cy(:)));

figure(fnum.Number+2),clf,hold on
%plot(cy(:,1),'b-')
%plot(cy(:,2),'r--')
for p=1:N-1
    line([p p+1]',[cy(p,1) cy(p+1,1)]','LineStyle','-','Color',cm(ceil(p*Nc/N),:))
    line([p p+1]',[cy(p,2) cy(p+1,2)]','LineStyle','--','Color',cm(ceil(p*Nc/N),:))
end
axis([0 N -my my])
set(gca,'XTick',0); set(gca,'YTick',0);

