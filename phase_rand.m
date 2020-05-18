function x_null = phase_rand(x, Nscram, ComScram)

% function [truecorr, pvals, nullcorrs] = PHASE_RAND_CORR(x,y,nscram, tail)
%
% This function calculates the correlation between 
%     [1] a phase-scrambled version of each column of 'x'
% and 
%     [2] the intact vector 'y'
% to produce a distribution of null correlations in which we have controlled
% for the power spectrum (and thus temporal autocorrelation) of the input time-series.
%
% INPUT
% x =        [Nsamp by K] matrix, each of whose K columns will be phase-scrambled,
%             and correlated against the vector input 'y'.
% y =        [Nsamp by 1] input vector [will be left intact]
% Nscram =   [integer] number of phase-scramble correlation value to compute (default: 1000)
% tail =     flag indicating distributional tail to examine for stats:
%            -1 --> left-tail, 0 --> two-tail, +1 --> right-tailed
% permutation = [boolean] true: permute existing phases, false: generate random phases
%
% OUTPUT
% truecorr = [1 by K] vector of Pearson correlation between intact 'x' and intact 'y'
% pvals =    [1 by K] vector of corresponding p-values for the values of truecorr
%                  based on comparison with null distributions 
% nullcorrs = [Nscram by K] matrix of "null" Pearson correlations between
%                   phase-scrambled columns of 'x and intact y
%
% TODO: 
%   implement some tapering to avoid high-freq artifacts in the fft
%   (but the procedure could be potentially probematic)
%   
%
% Author: CJ Honey
% Version: 0.1, April 2010  (phase_scram_corr_Nvs1)
% Version: 0.2, March 2011    -- fixed bug with zero-th phase component;
%                             -- randomizes rather than scrambles phase 
%                             based on feedback from Jochen Weber
%

if nargin < 3; ComScram = 0; end
if nargin < 2; Nscram = 1000; end

[Nsamp K] = size(x);  %extract number of samples and number of signals

x = x - repmat(mean(x,1), Nsamp,1);  %remove the mean of each column of X
x = x./sqrt(repmat(dot(x,x), Nsamp, 1)/(Nsamp-1)); %divide by the standard deviation of each column

%transform the vectors-to-be-scrambled to the frequency domain
Fx = fft(x); 

% identify indices of positive and negative frequency components of the fft
% we need to know these so that we can symmetrize phase of neg and pos freq
if mod(Nsamp,2) == 0
    posfreqs = 2:(Nsamp/2);
    negfreqs = Nsamp : -1 : (Nsamp/2)+2;
else
    posfreqs = 2:(Nsamp+1)/2;
    negfreqs = Nsamp : -1 : (Nsamp+1)/2 + 1;
end

x_amp = abs(Fx);  %get the amplitude of the Fourier components

x_phase = atan2(imag(Fx), real(Fx)); %get the phases of the Fourier components [NB: must use 'atan2', not 'atan' to get the sign of the angle right]

J = sqrt(-1);  %define the vertical vector in the complex plane

sym_phase = zeros(Nsamp,K);   %will contain symmetrized randomized phases for each bootstrap


for n = 1:Nscram

    if ComScram
        r = repmat(rand(Nsamp,1),1,K);
    else
        r = rand(Nsamp,K);
    end
    [~,rp] = sort(r);

    x_phase=x_phase(rp);
    sym_phase(posfreqs,:) = x_phase(1:length(posfreqs),:);
    sym_phase(negfreqs,:) = -x_phase(1:length(posfreqs),:);
    
    
    z = x_amp.*exp(J.*sym_phase); %generate (symmetric)-phase-scrambled Fourier components
    x_null(:,:,n) = ifft(z); %invert the fft to generate a phase-scrambled version of x
end






