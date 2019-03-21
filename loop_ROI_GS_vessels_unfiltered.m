%% Second test of multiscale PID (Partial Information Decomposition)
% analysis at varying scale tau 
% uses the PID code from here 
% Faes, L., Marinazzo, D., & Stramaglia, S. (2017).
% Multiscale Information Decomposition: Exact Computation for Multivariate Gaussian Processes.
% Entropy, 19(8), 408. doi:10.3390/e19080408

clear; close all; clc;
subs=dir('/home/daniele/unfiltered_timeseries/'); %change this
subs(1:2)=[];
nsubs=length(subs);
nrois=278;
% Parameters
maxscale=30;
tauv=(1:maxscale)';
nscales=length(tauv);
ncoeff=nscales; % if FIR, set the number of coeff


Til_j=zeros(nsubs,nrois,nscales); %Joint transfer
Ti_j=zeros(nsubs,nrois,nscales); % Individual transfer
Tl_j=zeros(nsubs,nrois,nscales);
Ij_il=zeros(nsubs,nrois,nscales); %Interaction transfer (NET SYNERGY)

% Partial Information Decomposition
Ril_j=zeros(nsubs,nrois,nscales); % Redundant transfer
Ui_j=zeros(nsubs,nrois,nscales); % Unique transfer
Ul_j=zeros(nsubs,nrois,nscales); % Unique transfer
Sil_j=zeros(nsubs,nrois,nscales); %Synergistic transfer

% model order
pmax=10; % pmax to scan for BIC criterion

Yo=zeros(3,1200); %three time series (two drivers, one target), and lenght of time series (1200 for HCP data).
jj=1; % index of target
ii=2; % index of first driver
kk=3; % index of second driver
%% MSTE Analysis
for isubs=nsubs:-1:1
    tic;
    disp(isubs);
    load(['/home/daniele/unfiltered_timeseries/' subs(isubs).name '/timeseries.mat']);
    Yo(3,:)=double(timeseries(279,:)); % Global Signal
    Yo(2,:)=double(timeseries(280,:)); % BOLD in vessels
    for iroi=1:nrois
        Yo(1,:)=double(timeseries(iroi,:));
        
        
        [M,N]=size(Yo);
        Y=zscore(Yo,[],2);
        
        %%%%% identification
        [p_aic,p_bic,aic,bic] = mos_idMVAR(Y,pmax,0); % finds model order
        p=p_bic; % use bayesian Information Criterion
        %figure; plot(bic)
        [Am,Su]=idMVAR(Y,p,0);
        
        % Stability check
        E=eye(M*p);AA=[Am;E(1:end-M,:)];lambda=eig(AA);lambdamaxo=max(abs(lambda));
        if lambdamaxo>=1
            warning('Non-stable VAR process');
        end
        
        parfor s=1:nscales
            %disp(['scale ' int2str(s)]);
            tau=tauv(s);
            
            % MA parameters resulting from the change of scale
            if tau==1
                q=0; b=1;
            else
                q=ncoeff; % number of filter coeffs
                ft=1/(2*tau); %cutoff frequency
                Wn=2*ft; %normalized cutoff frequency (fNyquist=1)
                b=fir1(q,Wn,'noscale'); %Hamming window, linear phase (symmetry of b coeffs)
            end
            Bk=zeros(M,M,q+1);
            for l=1:q+1
                Bk(:,:,l)=b(l)*eye(M);
            end
            Bm=[];
            for kp=1:q+1
                Bm=[Bm Bk(:,:,kp)];
            end
            B0=Bm(1:M,1:M);
            Bm=Bm(1:M,M+1:end);
            
            
            % ISS parameters
            [A,C,K,V,Vy] = varma2iss(Am,Bm,Su,B0); % max(abs(eig(A-K*C)))
            
            %%% parameters after downsampling
            [Ad,Kd,Vd] = iss_ds(A,C,K,V,tau);
            Cd=C; %Rtau=R;
            
            [VR, lambda0] = iss_PCOV(Ad,Cd,Kd,Vd,jj);
            Sj=lambda0(jj,jj);
            Sj_j=VR;
            
            tmp = iss_PCOV(Ad,Cd,Kd,Vd,[jj ii]);
            Sj_ji=tmp(1,1);
            
            tmp = iss_PCOV(Ad,Cd,Kd,Vd,[jj kk]);
            Sj_jl=tmp(1,1);
            
            tmp = iss_PCOV(Ad,Cd,Kd,Vd,[jj ii kk]);
            Sj_ijl=tmp(1,1);
            
            % Interaction Information Decomposition
            Til_j(isubs,iroi,s)=0.5*log(Sj_j/Sj_ijl); %Joint transfer
            Ti_j(isubs,iroi,s)=0.5*log(Sj_j/Sj_ji); % Individual transfer
            Tl_j(isubs,iroi,s)=0.5*log(Sj_j/Sj_jl);
            Ij_il(isubs,iroi,s)=-Ti_j(isubs,iroi,s)-Tl_j(isubs,iroi,s)+Til_j(isubs,iroi,s); %Interaction transfer (NET SYNERGY)
            
            % Partial Information Decomposition
            Ril_j(isubs,iroi,s)=min(Ti_j(isubs,iroi,s),Tl_j(isubs,iroi,s)); % Redundant transfer
            Ui_j(isubs,iroi,s)=Ti_j(isubs,iroi,s)-Ril_j(isubs,iroi,s); % Unique transfer
            Ul_j(isubs,iroi,s)=Tl_j(isubs,iroi,s)-Ril_j(isubs,iroi,s); % Unique transfer
            Sil_j(isubs,iroi,s)=Til_j(isubs,iroi,s)-Ui_j(isubs,iroi,s)-Ul_j(isubs,iroi,s)-Ril_j(isubs,iroi,s); %Synergistic transfer
        end
    end
    toc;
end
save MSID_ROI_GS_Vessels_unfiltered_par10 Til_j Ti_j Tl_j Ij_il Ril_j Ui_j Ul_j Sil_j