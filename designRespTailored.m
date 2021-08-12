% This fuprbp.Nction computes tailored kT-point pulses as described in
% Christoph S. Aigner, Sebastian Dietrich, Tobias Schaeffter and Sebastian
% Schmitter, Calibration-free pTx of the human heart at 7T via 3D universal 
% pulses, submitted to Magn. Reson. Med. 2021
%
% The optimization of the kT-points is performed using code by Zhipeng Cao 
% and Will Grissom (https://bitbucket.org/wgrissom/acptx/) who have given 
% permission for iprbp.Nclusion within this package. Please cite appropriately.
% 
% usage: wvfrms = designTailored(pulseType, numkTpoints, numPhaseInit, 
%             lambdavec, phsinitmode, b_evalAllDatasets, numTailored, prbp) 
%                   pulseType   = string that defines the pulse type
%                   numkTpoints = array of kT points
%                   numPhaseInit= array of phase indices
%                   lambdavec   = array of regularization parameters
%                   phsinitmode = phase init for the optimization 
%                   b_evalAllDatasets  = flag for evaluation of the results
%                   numTailored = just compute one tailored pulse
%                   prbp        = struct with most problem related params
%                   wvfrms      = struct that contains the waveforms
% pulseType can be 'default', 'tailored' or 'UP' and phsinitmode can be
% 'zerophase', 'defaultphase', 'quadmode' and 'randphase'.
% wvfrms contains the variables wvfrms.k and wvfrms.rf. 
%
% Created by Christoph S. Aigner, PTB, June 2021.
% Email: christoph.aigner@ptb.de

function [wvfrms, prbp] = designRespTailored(pulseType, numkTpoints, numPhaseInit, lambdavec, phsinitmode, b_evalAllDatasets, numTailored, prbp) 
    load('kTrandphases.mat');
        
    numberofmaps = 1; % how many different B1 maps are used in the optimization
    rrespstate = prbp.respstate;
    if prbp.respstate == 3 % resp robust
        numberofmaps = 3;
        rrespstate = 0;
    end
    
    for c_kTpoints = numkTpoints
        for c_diffrand = numPhaseInit
            for c_lambdaexp=1:length(lambdavec)
                for c_dat = 1:numTailored
                    disp(['design ',pulseType,num2str(c_dat),'-',...
                          num2str(c_kTpoints),'kT, phaseinit=',...
                          num2str(c_diffrand), ', regularization=',...
                          num2str(lambdavec(c_lambdaexp))]);

                    Nm = 1; %tailored design is for 1 B1+ dataset only
                    prbp.Nm = Nm;

                    % initialize maps (to be sure to have the right size)
                    Nx = size(prbp.rrmaps{1}.b1,1); % # x
                    Ny = size(prbp.rrmaps{1}.b1,2); % # y
                    Nz = size(prbp.rrmaps{1}.b1,3); % # slices
                    prbp.Nc = size(prbp.rrmaps{1}.b1,4); % number of physical coils (Ncoils)
                    maps.numberofmaps = Nm;
                    maps.fov = [ 31.2500   31.2500   25.0000];
                    maps.numberofmaps = numberofmaps;

                    maps.b1   = zeros(Nx,Ny,Nz*numberofmaps,prbp.Nc);
                    maps.mask = zeros(Nx,Ny,Nz*numberofmaps);
                    maps.b0   = zeros(Nx,Ny,Nz*numberofmaps);

                    for count_resp=1:numberofmaps
                        maps.b1(:,:,(1:Nz)+(count_resp-1)*Nz,:)= prbp.rrmaps{count_resp+rrespstate}.b1(:,:,1:end,:);
                        maps.mask(:,:,(1:Nz)+(count_resp-1)*Nz)= prbp.rrmaps{count_resp+rrespstate}.mask(end:-1:1,end:-1:1,1:end);
                        maps.b0(:,:,(1:Nz)+(count_resp-1)*Nz)  = prbp.rrmaps{count_resp+rrespstate}.mask(end:-1:1,end:-1:1,1:end)*0;
                    end

                    maps.mask = logical(maps.mask);

                     % set initial target phase to zero or default phase mode 
                    switch phsinitmode
                        case 'zerophase'
                            disp('zero phase initial')
                            maps.phsinit = zeros(size(maps.mask)); 
                        case 'defaultphase'
                            disp('default phase initial')
                            maps.phsinit = angle(sum(maps.b1,4));%default phase
                        case 'quadmode' %quad mode does not perform in the body
                            bcb1 = 0;
                            for ii = 1:prbp.Nc
                               bcb1 = bcb1 + maps.b1(:,:,:,ii)*...
                                             exp(1i*(ii-1)*2*pi/prbp.Nc).*...
                                             exp(-1i*angle(maps.b1(:,:,:,1)));
                            end
                            maps.phsinit = angle(bcb1);
                        case 'randphase'
                            bcb1 = 0;
                            for ii = 1:prbp.Nc
                                bcb1 = bcb1 + maps.b1(:,:,:,ii)*...
                                             exp(1i*randphases(c_diffrand,ii)).*...
                                             exp(-1i*angle(maps.b1(:,:,:,1)));
                            end
                            maps.phsinit = angle(bcb1);
                    end

                    % Algorithm and problem parameters
                    prbp.ndims = ndims(maps.mask);     % # spatial dimensions 
                    prbp.kmaxdistance = [Inf Inf Inf]; % maximum kT-point location
                    if lambdavec(c_lambdaexp)==0
                    	prbp.beta = 10^-6;                 % initial RF regularization
                        prbp.betaadjust = 1;               % automatically adjust RF regularization parameter                        
                    else
                    	prbp.beta = lambdavec(c_lambdaexp);% initial RF regularization
                        prbp.betaadjust = 0;               % automatically adjust RF regularization parameter
                    end


                    prbp.dimxyz = size(maps.b0);
                    prbp.filtertype = 'Lin phase';     % alternately add kT-points on either size of the DC point
                    prbp.trajres = 2;                  % maximum spatial frequency of OMP search grid (was 2 for most August 2015 results)
                    prbp.Npulse = c_kTpoints;          % number of kT-points subpulses
                    algp.nthreads  = 10;               % number of compute threads (for mex only)
                    algp.computemethod = 'mex';
                    algp.prbp.Ncgiters = 3; 
                    algp.cgtol = 0.9999;
                    prbp.Ncred = inf;

                    % Run the kT point design 
                    %   m ... STA solution 
                    %   wvfrms ... optimized RF and gradient blips)        
                    [all_m, wvfrms] = dzktpts(algp,prbp,maps); 
                    rfw = wvfrms.rf;

                    %evaluate the optimized results
                    farmse = sqrt(mean((abs(all_m.images(maps.mask))/pi*180 - prbp.delta_tip).^2));
                    rfrms = norm(rfw);
                    fprintf('Flip angle RMSE: %.4f, RMS RF power: %.4f.\n\n',farmse,rfrms);

                    %save the results for later
                    farmse_all(c_kTpoints, c_diffrand, c_lambdaexp)    = farmse;
                    rfrms_all(c_kTpoints, c_diffrand, c_lambdaexp)     = rfrms;
                    waveforms_all{c_kTpoints, c_diffrand, c_lambdaexp} = wvfrms;

                    if b_evalAllDatasets == true
%                         evalAllDatasets('tailored',wvfrms, numTailored, numkTpoints, c_diffrand, c_lambdaexp, lambdavec, c_kTpoints, prbp); drawnow;
                        prbp = predictRespkTpoints(maps, wvfrms, prbp);
                                   
                    end
                end
            end

            %plot l-curve if more than one runs with different beta were done
            if (c_lambdaexp == length(lambdavec)) && c_lambdaexp > 1
                figure;hold all;
                loglog(squeeze(farmse_all(c_kTpoints,:,1:end)).',squeeze(rfrms_all(c_kTpoints,:,1:end)).');
                xlabel('log_{10}(FA RMSE / deg)');
                ylabel('log_{10}(RF RMSE / a.u.)');
                sgtitle (['L-curve for ',pulseType,num2str(c_dat),'-',num2str(c_kTpoints),'kT phaseinit=',num2str(numPhaseInit)]);
            end
        end
    end
end