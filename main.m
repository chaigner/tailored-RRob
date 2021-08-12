% This script loads upt to 10 channel-wise invivo 3D B1+ datasets of the 
% human body during deep breathing at 7T and computes and evaluates 
% tailored non-selective respiration specific and respiration robust 
% kT-points pTx pulses in the human heart.
%
% Christoph S. Aigner, Sebastian Dietrich, Tobias Schaeffter and Sebastian
% Schmitter, Respiration induced B1+ changes and their impact on universal 
% and tailored 3D kT point pulses in 7T body imaging, submitted to Magn. 
% Reson. Med. 2021
%
% The 10 deep breathing channel-wise invivo B1+ datasets of the human body 
% at 7T are available at: 
% TBD
%
% The optimization of the kT-points is performed using code by Will Grissom
% and Zhipeng Cao (https://bitbucket.org/wgrissom/acptx/) who have given 
% permission for inclusion within this package. Please cite appropriately.
% 
% Created by Christoph S. Aigner, PTB, Agust 2021.
% Email: christoph.aigner@ptb.de
%
% This code is free under the terms of the GPL-3.0 license.

addpath ktutil % add the directory of the source code 
load('cmap.mat');
load('kTrandphases.mat');

% parameters
% allIndices     ... all B1+ datasets (library + test-cases)
% libraryIndices ... datasets used in the optimization (library)
prbp.pathDat        = 'RR_DB_B1R'; % set the folder that contains the in vivo B1+ datasets
% prbp.allIndices     = 1:31; % there are 31 B1+ datasets in total
% prbp.libraryIndices = 1:22; % the paper used the first 22 B1+ datasets
% prbp.allmaps        = cell(1, length(prbp.allIndices));     % pre-allocate the cell 
% prbp.librarymaps    = cell(1, length(prbp.libraryIndices)); % pre-allocate the cell 
prbp.dt             = 10e-6; % dwell time in sec
prbp.Nsubpts        = 10;    % # of time points for RF subpulses 
prbp.nblippts       = 20;    % # of time points for gradient blips
prbp.delta_tip      = 10;    % flip angle in degrees

datasets = 1:10;

%% perform the respiration resolved (robust) 3D kT point pulses 


for c_datasets = 1:length(datasets) %if you want to loop over all datas
% for c_datasets = 8
    
    load(['RR_DB_B1R\RRB1R' num2str(datasets(c_datasets)),'.mat'])
    prbp.rrmaps = rrmaps;
    
    for respstate = 0:3
        % respstate = 0,1 or 2 ... respiration resolved 3D kT point design
        % respstate = 3 ... respiration resolved 3D kT point design

        prbp.respstate = respstate;
        prbp.c_datasets = c_datasets;

        % tailored design
        % default parameters:
        %   numkTpoints       = 4;      % number of kT points; tested for 1:5
        %   numPhaseInit      = 165;    % 1-200, #165 performed best for the library
        %   lambdavec         = [4.64];  % result of the L curve optimization; 10^0-10^7
        %   phsinitmode       = 'randphase'; % performed best
        %   b_evalAllDatasets = true;   % evaluate the tailored pulse in allmaps
        %   numTailored       = 1;      % just compute one tailored pulse

        %do the tailored design and evaluate the tailored pulse in allmaps
        [wvfrms, prbp] = designRespTailored('tailored', 4, 165, 0, 'randphase', 1, 1, prbp); 
    end

    % plot the results
    figure;
    for c_pulse = 1:4
        plot(squeeze(prbp.CV_post_all(c_datasets,c_pulse,:))); hold all;
        h = plot(squeeze(prbp.CV_post_all(c_datasets,c_pulse,:)),'k*'); hold all;
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    end
    set(gca, ...
      'Box'         , 'off'     , ...
      'TickDir'     , 'out'     , ...
      'TickLength'  , [.02 .02] , ...
      'XMinorTick'  , 'off'      , ...
      'YMinorTick'  , 'on'      , ...
      'YGrid'       , 'on'      , ...
      'XColor'      , [.3 .3 .3], ...
      'YColor'      , [.3 .3 .3], ...
      'YTick'       , 0:0.025:1, ...
      'LineWidth'   , 2         );

    xticks([1 2 3])
    xticklabels({'exhale','intermediate','inhale'})
    xlabel('respiration resolved DB B1+ maps');
    ylabel('CV in % (ROI)');
    legend({'tailored-RSpec-exhale','tailored-RSpec-intermediate','tailored-RSpec-inhale','tailored-RRob'})
    title(['Coeficient of Variation of the flip angle map in the 3D heart ROI in subject ',num2str(c_datasets)])
            
end
