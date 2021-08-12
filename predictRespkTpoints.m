%% Script to predict and evaluate respiration resolved (robust) 3D kT point pulses in the STA regime  
%
% CSA, 2021.01.12
% 
% input: - maps.b1      ... B1 map    [Nx,Ny,Nz,Nc] 
%        - maps.masks   ... ROI       [Nx,Ny,Nz] 
%        - maps.B0 maps ... B0 map    [Nx,Ny,Nz] 
%        - maps.fov     ... FOV in cm [1x3]
%        - wvfrms       ... optimized RF and gradient blips
%        - prbp.dt      ... temporal discretization time in seconds
%        - prbp.nblippts... timepoints of the gradient blip
%        - prbp.Nsubpts ... timepoints of the RF
%
% TBD:  - check again if B0 works

% CV_pre_all = zeros(1,3);
% CV_post_all= zeros(1,3);
function prbp = predictRespkTpoints(maps, wvfrms, prbp)
for offset = 1:3
    numberofmaps = 1;
    % offset=0;
    maps.fov(3)=maps.fov(3);
    maps.numberofmaps = numberofmaps;

    maps.b1   = prbp.rrmaps{(offset)}.b1;
    maps.mask = prbp.rrmaps{(offset)}.mask(end:-1:1,end:-1:1,:);
    maps.b0   = prbp.rrmaps{(offset)}.mask*0;

    optsdef.COLORMAP='thermal';

    %get general parameters
    rfw     = wvfrms.rf; % optimized complex RF weights 
    k       = wvfrms.k;  % optimized k-spce locations
    B1in    = double(maps.b1);   % load B1 maps
    roi     = logical(maps.mask); % load mask to evaluate the B1 prediction
    fov     = maps.fov;  % Field of View in each dim, cm
    f0      = maps.b0(:);% get the B0 map 
    dt      = prbp.dt;   % temporal discretization time in seconds (for B0 only)

    Nc      = size(maps.b1,4); % # Tx channels
    dimxyz  = size(maps.mask); % dimension of design grid
    Ns      = prod(dimxyz);    % total # pixels
    Npulset = size(rfw,1);     % temporal # samples
    gambar  = 4257;            % gamma/2pi in Hz/T
    gam     = gambar*2*pi;     % gamma in radians/g

    %compute the spatial grid using fov in cm
    [xx,yy,zz]=ndgrid(-fov(1)/2:fov(1)/dimxyz(1):fov(1)/2-fov(1)/dimxyz(1), ...
        -fov(2)/2:fov(2)/dimxyz(2):fov(2)/2-fov(2)/dimxyz(2), ...
        -fov(3)/2:fov(3)/dimxyz(3):fov(3)/2-fov(3)/dimxyz(3));
    xx = [xx(:) yy(:) zz(:)];

    %% perform the Bloch prediction in the small tip angle regime
    % create the RF vector and the prediction m1rung for one kT-point 
    rfss = [ones(prbp.Nsubpts,1);zeros(prbp.nblippts,1)];
    Nrp = length(rfss); % number of samples in one rung
    tr = 0:dt:(length(rfss)-1)*dt; % time vector for one rung
    A = 1i*gam*dt*exp(1i*2*pi*f0*tr);
    m1rung = A*rfss;

    % construct design matrix by modifying one rung excitation
    % patterns with appropriate k-space locations, sensitivities, and
    % off-resonance time offset
    sensd = reshape(B1in,[Ns Nc]); 
    A = zeros(Ns,Nc*Npulset);      %initialize A 
        for ii = 1:Npulset
            % blip-induced phase shift
            kphs = xx*k(ii,:)';
            % off res-induced phase shift - account for phase accrual to
            % end of pulse
            totphs = exp(1i*2*pi*(f0*((ii-1)*Nrp - Npulset*Nrp)*dt+kphs));
            tmp = m1rung.*totphs;
            for kk = 1:Nc
                % apply sens, stick it in the design matrix
                A(:,(kk-1)*Npulset+ii) = sensd(:,kk).*tmp;
            end
        end

    % get excitation patterns and update target phase
    m = reshape(A * rfw(:),[Ns 1]);
    phs = angle(m);

    images = (reshape(m,size(maps.mask)))*180/pi; %get the prediction in degree
    % images=sum(abs(B1in),4);

    %% compare the prediction to the default state and compute the CV
    % using mainly code from show_shim_prediction_TXfct.m

    %get the B1 map (before/after)
    b1pat_pre = (sum(maps.b1,4)); %get the B1 map (default state)
    %scale the default state to achieve a median of the desired FA in the ROI
    b1pat_pre_scaling=1/median(abs(b1pat_pre(roi)))*prbp.delta_tip;
    b1pat_pre = b1pat_pre*b1pat_pre_scaling;
    b1pat_post = images; %get the B1 prediction (optimized state)
    b1pat_both = cat(3,b1pat_pre,b1pat_post);

    %compute the efficiency
    b1sumofmag = sum(abs(maps.b1),4)*b1pat_pre_scaling; %compute the SOS 
    Eff_pre = abs(b1pat_pre)./b1sumofmag;
    Eff_post = abs(b1pat_post)./b1sumofmag;
    Eff_both = cat(3,Eff_pre,Eff_post);

    %plot the prediction before and after shimming
    roi_bot = cat(3,roi,roi);

    tmp_pre = abs(b1pat_pre(~~roi));
    % min(tmp_pre)
    tmp_post = abs(b1pat_post(~~roi));

    % figure
    % nhist(tmp_post,'samebins','numbers','minx',2.5,'maxx',25,'binfactor',0.5);

    CV_pre = std(tmp_pre(:))/mean(tmp_pre(:));
    CV_post = std(tmp_post(:))/mean(tmp_post(:));

    MeanEff_pre = mean(Eff_pre(~~roi));
    MeanEff_post = mean(Eff_post(~~roi));

    lNoOfSlices = size(maps.b1,3);
    
    prbp.CV_post_all(prbp.c_datasets,prbp.respstate+1,offset) = CV_post;
    
end
