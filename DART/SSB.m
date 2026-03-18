classdef SSB
    methods (Static)
        
        function varargout = Estimate_BM4model(DATA,IDXp2bu,SAT,WeightingScheme)
            
            %% Settings
            defval('BSsig0',0.1) %Bin size used to bin sigma0 values [dB]
            defval('BSswh',0.25) %Bin size used to bin SWH values [m]
            
            defval('MakePlots',false);
            defval('WeightingScheme','BasedOnSEM')
            
            %% Exclude points not to be used in estimating the SSB model
            DATA = structfun(@(x) double(x(IDXp2bu,:)), DATA, 'UniformOutput', false);
            
            %% Bin the observations by sigma0 and SWH
            [~,EDGc,EDGr,IDXc,IDXr] = histcounts2(DATA.sigma0,DATA.SWH,floor(min(DATA.sigma0)):BSsig0:ceil(max(DATA.sigma0)),floor(min(DATA.SWH)):BSswh:ceil(max(DATA.SWH)));
            IDXlin                  = sub2ind([numel(EDGr)-1 numel(EDGc)-1],IDXr,IDXc);
            [BinX,BinY]             = meshgrid(EDGc(1:end-1) + 0.5*BSswh,EDGr(1:end-1) + 0.5*BSsig0);
            
            %% Compute summary statistics by IDXlin
            [BIN.IDX,BIN.Mean,BIN.NrP,BIN.RMS,BIN.SD] = grpstats(DATA.SSHi-DATA.CorrDT-DATA.DTU18MSL-DATA.SHcorr,IDXlin,{'gname','mean','numel','rms','std'});
            BIN.IDX                                   = str2double(BIN.IDX);
            BIN.Sig0                                  = BinX(BIN.IDX);
            BIN.SWH                                   = BinY(BIN.IDX);
            BIN.SEM                                   = BIN.SD./sqrt(BIN.NrP); %standard error of the mean
            
            %% Find cells for which number of points used to compute mean >= 10 AND RMS >= 0.03 (to avoid over-weighting)
            BIN = structfun(@(x) (x(BIN.NrP >= 10 & BIN.RMS >= 0.03,:)), BIN, 'UniformOutput', false);
            
            if MakePlots
                % CreateScatterPlot(DATA.sigma0,DATA.SWH,'Sigma0 [dB]','Significant Wave Height [m]',sprintf('Sigma0 versus SWH for %s',SAT),0.75)
                
                DUM          = nan(size(BinX));
                DUM(BIN.IDX) = BIN.Mean;
                figure,imagescnan(BinX(1,:),BinY(:,1),DUM,'NanColor',[.5 .5 .5]),hold on,caxis([-0.4 0.4]),axis xy
                colorbar
                title(sprintf('Direct estimate SBB [meter] for mission %s',SAT))
                xlabel('Sigma0 [dB]')
                ylabel('Significant Wave Height [m]')
                
                DUM(BIN.IDX) = BIN.NrP;
                figure,imagescnan(BinX(1,:),BinY(:,1),DUM,'NanColor',[.5 .5 .5]),hold on,caxis([0 1000]),axis xy
                colorbar
                title(sprintf('Number of points used to compute mean for mission %s',SAT))
                xlabel('Sigma0 [dB]')
                ylabel('Significant Wave Height [m]')
                
                DUM(BIN.IDX) = BIN.RMS;
                figure,imagescnan(BinX(1,:),BinY(:,1),DUM,'NanColor',[.5 .5 .5]),hold on,caxis([0 0.4]),axis xy
                colorbar
                title(sprintf('RMS [meter] for mission %s',SAT))
                xlabel('Sigma0 [dB]')
                ylabel('Significant Wave Height [m]')
            end
            
            %% Estimate BM4 model parameters
            A           = [ones(numel(BIN.IDX),1),BIN.SWH,BIN.SWH.^2,BIN.SWH.*BIN.Sig0,BIN.SWH.*BIN.Sig0.^2];
            if strcmpi(WeightingScheme,'BasedOnSEM')
                WEIGHTS = 1./BIN.SEM;
            elseif strcmpi(WeightingScheme,'UnitWeights')
                WEIGHTS = ones(numel(BIN.IDX),1);
            end
            x_hat       = lscov(A,BIN.Mean,WEIGHTS);
            
            if MakePlots
                DUM          = nan(size(BinX));
                DUM(BIN.IDX) = WEIGHTS./max(WEIGHTS);
                figure,imagescnan(BinX(1,:),BinY(:,1),DUM,'NanColor',[.5 .5 .5]),hold on,caxis([0 1]),axis xy
                colorbar
                title(sprintf('Normalized weights for mission %s using %s',SAT,WeightingScheme))
                xlabel('Sigma0 [dB]')
                ylabel('Significant Wave Height [m]')
                
                As = [zeros(numel(BinX),1),BinY(:),BinY(:).^2,BinY(:).*BinX(:),BinY(:).*BinX(:).^2];
                figure,imagesc(BinX(1,:),BinY(:,1),reshape(As*x_hat,size(BinX))),hold on,caxis([-0.4 0.4]),axis xy
                colorbar
                title(sprintf('Modelled SSB [meter] for mission %s using %s, note that bias is removed!',SAT,WeightingScheme))
                xlabel('Sigma0 [dB]')
                ylabel('Significant Wave Height [m]')
            end
            
            %The approach of Scharroo will result in problems at the edges of the
            %SWH/Vwnd domain, disregard them and use estimate BM4-model
            % GKernel = fspecial('gaussian', size(XI), 1);
            % TEST    = reshape(Ysynth,size(XI));
            % A       = [ones(size(BE,1),1),zeros(size(BE,1),4)];
            % TEST(sub2ind(size(XI),SortBINIDs(BE(:,2),2),SortBINIDs(BE(:,1),1))) = BE(:,5)-A*x_hat;
            % % would generate a Gaussian low-pass filter of size 3X3 and Sigma =1.
            % % Then you can filter your image with 'imfilter' function:
            % J = imfilter(TEST, GKernel);
            % disp('Note that Scharroo uses additional weighting when computing Gaussian smoothed results')
            
            As        = [zeros(numel(DATA.TIME),1),DATA.SWH,DATA.SWH.^2,DATA.SWH.*DATA.sigma0,DATA.SWH.*DATA.sigma0.^2];
            SSB       = As*x_hat;
            
            varns     = {x_hat,SSB};
            varargout = varns(1:nargout);
            
        end
        
        function SSB = ComputeSSBcorr(DATA,IDXp2bu,x_hat)
            
            %% Exclude points not to be used in estimating the SSB model
            DATA = structfun(@(x) double(x(IDXp2bu,:)), DATA, 'UniformOutput', false);

            %% Compute SSB corrections
            As   = [zeros(numel(DATA.TIME),1),DATA.SWH,DATA.SWH.^2,DATA.SWH.*DATA.sigma0,DATA.SWH.*DATA.sigma0.^2];
            SSB  = As*x_hat;

        end
    end
end