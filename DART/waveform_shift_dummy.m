function [CS1b] = waveform_shift_dummy(CS1b,n_multilook,mission)
    % function to shift neighbouring FF-SAR waveforms onto the same tracker
    % distances before multilooking, tryout implementation, last
    % n<n_multilooks waveform not being aligned
    [CONST,TR,FFSAR_processing_settings] = FFSAR_LoadCommonSettings;
    DDAcf   = DDA_ConfigFile(mission,'SAR');
    
    Nmax = numel(CS1b.MEA.win_delay);
    n = n_multilook;
    tau = CS1b.MEA.win_delay;
    Rtrk = tau*CONST.c/2;
    wav = CS1b.SAR.dataIQ;

    DDAcf.Ns = size(CS1b.SAR.dataIQ,1);        %Nr of bins/samples in any waveform
    os_ZP    = DDAcf.Ns/DDAcf.Np;          %Zero-Padding Oversampling Factor
    % RefBin   = (DDAcf.RefBin-1)*os_ZP + 1;
    dt       = 1/(os_ZP*DDAcf.B);          %Waveform sampling interval [s]
    dr       = dt*CONST.c/2; 

    

    for i = 1:floor(Nmax/n_multilook)
        
        wav_ind = (n*(i-1)+1):n*i;
        center_wav_ind = (i-1)*n + n/2;
        
        % ###########################################
        if false%i==floor(Nmax/n_multilook/2)-5
            figure;
            plot(Rtrk - Rtrk(1));

            figure;
            subplot(1,3,1);
            imagesc(abs(wav(:,wav_ind)').^2);

            subplot(1,3,2);
            mean_wav_before = mean(abs(wav(:,wav_ind)).^2,2);
            plot(mean_wav_before);

            subplot(1,3,3);
            [maxval,ind] = max(abs(wav(:,wav_ind)).^2,[],1);
            plot(ind,wav_ind);hold on
            for j = 1:4
                scatter(median(ind((j-1)*n/4+1:j*n/4)),wav_ind(1)+j*n/4-n/8,'ro');hold on
            end
            xlim([40,50])
            grid on;
        end
        % ###########################################
        bin_shifts = (Rtrk(wav_ind)-Rtrk(center_wav_ind))/dr;      % take the normalized differences in tracker ranges along the waveforms [#bins]
        shift_phasor = exp(-2*pi*1i*((1:DDAcf.Ns)-1)'/DDAcf.Ns*bin_shifts);

        wav(:,wav_ind) = fft(wav(:,wav_ind),[],1);
        wav(:,wav_ind) = wav(:,wav_ind).*shift_phasor;
        wav(:,wav_ind) = ifft(wav(:,wav_ind),[],1);
        
        CS1b.MEA.win_delay(wav_ind) = Rtrk(center_wav_ind)*2/CONST.c;
        CS1b.MEA.ref_range(wav_ind) = Rtrk(center_wav_ind)*2/CONST.c;
        CS1b.SAR.data(:,wav_ind) = abs(wav(:,wav_ind)).^2;
        
        % ###########################################    
        if false%i==floor(Nmax/n_multilook/2)-5
            figure;
            subplot(1,3,1);
            imagesc(abs(wav(:,wav_ind)').^2);

            subplot(1,3,2);
            mean_wav_after = mean(abs(wav(:,wav_ind)).^2,2);
            plot(mean_wav_before); hold on
            plot(mean_wav_after);
            grid on;
            legend('before', 'after')

            subplot(1,3,3);
            [maxval,ind] = max(abs(wav(:,wav_ind)).^2,[],1);
            plot(ind,wav_ind);hold on
            for j = 1:4
                scatter(median(ind((j-1)*n/4+1:j*n/4)),wav_ind(1)+j*n/4-n/8,'ro');hold on
            end
            
            xlim([40,50])
            grid on;
        end
    end
    %%
    % bin_shift = 3;
    % wav_double = wav(:,1:2);
    % 
    % figure;
    % plot(abs(wav_double(:,1)).^2); hold on
    % 
    % wav_double_fft = fft(wav_double,[],1);
    % 
    % % positive shift = shifting towards higher ranges
    % % negative shift = shifting towards lower ranges
    % 
    % range_corr = exp(-2*pi*1i*((1:DDAcf.Ns)-1)'/DDAcf.Ns*bin_shift );
    % 
    % wav_double_fft = wav_double_fft.*range_corr;
    % 
    % wav_double_shifted = ifft(wav_double_fft,[],1);
    % 
    % 
    % plot(abs(wav_double_shifted(:,1)).^2); hold on
    % legend('original','shifted')
    % % plot(imag(wav_double_shifted(:,1))); hold on
    % 
    % % figure;
    % % plot(real(wav_double_shifted(3:end,1)) - wav_double(1:end-2,1)); hold on
    % 
    % sum(abs(wav_double_shifted(:,1)).^2)
    % 
    % 
    % %% undo shifting
    % 
    % wav_double_fft = fft(wav_double_shifted,[],1);
    % 
    % % positive shift = shifting towards higher ranges
    % % negative shift = shifting towards lower ranges
    % 
    % range_corr = exp(-2*pi*1i*((1:DDAcf.Ns)-1)'/DDAcf.Ns*-bin_shift );
    % 
    % wav_double_fft = wav_double_fft.*range_corr;
    % 
    % wav_double_shifted = ifft(wav_double_fft,[],1);
    % 
    % plot(abs(wav_double_shifted(:,1)).^2,'bo'); hold on
    % legend('unshifted')
    % 
    % sum(abs(wav_double_shifted(:,1)).^2)