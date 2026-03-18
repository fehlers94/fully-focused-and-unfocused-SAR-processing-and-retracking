function CS1b = Compute_MultiLooked_Waveforms(CS1b)

%COMPUTE_MULTILOOKED_WAVEFORMS computes the multi-looked waveform at each
%surface location

%Input:
%CS1b: range compressed stacks

%Output:
%CS1b.S**.data:              multi-looked waveforms
%CS1b.S**.N_averaged_echoes: # of averaged echoes

%% Compute multi-looked waveforms
switch CS1b.GEO.OPERATION_MODE
    case 'SIR_L1B_SAR'
        for i=find(~isnan(CS1b.GEO.LAT))'
            CS1b.SAR.data(:,i)            = nanmean(CS1b.SAR.Stack(:,CS1b.GEO.IDXstack{i}),2);
            CS1b.SAR.N_averaged_echoes(i) = numel(CS1b.GEO.IDXstack{i});
        end
    case 'SIR_L1B_SARIN'
        for i=find(~isnan(CS1b.GEO.LAT))'
            %Average over stack per channel
            data_ch1                      = nanmean(CS1b.SIN.Stack_ch1(:,CS1b.GEO.IDXstack{i}),2);
            data_ch2                      = nanmean(CS1b.SIN.Stack_ch2(:,CS1b.GEO.IDXstack{i}),2);
            
            %Average over channels
            CS1b.SIN.data(:,i)            = nanmean([data_ch1 data_ch2],2);
            
            CS1b.SIN.N_averaged_echoes(i) = numel(CS1b.GEO.IDXstack{i});
        end
    otherwise
        error('OPERATION_MODE: %s not recognized')
end

end
    
 