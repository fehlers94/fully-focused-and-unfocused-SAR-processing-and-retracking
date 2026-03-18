function retr_corr = threshold_retracker(wav,ref_bin,binsize)
    threshold = 0.75;

    % go through all waveforms
    retr_corr=zeros(length(wav(1,:)),1);
    for i=1:length(wav(1,:))

        % find first value above 0.5*max(power)
        q=find(wav(:,i) > threshold*max(wav(:,i)));
        p=q(1)-1; % first index with > 0.5*max(power)
        q=q(1); % last index below < 0.5*max(power)
        
        % get retracker correction
        % this gets the 'decimal bin' at 0.5*max(power)
        if p==0
            retr_corr(i)=NaN;
        else
            retr_bin=interp1([wav(p,i) wav(q,i)],[p q],threshold*max(wav(:,i)));
            retr_corr(i)=(retr_bin-ref_bin)*binsize;
        end

    end
end