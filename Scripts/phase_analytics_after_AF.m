%                 %% debug
%                 figure;
%                 subplot(1,2,1)
%                 [~,I_wav]=max(abs(obj.wav_select(floor(length(obj.wav_select(:,1))/2),:)));
%                 phase = angle(obj.wav_select(:,I_wav));
%                 phase = mod(phase+pi,2*pi);
%                 plot(obj.t_pulse,phase);hold on
%                 plot(obj.t_select_corr(1:end-1),-50*diff(Rtrk_select));hold on
%                 plot(obj.t_select(1:end-1),diff(Rtrk_select));hold on
%                 
%                 subplot(1,2,2)
%                 plot(mod(diff(phase),2*pi)/(2*pi),'ro')
%                 %%
%                 figure;
%                 scatter(diff(Rtrk),mod(diff(phase(1:end)),2*pi),'.'); hold on
%                 grid()
%                 xlabel('\Delta Rtrk in m')
%                 ylabel('\Delta phase')
%                 
%                 r = -(1:50)*0.00295
%                 
%                 plot(r,mod(-P(1)*(r+0.0011),2*pi),'o')
%                 
%                 %% fit line behaviour
%                 x = [0.0813 0.09315 0.09615 0.0991]
%                 y = 2*pi*[0.36-1 0.43 0.7 0.97]
%                 figure; plot(x,y,'ro'); hold on
%                 P = polyfit(x,y,1);
%                 plot(x,P(1)*x+P(2))
%                 
%                 %% try to correct the phase in place:
%                 phase_corr = mod(phase + P(1)*Rtrk,2*pi)
%                 figure;
%                 plot(phase);hold on
%                 plot(phase_corr);hold on
%                 