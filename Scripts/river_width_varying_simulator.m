
%% first, assuming the case of symmetric widening of the river extent:

w = [50, 100, 300, 500, 1000, 2000]
%dx = -100:1:100
Hs = 1348000
RE = 6371000

fig = figure('units','inch','position',[0,0,5,3.5]);
subplot(1,2,1)
leg = {}
sym = ({'.','o','x','+','*','s','d','v','^','<','>','p','h'})
for i = 1:numel(w)
    dx = -w(i)/3:10:w(i)/3
    dH = (1/(2*Hs)+1/(2*RE))*(dx*w(i)+dx.^2);
    dw = 2*dx;
    ax = plot(dw(1:5:end),100*dH(1:5:end),[sym{i} '-']);hold on
    leg{i} = ['w_c^*=' num2str(w(i)) ' m']
end

legend(leg,'NumColumns',2,'Location','northoutside')
xlabel('                                                       error in cross-track projected river width \Deltaw (m)')
ylabel({'error in WSE \DeltaH_r (cm)'})
grid on
set(gcf,'color','w');
xlim([-600,600])

subplot(1,2,2)
for i = 1:numel(w)
    dx = -w(i)/3:10:w(i)/3
    dH = (1/(2*Hs)+1/(2*RE))*(dx*w(i)+dx.^2);
    dw = 2*dx;
    ax = plot(dw,100*dH,[sym{i} '-']);hold on
    leg{i} = ['w_c^*=' num2str(w(i)) ' m']
end
grid on
set(gcf,'color','w');
xlim([-99,99])
ylim([-4.9,4.9])
legend(leg,'NumColumns',2,'Location','northoutside')

exportgraphics(gcf,'systematic_bias_river_width.pdf') 


%% second, assuming the case of random errors in river cross-track position: checking the formula

w = 100;%[50, 100, 300, 500, 1000, 2000]
alpha = 0;
wc = w/cos(alpha)
x_r = 7000;
%dx = -100:1:100
Hs = 1348000
RE = 6371000
K = (1/(2*Hs)+1/(2*RE))
sigma_x = 10

% assume that dx_ra and dx_rb are gaussian random numbers 
dxa = normrnd(0,sigma_x,100000000,1);
dxb = normrnd(0,sigma_x,100000000,1);

% for random uncorrelated case:
DH = K/2*( 2*(x_r-wc/2).*dxa + 2*(x_r+wc/2).*dxb + dxa.^2 + dxb.^2); % we can possibly neglect the squared terms

% for the symmetrically retreating/expanding case
% dxa = -20;%
% dxb = -20;
% DH = K/2*( 2*(x_r-wc/2).*dxa - 2*(x_r+wc/2).*dxa + dxa.^2 + dxa.^2);

mean(DH)
K*sigma_x^2
std(DH)
sqrt(2*K^2*(x_r^2+wc^2/4)*sigma_x^2)


%% plot the results 
w = [50, 100, 300, 500, 1000, 2000]
x_r = [0 1000 2000 3000 4000 5000 6000 7000]

fig = figure('units','inch','position',[0,0,5,3.5]);
subplot(1,2,1)
leg = {}
sym = ({'.','o','x','+','*','s','d','v','^','<','>','p','h'})
for i = 1:numel(w)
    sigmaDH = K*sigma_x*sqrt(2*x_r.^2 + w(i).^2/2)
    
    ax = plot(x_r/1000,100*sigmaDH,[sym{i} '-']);hold on
    leg{i} = ['w_c^*=' num2str(w(i)) ' m']
end

legend(leg,'NumColumns',2,'Location','northoutside')
xlabel('cross-track distance x_r (km)')
ylabel({'random error of WSE \sigma_{\DeltaH} (cm)'})
% grid on
% set(gcf,'color','w');
% xlim([-600,600])

% subplot(1,2,2)
% for i = 1:numel(w)
%     dx = -w(i)/3:10:w(i)/3
%     dH = (1/(2*Hs)+1/(2*RE))*(dx*w(i)+dx.^2);
%     dw = 2*dx;
%     ax = plot(dw,100*dH,[sym{i} '-']);hold on
%     leg{i} = ['w_c^*=' num2str(w(i)) ' m']
% end
grid on
set(gcf,'color','w');
%xlim([-99,99])
%ylim([-4.9,4.9])
legend(leg,'NumColumns',2,'Location','northoutside')

exportgraphics(gcf,'random_error_cross_track_dist.pdf') 
