clear all
close all
clc

%% definig some variables from Ray et al. (2015) as in SAMOSAfun
DDAcf = DDA_ConfigFile('S6A','SAR');

DDAcf.PRI = median(CS1b.MEA.PRI);
DDAcf.fp = 1/DDAcf.PRI;
h = median(CS1b.MEA.tracker_range);
V = median(CS1b.GEO.V.V);
dHdt = median(CS1b.GEO.H_rate);
vr = dHdt;
vt = sqrt(V.^2 - vr.^2);
T = 3                       % coherent integration time

os_ZP     = DDAcf.os_ZP;                             %Zero-Padding Oversampling Factor
B_os      = DDAcf.B * os_ZP;
%RefBin    = (DDAcf.RefBin-1)*os_ZP + 1;                    %Reference bin
alpha     = 1%+h/Re;                                        %Geometric parameter
t_res     = 1/B_os;                                     %Time resolution [s]
t_s       = t_res;                                   %Waveform sampling interval [s] (~= t_res because of zero-padding)
Lx        = (DDAcf.c*h*DDAcf.fp)/(2*vt*DDAcf.fc*DDAcf.Nb); %Along-track resolution (Ray et al. (2015), Eq. 14)
Lx_FF        = (DDAcf.c*h)/(2*vt*DDAcf.fc*T); %Alternative for fully focused SAR, since the PRI*Number of pulses per burst was assumed as integration time, now we choose T
Ly        = sqrt((DDAcf.c*h)/(alpha*B_os));             %Ray et al. (2015), Eq. 21
Lz        = DDAcf.c/(2*B_os);                           %Range resolution (Ray et al. (2015), Eq. 21)
gamma_x   = 8*log(2)/DDAcf.theta_x^2;                      %Along-track antenna gain parameter (Dinardo et al. (2018), Eq. 27)
gamma_y   = 8*log(2)/DDAcf.theta_y^2;                      %Cross-track antenna gain parameter (Dinardo et al. (2018), Eq. 27)
L_gamma   = alpha*h/(2*gamma_y);                           %(Dinardo et al. (2018), Eq. 27)

%% plotting and looking at Ckl2

osf = 0.025

l = osf*(-64/(2*osf) : 64/(2*osf));
x = 0;
y = 0;
z = 0;

u = x/Lx - l;
kl = Lx^2/Ly^2*(u.^2 + 2*l.*u) + y^2/Ly^2 - z/Lz;

k = osf*(-DDAcf.os_ZP*DDAcf.Np/(2*osf) : DDAcf.os_ZP*DDAcf.Np/(2*osf))';

Ckl2 = DDAcf.Nb^2 .* sinc(u).^2 .* sinc(k - kl).^2;

figure;
subplot(2,1,2)
imagesc(log10(Ckl2));
caxis([-5,1])
subplot(2,1,1)
plot(sum(Ckl2,1))

% the remaining parabola is reasonable, because Ckl describes where we see
% a point target at x,y,z in the doppler beams from one burst above it.
% Now, we need to adjust our point of view, hence, place the point target
% at along track (x) positions that are always a distance d from the focus
% point at l.

%% plotting and looking at Ckl2 with all beams towards a single target

osf = 0.025;
l = osf*(-64/(2*osf) : 64/(2*osf)); % how many doppler beams are used will determine the effective integration time and the range widening!
k = osf*(-(DDAcf.os_ZP*DDAcf.Np-1)/(2*osf) : DDAcf.os_ZP*DDAcf.Np/(2*osf))';

d = 1000; % distance focus point to target in m
x = d + Lx*l;
y = 0;
z = 0;

u = x/Lx - l;
kl = Lx^2/Ly^2*(u.^2 + 2*l.*u) + y^2/Ly^2 - z/Lz;

k = osf*(-DDAcf.os_ZP*DDAcf.Np/(2*osf) : DDAcf.os_ZP*DDAcf.Np/(2*osf))';

Ckl2 = DDAcf.Nb^2 .* sinc(u).^2 .* sinc(k - kl).^2;

figure;
subplot(2,1,2)
imagesc(log10(Ckl2));
caxis([-5,1])
subplot(2,1,1)
plot(sum(Ckl2,1))


%% plotting and looking at Ckl2 in an IRF point of view (antenna gain pattern missing)

osf = 0.5;
l = osf*(-64/(2*osf) : 64/(2*osf)); % how many doppler beams are used will determine the effective integration time and the range widening!
%l = osf*(-128/(2*osf) : 128/(2*osf));
k = osf*(-(DDAcf.os_ZP*DDAcf.Np-1)/(2*osf) : DDAcf.os_ZP*DDAcf.Np/(2*osf))';

Nd = 2000;
ds = 2*((1:Nd)-Nd/2);
IRF = zeros(size(k,1),Nd);

for i = 1:Nd
    %x = 0;
    d = ds(i);
    x = d + Lx*l;
    y = 0;
    z = 0;

    u = x/Lx - l;
    kl = Lx^2/Ly^2*(u.^2 + 2*l.*u) + y^2/Ly^2 - z/Lz;

    Ckl2 = DDAcf.Nb^2 .* sinc(u).^2 .* sinc(k - kl).^2;

%     figure;
%     subplot(2,1,2)
%     imagesc(log10(Ckl2));
%     caxis([-5,1])
%     subplot(2,1,1)
%     plot(sum(Ckl2,1))
    IRF(:,i) = sum(Ckl2,2);
end

figure;
subplot(2,1,2)
imagesc(ds,k,log10(IRF./max(IRF(:))));
caxis([-7,0])
colormap('pink')
subplot(2,1,1)
plot(sum(IRF,1))

%% now, try to adjust this result so it yields something like the FF-SAR IRF with the same integration time by adjusting Lx

osf = 0.5;
l = osf*(-64/(2*osf) : 64/(2*osf)); % how many doppler beams are used will determine the effective integration time and the range widening!
%l = osf*(-128/(2*osf) : 128/(2*osf));
k = osf*(-(DDAcf.os_ZP*DDAcf.Np-1)/(2*osf) : DDAcf.os_ZP*DDAcf.Np/(2*osf))';

Nd = 2000;
ds = 2*((1:Nd)-Nd/2);
IRF = zeros(size(k,1),Nd);

for i = 1:Nd
    %x = 0;
    d = ds(i);
    x = d + Lx*l; % hier soll noch immer Lx stehen, denn der gedachte focal point (fuer den Doppler beam) liegt ja bei Lx*l, sodass wir immer um d entfernt neben den focal point schauen wollen
    y = 0;
    z = 0;

    %u = x/Lx_FF - l;
    u = x/Lx - l;
    kl = Lx^2/Ly^2*(u.^2 + 2*l.*u) + y^2/Ly^2 - z/Lz;
    
    %u_FF = x/Lx_FF - l;
    Ckl2 = DDAcf.Nb^2 .* sinc(d/Lx_FF).^2 .* sinc(k - kl).^2; % here is the only change to make it visually fit. Doesn't mean anything...

%     figure;
%     subplot(2,1,2)
%     imagesc(log10(Ckl2));
%     %caxis([-5,1])
%     subplot(2,1,1)
%     plot(sum(Ckl2,1))
    IRF(:,i) = sum(Ckl2,2);
end

figure;
subplot(2,1,2)
imagesc(ds,k,log10(IRF./max(IRF(:))));
caxis([-7,0])
colormap('pink')
subplot(2,1,1)
plot(sum(IRF,1))

