%% Suggested solution: Ultrasound signals from moving targets
% 
% <latex>
% % Original version JK, 03.2001
% % Updated by TH, 02.2010. Removed all references to Readecho and added
% % the inline latex code. Use the publish-command in Matlab to make a nice
% % latex-file, or preferably publish_to_latex(this_file_name)
% </latex>
% 

%% Part 1: Analytic expression for velocity
%
% <latex>
% Finding $z$ using the rule of cosines, solving for z, and using $R<<L$
% to approximate:
% \begin{align}
% L^2 &= z^2+R^2-2zR\cos(\theta)\\
% 0&=z^2-2R\cos(\theta)z + (R^2-L^2)\\
% z &= R \cos(\theta) \pm \sqrt{R^2(\cos(\theta)-1) + L^2}\\
% &\approx R\cos(\theta) \pm L\\
% \intertext{$\theta$ is found from the period $T$ and the time $t$:}
% \theta&=\frac{2\pi(t-t_0)}{T}
% \end{align}
% Finding the velocity of the piston by differentiation of $z$:
% \begin{align}
% z'(t)&=-2\pi R/T \sin(2\pi(t-t_0)/T)\\
% \intertext{The velocity of a point $r$ relative to the probe surface
%   is the opposite of $z'(t)$:} 
% r'(t)&=-z'(t)\\
% &=2\pi R/T \sin(2\pi(t-t_0)/T)
% \end{align}
% At $t_0$ the probe is closest to the phantom, giving the minimum value
% of $r$ and the maximum value of $z$.
% </latex>
%

%% Part 2: M-mode
%
% Start by loading the ultrasound data with settings:

load slowmotion;

%%
%The iq matrix is organized as samples x beams x frames. Find a beam near
%the center of the image, and extract all frames to make an M-mode. Use
%"squeeze" to remove singleton dimension.

iqmm=squeeze(iq(:,round(end/2),:));
 
%%
%Make time axis and depth axis:
timeaxis=[1:s.iq.FramesIQ]/s.Framerate_fps;
depthaxis=[1:s.iq.SamplesIQ]*s.iq.DepthIncrementIQ_m*100; %cm

%%
%Detect, calculate power, log compress, and show M-mode:

figure(1);clf;
gain=-40;dyn=40;
iqmmlog=imagelog(abs(iqmm).^2,gain,dyn);
image(timeaxis,depthaxis,iqmmlog);colormap(gray(64));
xlabel('Time [s]');
ylabel('Depth [cm]');
title('M-mode of IQ-beams');

%%
%Find three points: first peak, first trough and second peak using
%[x,y]=ginput(3);
x=[0.0708;0.528;0.979];y=[2.83;4.17;2.83];

%%
%Calculate R, T and t0:
R=(y(2)-y(1))/2;
T=x(3)-x(1);
t0=x(1);

%%
%Calculate the curve r(t) passing through the three points and plot on
%top of the M-mode
r0=(y(1)+y(2))/2; 
r=r0-R*cos(2*pi*(timeaxis-t0)/T);

hold on;
plot(timeaxis,r,'y-');

%%
%Show your Matlab magic and plot T, R og t0 over M-mode:
plot(x([1,3]),y([1,3]),'w+-');
h=text(mean(x([1 3])),y(1),['T= ', num2str(T,3),' s'],'color','w');
set(h,'verticalalignment','bottom','horizontalalignment','center');

plot([ min(xlim) x(1)],y([1,1]),'w+-');
h=text(mean([min(xlim) x(1)]),y(1),['t_0= ',num2str(t0,3),' s'],'color','w');
set(h,'verticalalignment','bottom','horizontalalignment','center');

plot(x([2 2]),y([2 1]),'w+-');
h=text(x(2),mean(y([1 2])),['2R = ',num2str(2*R*10,3),' mm'],'color','w');
set(h,'verticalalignment','middle','horizontalalignment','left');

disp(['Measured Values: R=' num2str(R) 'mm, T=' num2str(T) 's, t_0=' ...
      num2str(t0) 's. The analytical curve fits well with the movement as observed from the M-mode'])

%% Part 3: Time shift in RF data 
%
%Plot speed versus frame number

figure(2);clf;

%Velocity, r'(t):
rt=2*pi*R/T*sin(2*pi*(timeaxis-t0)/T);
plot(rt,'k-');
xlabel('Frameno');ylabel('cm/s');

%%
%Find a frame where the velocity is approximately one third of the maximum velocity
maxVel=max(rt);
frameno=min(find(rt>maxVel/3));
velocity=rt(frameno);
hold on;
plot(frameno*[1 1],ylim,'r--');
title(['Frameno: ',num2str(frameno)]);

%%
%Convert two subsequent beams in the IQ matrix to RF data with sampling frequency 200MHz
fs=200e6;Ts=1/fs;
demod=s.iq.fDemodIQ_Hz; %IQ demodulation freq.
frsIQ=s.iq.frsIQ_Hz; %IQ sampling freq
xrf=iq2rf(iqmm(:,frameno),demod,frsIQ,fs);
yrf=iq2rf(iqmm(:,frameno+1),demod,frsIQ,fs);

%%
%Time axis for the RF signal ("fast time") 
RFtimeaxis=[1:length(xrf)]*Ts*1e6; %us

%%
%Plot the two RF signals in the same axis, and zoom in
figure(3);clf;
subplot(2,1,1); hold on;
plot(RFtimeaxis,xrf,'r-');
plot(RFtimeaxis,yrf,'k--');
xlabel('Time [\mus]');ylabel('RF-value');title('RF-signals');
legend(['x_{RF} (Frame ',num2str(frameno),')'],...
       ['y_{RF} (Frame ',num2str(frameno+1),')']);

%%
%Show which part you are zooming into
plot([50 50],ylim,'k:');
plot([51 51],ylim,'k:');

%%
%Plot the zoomed area in a separate subplot
subplot(2,1,2);hold on;
plot(RFtimeaxis,xrf,'k-');
plot(RFtimeaxis,yrf,'k--');
xlabel('Time [\mus]');ylabel('RF-value');
set(gca,'ytick',0,'ygrid','on')
xlim([50 51]);

%%
%Indicate time shift between RF signals, found from ginput(2)
x=[50.530;50.605]; %us
rf_lag=diff(x)*1e3;%ns
%Fancy index-shuffling to plot a nice legend
x=x([1,1,1,2,2,2]);
y=[0 800 900];y=y([1,3,2,2,3,1]);

plot(x,y,'k-');
text(mean(x),max(y),['\Deltat=',num2str(rf_lag,3),' ns'],...
    'verticalalignment','bottom',...
    'horizontalalignment','center');

disp(['The time lag between the two beams was found to be ' num2str(rf_lag,3) ' ns.']);

% %% Task 4: Cross correlation
% %
% %Select 5000 samples along the beams, and calculate the cross correlation
% %with lags ranging from -100 to 100 samples
% 
% sumindex=[15000:20000];
% lag=[-100:100];
% ryx=zeros(size(lag));
% for n=1:length(lag),
%     ryx(n)=sum(yrf(sumindex).*xrf(sumindex+lag(n)));
% end;
% 
% %%
% %Normalize
% ryx=ryx/max(ryx);
% 
% %%
% %Time lag axis [ns]
% lagaxis=lag*Ts*1e9;
% 
% %%
% %Plot the cross correlation, and find the lag with maximum
% %correlation to determine the time shift between the two beams.
% figure(4);clf;hold on;
% plot(lagaxis,ryx,'k-');
% 
% title(['Cross-correlation. File: slowmotion.mat',...
%        '   Frameno:',num2str(frameno),'   Velocity: ',num2str(velocity,2),' cm/s']);
% ylabel('r_{yx}(l)'); xlabel(['Time lag [ns]']);
% ylim([-1 1.8]);
% 
% %Find max cross correlation
% [max_ryx,index]=max(ryx);
% maxlag=lagaxis(index);
% plot(maxlag+[0 0],ylim,'r');
% 
% %%
% %Compare with the expected lag and the lag found from inspecting the RF data: 
% expected_lag=2*velocity/100/(s.Framerate_fps*1540)*1e9;%ns
% plot(expected_lag*[1 1],ylim,'k--',rf_lag*[1 1],ylim,'b:');
% 
% legend('r_{yx}',...
%        ['Max. cross correlation: ',num2str(maxlag,3),' ns'],...
%        ['Expected lag: ',num2str(expected_lag,2),' ns'],... 
%        ['Lag found from RF-signal: ',num2str(rf_lag,3),' ns'])
% 
% %%
% %The plot of the cross correlation of the two beams peaks at lag 70ns. This
% %is approximately the same as found from plotting the RF signal, but more
% %than what was found from the analytical equations. Note that the peak at
% %-320ns is almost as strong as the one detected. (Do you understand why?)
% %
% %Now calculate the phase shift:
% %
% %<latex>
% %
% %\begin{align}
% %\frac{\Delta\theta}{2\pi}&=\frac{\Deltat}{T}=\Deltatf\\
% %&\Downarrow\\
% %\Deltat&=2\pi\Deltatf
% %\end{align}
% %
% % With center frequency $f_0=2.5MHz$ we get these phase shifts:
% %\begin{align}
% %\Deltat=58ns &\Rightarrow \Delta\theta = 0.91 rad\\
% %\Deltat=70ns &\Rightarrow \Delta\theta = 1.10 rad\\
% %\Deltat=75ns &\Rightarrow \Delta\theta = 1.18 rad\\
% %\end{align}
% %</latex>
% 


%% Task 5: The autocorrelation method
%
% Now we will use the autocorrelation method to estimate the phase shift
% in IQ data.

%%
%Make two matrixes with lag 1 in time
iq1=iqmm(:,frameno);
iq2=iqmm(:,frameno+1);

%%
%Calculate the autocorrelation for all depth samples:
Rxx=conj(iq1).*iq2;
figure(4), plot( depthaxis, angle( Rxx), 'k-' );
xlabel('Depth [cm]');
ylabel('Phase shift');
title('Phase shift vs depth');

% Average 20 samples near center of beam
depthindex=round(s.iq.SamplesIQ/2)+[-10:10]; 
avgRxx = angle(mean( Rxx(depthindex,:), 1) );
disp(['Phase shift in the IQ-data:', num2str(avgRxx,3),' rad']);

%%
% This means that the phase shift estimated by the autocorrelation method
% is approximately the same as the one obtained by reading the time
% shift in the RF signal.
% <latex> Phase shift ($\Delta\theta$) is the relation between time
% displacement and pulse period, $T = 1/f_0$. Displacement of one period
% corresponds to a phase shift of $2\pi$ radians. We therefore, get the
% following relationship between phase shift and velocity:
% \begin{align}
%  \frac{\Delta\theta}{2\pi} &= \frac{\Delta t}{T}\\
%  &\Downarrow\\
%  \Delta\theta &= 2\pi \frac{\Delta t}{T}\\
%  &=2\pi\Delta t f\\
%  &= \frac{4\pi\Delta r f }{c}\\
%  &= \frac{4\pi v f}{c FR}\\
%  &\Downarrow\\
%  v&=\frac{\Delta\theta c FR}{4\pi f}
% \end{align}
% <\latex>
% 

%%
%Phase shift and velocity of the whole sequence
%
%Find the real velocity from analytic calculation
v_real=2*pi*R/T*sin(2*pi*(timeaxis-t0)/T);

%%
%Make two matrices with lag 1 in time
frames=s.iq.FramesIQ;

iq1=iqmm(depthindex,1:frames-1);
iq2=iqmm(depthindex,2:frames);

%%
%Calculate the autocorrelation, and average in the depth direction
Rxx=mean(conj(iq1).*iq2);

%Estimate the velocity from the angle of the autocorrelation with lag 1
c=1540;
f0=s.iq.TxFreqIQ_Hz;
framerate=s.Framerate_fps;
v_est=c*framerate/(4*pi*f0)*angle(Rxx);

figure(5);clf;hold on;
plot(timeaxis,v_real,'k-');
plot(timeaxis(2:end),v_est*100,'k:');
xlabel('Time [s]');ylabel('Velocity [cm/s]');
title('Slowmotion (middle beam)');
legend('Calculated velocity','Measured velocity'); 

%%
%The plots of analytical and estimated velocities correspond quite
%well. There is a little bit of ring movement in the piston and the rotation
%in the motor is a little bit uneven with low rotational velocity. This can
%explain that the measured velocity is a little shifted compared to the
%calculated velocity.

%%Part 6: Aliasing
%

load fastmotion;

%%
%coordinates of two peaks and the trough inbetween
x=[0.0419;0.356;0.675];y=[2.75;4.12;2.78];
R=(y(2)-y(1))/2;
T=x(3)-x(1);
t0=x(1);%R in cm, T and t0 in seconds

%Time axis:
framerate=s.Framerate_fps;
frames=s.iq.FramesIQ;
timeaxis=[1:frames]/framerate;

%Find the real velocity from analytic calculation
v_real=2*pi*R/T*sin(2*pi*(timeaxis-t0)/T);

%Find a beam near the center
iqmm=squeeze(iq(:,round(end/2),:));

%%
%Choose approximately 20 samples near the highest amplitude:
%Make two matrixes with lag 1 in time
iq1=iqmm(depthindex,1:frames-1);
iq2=iqmm(depthindex,2:frames);

%%
%Calculate the autocorrelation with averaging of the depth samples:
Rxx=mean(conj(iq1).*iq2);
disp(['Phase shift in the IQ-data:', num2str(angle(Rxx),3),' rad']);

%Estimate the velocity from the angle of the autocorrelation with lag 1
c=1540;
f0=s.iq.TxFreqIQ_Hz;
framerate=s.Framerate_fps;
v_est=c*framerate/(4*pi*f0)*angle(Rxx);

figure(6);clf;hold on;
plot(timeaxis,v_real,'k-');
plot(timeaxis(2:end),v_est*100,'k:');
xlabel('Time [s]');ylabel('Velocity [cm/s]');
title('Fastmotion (middle beam)');
legend('Calculated velocity','Measured velocity'); 

%%
%When the velocity is high, the phase shift between subsequent samples
%become greater than $\pi$ radians. Because the phase angle to the complex
%number is defined between $-\pi$ and $\pi$, the phase angle over $\pi$
%radians "folds over" to the negative phase angle.

% %%Part 7: Angle between velocity and beam
% %
% %Going back to the slowmotion clip
% load slowmotion;
% 
% %%
% %Find three points: first peak, first trough and second peak using
% %[x,y]=ginput(3);
% x=[0.0708;0.528;0.979];y=[2.83;4.17;2.83];
% 
% %%
% %Calculate R, T and t0:
% R=(y(2)-y(1))/2;
% T=x(3)-x(1);
% t0=x(1);
% 
% %Time axis:
% framerate=s.Framerate_fps;
% frames=s.iq.FramesIQ;
% timeaxis=[1:frames]/framerate;
% 
% %%
% %Finding the angle with the beam:
% beamno=1;
% startangle=s.iq.StartAngleIQ_rad;
% angleincrement=s.iq.AngleIncrementIQ_rad;
% phi=startangle+((beamno-1)*angleincrement);
% 
% %%
% %Analytic expression of velocity, with angle correction
% v_real=2*pi*R/T*sin(2*pi*(timeaxis-t0)/T);
% v_anglecorrected=v_real*cos(phi);
% 
% %Get an M-mode from the first beam:
% iqmm=squeeze(iq(:,1,:));
% 
% %%
% %Choose approximately 20 samples in the middle of the beam:
% samples=s.iq.SamplesIQ;
% depthindex=round(samples/2)+[-10:10]; 
% 
% %%
% %Make two matrixes with lag 1 in time
% iq1=iqmm(depthindex,1:frames-1);
% iq2=iqmm(depthindex,2:frames);
% 
% %%
% %Calculate the autocorrelation with averaging of the depth samples:
% Rxx=mean(conj(iq1).*iq2);
% 
% %Velg ut ca. 20 dybdesamples midt p? str?len
% samples=s.iq.SamplesIQ;
% depthindex=round(samples/2)+[-10:10]; 
% 
% %Calculate the velocity:
% c=1540;
% f0=s.iq.TxFreqIQ_Hz;
% framerate=s.Framerate_fps;
% v_est=c*framerate/(4*pi*f0)*angle(Rxx);
% 
% %Plot the calculated and estimated velocity in the same figure
% figure(7);clf;hold on;
% plot(timeaxis,v_anglecorrected,'k-');
% plot(timeaxis(2:end),v_est*100,'k:');%one sample is "lost" in the
%                                      %estimation due to the lag
% xlabel('Time [s]'); ylabel('Velocity [cm/s]');
% title(['Slowmotion. Beam: ',num2str(beamno)]);
% legend('Computed angle corrected velocity','Measured velocity');
% ylim([-5 6]);
% 
% %%
% %The agreement between the measured velocity and the angle-corrected
% %velocity is quite good.
