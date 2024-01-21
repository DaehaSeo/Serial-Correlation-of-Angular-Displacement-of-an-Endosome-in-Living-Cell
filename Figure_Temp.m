%% Description:

%       Detailed pipeline of data processing to identify timing and degree of Au NR 3D rotation.
%       Figure showing specific region (Intensity, Azimuthal, Polar angle, time, Trajectory) for the manuscript figure 2 

% Lead Contact: Daeha Seo <livewire@dgist.ac.kr>
% Credit: Siwoo Jin <jhg9500@dgist.ac.kr>
% Reference: Siwoo Jin et al., Adv. Sci., 2024

%% Crop for manuscript figure 2
close all;

Intensity_domain = 2788:4065; %Manuscript figure 2, time domain for Intensity 
STFT_domain = 2788:4021; %Manuscript figure 2, angle domain for polar angle and azimthal angle
%figure(); 
%subplot(3,1,1); 

figure('pos',[0, 0, 1500, 900]);

%% Intensity time trace
subplot(5,1,1);
hold on;

title('Intensity  (Range Restrict)','FontSize',12,'FontWeight','bold');
plot(Time(Intensity_domain),Intensity(Intensity_domain), 'color','#FF9900','linewidth',1.5);
xlabel('Time(s)'); 
ylabel('Scattering intensity'); 
xlim([min(Time(Intensity_domain)),max(Time(Intensity_domain))]);
ylim([0, 30000]);
ax=gca;ax.LineWidth=1.2;

%% Angle for Phi and state time trace
subplot(5,1,2);
hold on;

title('Phi State (Range Restrict)','FontSize',12,'FontWeight','bold');
plot(Time(STFT_domain),Final_Phi(STFT_domain)*180/pi, 'color','#cccc00','linewidth',1.5); % change of radian to degree value.
plot(Time(STFT_domain),Final_Phi_State(STFT_domain)*180/pi, 'color','#99cc00','linewidth',3); % change of radian to degree value.
xlabel('Time(s)'); 
ylabel('Azimuthal angle (degree)');
xlim([min(Time(STFT_domain))-0.246682,max(Time(STFT_domain))+0.246682]); % Adding and subtracting the 0.24 value were based on calculations from the midpoint of cycle
ylim([-pi/2*180/pi, pi/2*180/pi]*1.3);
set(gca,'YTick',[-90:45:90]);
ax=gca;ax.LineWidth=1.2;


%% Angle for Theta and state time trace
subplot(5,1,3);
hold on;

title('Theta State (Range Restrict)','FontSize',12,'FontWeight','bold');
plot(Time(STFT_domain),Final_Theta(STFT_domain)*180/pi,'color','#99cc66','linewidth',1.5);
plot(Time(STFT_domain),Final_Theta_State(STFT_domain)*180/pi, 'color','#33cccc','linewidth',3);
xlabel('Time(s)'); 
ylabel('Polar angle (degree)');
xlim([min(Time(STFT_domain))-0.246682,max(Time(STFT_domain))+0.246682]); % Adding and subtracting the 0.24 value were based on calculations from the midpoint of cycle
ylim([0, pi/2*180/pi]*1.3);
set(gca,'YTick',[0:30:90]);
ax=gca;ax.LineWidth=1.2;


%% Diffusion coefficient time trace
subplot(5,1,4);
hold on;
title('Diffusion Coefficient (Range restrict)','FontSize',16,'FontWeight','bold');

%marking of change point for 3D rotation
for i = 1:length(icp)
    plot(STFT_Time([icp(i), icp(i)]-21),[min(D), max(D)],'color','k','LineStyle',':');
end

%marking of diffusion coefficient in the region between changes in angle using D cut off value,
%and distinguish whether it is in a 'Run' or 'Pause' state.
 for i = 1:(length(icp)-1)
     if D_PauseZeroMoveOne(i) == 1
         x = STFT_Time([Newipt(i), Newipt(i+1), Newipt(i+1), Newipt(i)]);
         y = [min(D), min(D), max(D), max(D)];
         p = patch(x,y,'k');
         set(p,'FaceAlpha',0.5,'EdgeAlpha',0);
     end
 end

plot(Time(STFT_domain),D(STFT_domain),'color','b','linewidth',2);
xlim(([min(Time(STFT_domain))-0.246682,max(Time(STFT_domain))+0.246682]));
ylim([min(D), max(D)]);

xlabel('Time(s)'); 
ylabel('Diffusion coefficient(\mum^2/s)');
set(gca,'FontSize',8,'FontWeight','bold');

%% Trajectory
subplot(5,1,5);
hold on;
title('Trajectory','FontSize',12,'FontWeight','bold');
plot(Track(Intensity_domain,1)*LengthPerPixel',Track(Intensity_domain,2)*LengthPerPixel','k-');
plot([10.5 10.5 11.5],[11 11 11],'k-','linewidth',3);
text(11, 11,'1\mum','horiz','center','vert','top'); 

pbaspect([1 1 1]);
xlim([6, 12]);
ylim([10, 16]);
axis off 




