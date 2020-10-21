clear all
%Input=load('out.txt');
Input = import(3); 
X=Input(:,1);               % Extract X/Y-coordinates
Y=Input(:,2);
U=Input(:,3);               % Extract U/V-velocities
V=Input(:,4);
Vel=U+V*i;

Kernel=repmat([-20:20],41,1).^2 ...    % Neighborhood radius squared
      +repmat([-20:20]',1,41).^2;
Kernel=exp(-4*Kernel/1000);           % Gaussian, Width^2=16
Kernel=Kernel/sum(Kernel(:));       % Normalise kernel
C=filter2(Kernel,ones(size(Vel)));  % Map edge effects
Vel=filter2(Kernel,Vel)./C;         % Apply kernel & compensate edges
Vel=U+V*i-Vel;                      % Subtract mean-flow from original
Vel=filter2([1,2,1;...              % Mild Gaussian smoothing...
             2,4,2;...
             1,2,1]/200,Vel);
Vel(:,1)=Vel(:,1)*4/3;              % Compensate edge effects
Vel(1,:)=Vel(1,:)*4/3;
Vel(:,end)=Vel(:,end)*4/3;
Vel(end,:)=Vel(end,:)*4/3;

U1=real(Vel);                        % Extract new U-component
V1=imag(Vel);                        % Extract new V-component

mag=(U1.^2+V1.^2).^0.5;
maxvalue=max(mag);


 
figure();
%quiver(X,Y,U1_new,V1_new,3,'k');
quiver(X,Y,U1,V1,3,'k');
% hold on;
 

% Vvor2=[];
% Uvor2=[];
% Xvornew=[];
% Yvornew=[];
% %K = [X Y U V U1_new V1_new];
% K = [X Y U V U1 V1];
% save 45hz_07hz_3.txt K -ascii





