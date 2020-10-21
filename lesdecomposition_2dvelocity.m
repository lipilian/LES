clear all
%Input=load('out.txt');
Input = import(3); 
X=Input(:,1);               % Extract X/Y-coordinates
Y=Input(:,2);
U=Input(:,3);               % Extract U/V-velocities
V=Input(:,4);
Vel=U+V*i;

Kernel=repmat([-3:3],7,1).^2 ...    % Neighborhood radius squared
      +repmat([-3:3]',1,7).^2;
Kernel=exp(-4*Kernel/200);           % Gaussian, Width^2=16
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

% Xlo1=1.163;
% Xlo2=1.4611;
% Ylo1=0.1923;
% Ylo2=0.1827;
% 
% slope=(Ylo2-Ylo1)/(Xlo2-Xlo1);
% 
% for i=1:length(U)
%     if Y(i)<Ylo1+0.03+slope*(X(i)-Xlo1);
%         if mag(i)>0.08*maxvalue
%             U1(i)=0;
%             V1(i)=0;
%         end
%     end
%     if mag(i)<0.01*maxvalue
%         U1(i)=0;
%         V1(i)=0;
%     end
%     if X(i)>2.34 && Y(i)>-0.666
%         U1(i)=0.2*U1(i);
%         V1(i)=0.2*V1(i);
%     end
%     if X(i)<0.25
%         U1(i)=0.2*U1(i);
%         V1(i)=0.2*V1(i);
%     end
% end
%     

Output.name='LES decomposition';
Output.type='vectors';
Output.dataset.X=X;
Output.dataset.Y=Y;
Output.dataset.U1=U1;
Output.dataset.V1=V1;

U1_old = U1; V1_old = V1;
%%%%%%%% try filter using magnitude

%%%if large magnitude near by do not filter center
% k = [];
% U1_new = zeros(size(mag));
% V1_new = zeros(size(mag));
% for i = 1:(size(mag))
%      if mag(i,1)> 0.0004
%         %k = [k i+2 i+1 i i-1 i-2];
%         k = [k i i-1 i-2];
%         for j = k(1,:)
%             U1_new(j,1)= U1_old(j,1);
%             V1_new(j,1)= V1_old(j,1); 
%         end
%             
%      end
% end


%%%%%% filter by U, V magnitude seperately
% h = []; l = [];
% for i = 1:size(mag)
%     if U1_old(i,1) >0.0002
%         h = [h i i-1];
%         for j = h(1,:)
%             U1_new(j,1) = U1_old(j,1);
%         end
%     end
%     if V1_old(i,1) >0.0002
%         l = [l i i-1];
%         for j = l(1,:)
%             V1_new(j,1) = V1_old(j,1);
%         end
%     end
% end  

% direct filter those have magnitude smaller than threshold
% for i = 1:size(mag)
%      if mag(i,1)< 0.0003
%          U1(i,1) = 0;
%          V1(i,1) = 0;
%      end
% end
 
figure();
%quiver(X,Y,U1_new,V1_new,3,'k');
quiver(X,Y,U1,V1,3,'k');
hold on;
 

Vvor2=[];
Uvor2=[];
Xvornew=[];
Yvornew=[];
%K = [X Y U V U1_new V1_new];
K = [X Y U V U1 V1];
save 45hz_07hz_3.txt K -ascii





% Xvor2=linspace(1.763,1.85,30);
% Yvor2=linspace(0.37,0.37,30);

% Xvor2=linspace(1.635,1.7,30);
% Yvor2=linspace(0.34,0.34,30);
% for i=1:length(Yvor2)
%     ll=abs((Y-Yvor2(i)*ones(length(U),1)).^2+(X-Xvor2(i)*ones(length(U),1)).^2);
%     X00=find(ll==min(ll));
%     Vvor2=[Vvor2 V1(X00)];
%     Uvor2=[Uvor2 U1(X00)];
%     Xvornew=[Xvornew X(X00)];
%     Yvornew=[Yvornew X(X00)];
% end
% 
% Xvornew=Xvornew-ones(1,length(Xvornew))*Xvornew(1);
% 
% figure;
% plot(Xvornew,Vvor2/0.488,'k','LineWidth',1);
% 
% hold on;
% 
% set(gca,'box','on','FontSize',25);
% set(gcf, 'PaperPosition', [0 0 9 9]); 
% set(gcf, 'PaperSize', [9 9]);
% set(gca, 'Position', [0.18 0.18 0.75 0.7]);
% xlabel('$r/d$','Interpreter','LaTex','Fontsize',35);
% ylabel('$v"/V_c$','Interpreter','LaTex','Fontsize',35);
% 
% Vvor2=[];
% Uvor2=[];
% Xvornew=[];
% Yvornew=[];
% % Xvor2=linspace(1.763,1.85,30);
% % Yvor2=linspace(0.37,0.37,30);
% 
% Xvor2=linspace(1.635,1.5,30);
% Yvor2=linspace(0.34,0.34,30);
% for i=1:length(Yvor2)
%     ll=abs((Y-Yvor2(i)*ones(length(U),1)).^2+(X-Xvor2(i)*ones(length(U),1)).^2);
%     X00=find(ll==min(ll));
%     Vvor2=[Vvor2 V1(X00)];
%     Uvor2=[Uvor2 U1(X00)];
%     Xvornew=[Xvornew X(X00)];
%     Yvornew=[Yvornew X(X00)];
% end
% Xvornew=Xvornew-ones(1,length(Xvornew))*Xvornew(1);
% 
% plot(Xvornew,Vvor2/0.488,'k','LineWidth',1);
% 
% hold off;
% 
% Vvor2=[];
% Uvor2=[];
% Xvornew=[];
% Yvornew=[];
% % Xvor2=linspace(1.763,1.85,30);
% % Yvor2=linspace(0.37,0.37,30);
% 
% Xvor2=linspace(1.635,1.635,30);
% Yvor2=linspace(0.34,0.515,30);
% for i=1:length(Yvor2)
%     ll=abs((Y-Yvor2(i)*ones(length(U),1)).^2+(X-Xvor2(i)*ones(length(U),1)).^2);
%     X00=find(ll==min(ll));
%     Vvor2=[Vvor2 V1(X00)];
%     Uvor2=[Uvor2 U1(X00)];
%     Xvornew=[Xvornew X(X00)];
%     Yvornew=[Yvornew Y(X00)];
% end
% Yvornew=Yvornew-ones(1,length(Yvornew))*Yvornew(1);
% 
% figure;
% plot(Yvornew,Uvor2/0.488,'k','LineWidth',1);
% 
% hold on;
% 
% set(gca,'box','on','FontSize',25);
% set(gcf, 'PaperPosition', [0 0 9 9]); 
% set(gcf, 'PaperSize', [9 9]);
% set(gca, 'Position', [0.18 0.18 0.75 0.7]);
% xlabel('$r/d$','Interpreter','LaTex','Fontsize',35);
% ylabel('$v"/V_c$','Interpreter','LaTex','Fontsize',35);
% 
% Vvor2=[];
% Uvor2=[];
% Xvornew=[];
% Yvornew=[];
% % Xvor2=linspace(1.763,1.85,30);
% % Yvor2=linspace(0.37,0.37,30);
% 
% Xvor2=linspace(1.635,1.635,30);
% Yvor2=linspace(0.34,0.225,30);
% for i=1:length(Yvor2)
%     ll=abs((Y-Yvor2(i)*ones(length(U),1)).^2+(X-Xvor2(i)*ones(length(U),1)).^2);
%     X00=find(ll==min(ll));
%     Vvor2=[Vvor2 V1(X00)];
%     Uvor2=[Uvor2 U1(X00)];
%     Xvornew=[Xvornew X(X00)];
%     Yvornew=[Yvornew Y(X00)];
% end
% Yvornew=Yvornew-ones(1,length(Yvornew))*Yvornew(1);
% 
% plot(Yvornew,Uvor2/0.488,'k','LineWidth',1);
% 
% % Ustream=[];
% % Vstream=[];
% 
% % maxy=max(Y);
% % miny=min(Y);
% % 
% % for i=1:255
% %     Ustream=[Ustream;(U1((i-1)*255+1:(i-1)*255+255))'];
% %     Vstream=[Vstream;(V1((i-1)*255+1:(i-1)*255+255))'];
% % end
% 
% % for jj=0:0.02:2.5
% %     for ii=0:0.02:2.5
% % streamline(X(1:255),linspace(maxy,miny,255),Ustream,Vstream,jj,ii);
% % hold on;
% %     end
% % end
% 
% % Clean up workspace
% clear C Kernel

%% function for import data 
function inst = import(f_n) %f_n as image number
    for i= f_n
        if i<10
            temp=importdata(['45hz_07hz00000' num2str(i) '.T000.D000.P000.H000.L.vec']);
        end
        if i>=10&&i<100
            temp=importdata(['45hz_07hz0000' num2str(i) '.T000.D000.P000.H000.L.vec']);
        end
        if i>=100&&i<1000
            temp=importdata(['45hz_07hz000' num2str(i) '.T000.D000.P000.H000.L.vec']);
        end
        if i>=1000&&i<10000
            temp=importdata(['8hznew00' num2str(i) '.T000.D000.P000.H000.L.vec']);
        end   
    end
    inst=temp.data; 
end
