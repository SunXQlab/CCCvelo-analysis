
% PDE_SDE model for simulating single cell gene expressions involved in
% cell-cell communications

base_save_path = 'E:/CCCvelo/simulationData/data/fan_add_noise_20250807/';
%% creat save path
sim = 1;
save_path = fullfile(base_save_path, ['slide_' num2str(sim)]);
status = mkdir(save_path);
%% cellular scale initilization
N=100; % lattice size

SC=rand(N,N); % spatial cells
TC=zeros(N,N); % latent time

for x=1:N
    for y=1:N
       if (x-10)^2+(y-10)^2>=80^2
           SC(x,y)=0;
       end
       
       if x<10 || y<10
           SC(x,y)=0;
       end
    end
end
SC(SC>0.6)=1;SC(SC<0.6)=0;

for x=1:N
    for y=1:N
       TC(x,y)=sqrt((x-10)^2+(y-10)^2)/80;
    end
end
TC(SC==0)=0;
%% visualization
tic
figure,
for x=1:N
    for y=1:N
        if  SC(x,y)==1        
           plot(x,y,'c o','MarkerSize',3,'MarkerFaceColor',[0.5 TC(x,y) TC(x,y)],'LineWidth',1,'MarkerEdgeColor',[0.5 TC(x,y) TC(x,y) ]);hold on;
        end
    end
end
title('Cell Distribution');
box on;
axis([0 100 0 100]);

% ColorBar
colormap([0.5*ones(1,11); 0:0.1:1; 0:0.1:1]'); 
c = colorbar;
ylabel(c, 'Latent time');
saveas(gcf, fullfile(save_path, 'spatial_distribution.png'));
%% molecular scale initilization

% spatially similar expression
R1=(rand(N,N)*0.1+TC).*SC; % Receptor gene expression was assumed to be fixed but not evolved.
R2=(rand(N,N)*0.1+1-TC).*SC; 

tf1=(rand(N,N)*0.1+0.5).*SC; 
tf2=(rand(N,N)*0.1+0.5).*SC; 
tf3=(rand(N,N)*0.1+0.5).*SC; 


%% interaction parameters in the network
global beta1 beta2 beta3 mu11 mu12 mu13 mu14 mu21 mu22 mu23 mu24 mu31 mu32 mu33 mu34 
global ktt11 ktt12 ktt13 ktt14 ktt21 ktt22 ktt23 ktt24 ktt31 ktt32 ktt33 ktt34
%LR links
b11=1;b21=0;b31=1;b41=0;b51=1;  % Li-R1
b12=0;b22=0;b32=1;b42=1;b52=0;  % Li-R2
%Maximal reaction rate
alpha11=1*rand;alpha21=0; % Ri-TF1
alpha12=1*rand;alpha22=1*rand; % Ri-TF2
alpha13=0;alpha23=1*rand; % Ri-TF3
%MM coefficient 
krt11=rand;krt21=1; % Ri-TF1
krt12=rand;krt22=rand; % Ri-TF2
krt13=1;krt23=rand; % Ri-TF3

beta1=rand;beta2=rand;beta3=rand; % degradation rate of TFs

paraM1=zeros(5,3); % Ri-TFi parameters
paraM1(1,1)=alpha11;paraM1(2,1)=alpha21; % Maximal reaction rate, i.e., v
paraM1(1,2)=alpha12;paraM1(2,2)=alpha22;
paraM1(1,3)=alpha13;paraM1(2,3)=alpha23;
paraM1(3,1)=krt11;paraM1(4,1)=krt21; % MM coefficient, i.e., k
paraM1(3,2)=krt12;paraM1(4,2)=krt22;
paraM1(3,3)=krt13;paraM1(4,3)=krt23;
paraM1(5,1)=beta1;paraM1(5,2)=beta2;paraM1(5,3)=beta3; % degradation rate of TFs
 
%Maximal reaction rate
mu11=0.5+0.1*randn;mu12=1*0.5 + 0.1*randn;mu13=0;mu14=0; % TF1-TGi 
mu21=0.5+0.1*randn;mu22=1*0.5+0.1*randn;mu23=0;mu24=0; % TF2-TGi
mu31=0;mu32=0;mu33=0.5+0.1*randn;mu34=-1*(0.5+0.1*randn); % TF3-TGi
%MM coefficient
ktt11=rand;ktt12=rand;ktt13=1;ktt14=1; % TF1-TGi
ktt21=rand;ktt22=rand;ktt23=1;ktt24=1; % TF2-TGi
ktt31=1;ktt32=1;ktt33=rand;ktt34=rand; % TF3-TGi

gama1=rand;gama2=rand;gama3=rand;gama4=rand; % degradation rate of TGs

paraM2=zeros(7,4); % TFi-TGi parameters
paraM2(1,1)=mu11;paraM2(2,1)=mu21;paraM2(3,1)=mu31;% Maximal reaction rate, i.e., v
paraM2(1,2)=mu12;paraM2(2,2)=mu22;paraM2(3,2)=mu32;
paraM2(1,3)=mu13;paraM2(2,3)=mu23;paraM2(3,3)=mu33;
paraM2(1,4)=mu14;paraM2(2,4)=mu24;paraM2(3,4)=mu34;
paraM2(4,1)=mu11;paraM2(5,1)=mu21;paraM2(6,1)=mu31;% Maximal reaction rate, i.e., v
paraM2(4,2)=ktt11;paraM2(5,2)=ktt21;paraM2(6,2)=ktt31;
paraM2(4,3)=ktt13;paraM2(5,3)=ktt23;paraM2(6,3)=ktt33;
paraM2(4,4)=ktt14;paraM2(5,4)=ktt24;paraM2(6,4)=ktt34;
paraM2(7,1)=gama1;paraM2(7,2)=gama2;paraM2(7,3)=gama3;paraM2(7,3)=gama4; % degradation rate of TFs

%% PDE numerical simulation 
D1=0.1*1e1; D2=0.2*1e1; D3=0.3*1e1; D4=0.1*1e1; D5=0.3*1e1; % diffusion coefficient of ligands
r11=abs(normrnd(0.5,0.1));r21=abs(normrnd(0.4,0.1));r31=abs(normrnd(0.3,0.1));r41=abs(normrnd(0.2,0.1));r51=abs(normrnd(0.1,0.1)); % prodcution/release rate of ligand by cells
% r12=abs(normrnd(0.1,0.01));r22=abs(normrnd(0.2,0.1));r32=abs(normrnd(0.3,0.1));r42=abs(normrnd(0.4,0.1));r52=abs(normrnd(0.5,0.1));
d1=0.1; d2=0.1; d3=0.1; d4=0.1; d5=0.1; % degradation coefficient

%% Steady-state of ligands (elliptic PDEs)
f1=r11*(SC>0);
f2=r21*(SC>0);
f3=r31*(SC>0);
f4=r41*(SC>0);
f5=r51*(SC>0);

L1=FDM_elliptic(D1,d1,f1);
L2=FDM_elliptic(D2,d2,f2);
L3=FDM_elliptic(D3,d3,f3);
L4=FDM_elliptic(D4,d4,f4);
L5=FDM_elliptic(D5,d5,f5);

%% Solve ODEs for simulating Dynamic expressions of TGs in cells

R1_activation=b11.*L1.*R1+b21.*L2.*R1+b31.*L3.*R1+b41.*L4.*R1+b51.*L5.*R1;
R2_activation=b12.*L1.*R2+b22.*L2.*R2+b32.*L3.*R2+b42.*L4.*R2+b52.*L5.*R2;
TF1=zeros(N,N);TF2=zeros(N,N);TF3=zeros(N,N);TG1=zeros(N,N);TG2=zeros(N,N);TG3=zeros(N,N);TG4=zeros(N,N);
            
f_TF1=tf1.*(alpha11.*R1_activation./(krt11+R1_activation)+alpha21.*R2_activation./(krt21+R2_activation));
f_TF2=tf2.*(alpha12.*R1_activation./(krt12+R1_activation)+alpha22.*R2_activation./(krt22+R2_activation));
f_TF3=tf3.*(alpha13.*R1_activation./(krt13+R1_activation)+alpha23.*R2_activation./(krt23+R2_activation));

y0=[0.1 0.1 0.1 0.1 0.1 0.1 0.25]*5;  
in_noise=0.1;

% save Gaussian White noise
GW_TF1=zeros(N,N);GW_TF2=zeros(N,N);GW_TF3=zeros(N,N);
GW_TG1=zeros(N,N);GW_TG2=zeros(N,N);GW_TG3=zeros(N,N);GW_TG4=zeros(N,N);

tic
for i=1:N
    for j=1:N
        if SC(i,j)~=0 && TC(i,j)~=0
            %fprintf('The value of i is %d.\n', i);
            %fprintf('The value of j is %d.\n', j);
            t=TC(i,j); tspan=[0 t];
            A=f_TF1(i,j);B=f_TF2(i,j);C=f_TF3(i,j);
            h=0.01;
            GWnoise=0.5*randn(length(h:h:max(tspan)),7);
            fprintf('The value of GWnoise is %.2f.\n', GWnoise);
            y=y0;
            for tt=h:h:max(tspan)
                tmp1 = sqrt(h)*in_noise*sqrt(abs(y)).*GWnoise(round(tt/h),:);
                %fprintf('The value of sqrt(abs(y)) is %d.\n', tmp1);
                y=y+odefcn(tt,y,A,B,C)*h+sqrt(h)*in_noise*sqrt(abs(y)).*GWnoise(round(tt/h),:);  %%% CLE  chemical langevan equation  %% in_noise is 1/sqrt(V), V=10000 % PMIDs: 21697125, 30936559, 31907445. 
                y(y<0)=0;
            end
            TF1(i,j)=y(1);TF2(i,j)=y(2);TF3(i,j)=y(3);
            TG1(i,j)=y(4);TG2(i,j)=y(5);TG3(i,j)=y(6);TG4(i,j)=y(7);
            GW_TF1(i,j) = GWnoise(1);GW_TF2(i,j) = GWnoise(2);GW_TF3(i,j) = GWnoise(3);
            GW_TG1(i,j) = GWnoise(4);GW_TG2(i,j) = GWnoise(5);GW_TG3(i,j) = GWnoise(6);GW_TG4(i,j) = GWnoise(7);
        end
    end
end 
toc
        
%% Spatial coordinates of cells
[x,y]=find(SC(:,:)==1);
SC_xy=[x,y];
Cell_xy=SC_xy;
%% Expression matrix
SC_num=sum(sum(SC(:,:)==1));
L1_2D=L1(:,:);L2_2D=L2(:,:);L3_2D=L3(:,:);L4_2D=L4(:,:);L5_2D=L5(:,:);
R1_2D=R1(:,:);R2_2D=R2(:,:);
TF1_2D=TF1(:,:);TF2_2D=TF2(:,:);TF3_2D=TF3(:,:);
TG1_2D=TG1(:,:);TG2_2D=TG2(:,:);TG3_2D=TG3(:,:);TG4_2D=TG4(:,:);
tf1_2D=tf1(:,:);tf2_2D=tf2(:,:);tf3_2D=tf3(:,:);
%% all result
EM=zeros(18,SC_num);
EM(1,:)=L1_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2))); % expression of L
EM(2,:)=L2_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
EM(3,:)=L3_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
EM(4,:)=L4_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
EM(5,:)=L5_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));

EM(6,:)=R1_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2))); % expression of R
EM(7,:)=R2_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));

EM(8,:)=tf1_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2))); % expression of TF
EM(9,:)=tf2_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
EM(10,:)=tf3_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));

EM(11,:)=TG1_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2))); % expression of TG
EM(12,:)=TG2_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
EM(13,:)=TG3_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
EM(14,:)=TG4_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));

EM(15,:)=TF1_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2))); % activation of TF
EM(16,:)=TF2_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
EM(17,:)=TF3_2D(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));

EM(18,:)=TC(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));  % latent time for each cell

%% GW results
GW=zeros(7,SC_num);
GW(1,:)=GW_TF1(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
GW(2,:)=GW_TF2(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
GW(3,:)=GW_TF3(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
GW(4,:)=GW_TG1(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
GW(5,:)=GW_TG2(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
GW(6,:)=GW_TG3(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
GW(7,:)=GW_TG4(sub2ind([100,100],SC_xy(:,1),SC_xy(:,2)));
%% plot gene expression vs latent time

figure,
set(gcf, 'PaperPosition', [0, 0, 35, 26]); 
set(gcf, 'PaperSize', [35, 26]); 

subplot(2, 2, 1);
y=EM(11,:);x=EM(18,:);
scatter(x, y, 5,'filled');  
title('TG1 expression');
xlabel('Latent time');
ylabel('TG1 expression');
%saveas(gcf, fullfile(save_path, 'TG1_expression_with_latent_time.png'));

subplot(2, 2, 2);
y=EM(12,:);x=EM(18,:);
scatter(x, y, 5, 'filled');  
title('TG2 expression');
xlabel('Latent time');
ylabel('TG2 expression');
%saveas(gcf, fullfile(save_path, 'TG2_expression_with_latent_time.png'));

subplot(2, 2, 3);
y=EM(13,:);x=EM(18,:);
scatter(x, y, 5,'filled'); 
title('TG3 expression');
xlabel('Latent time');
ylabel('TG3 expression');
%saveas(gcf, fullfile(save_path, 'TG3_expression_with_latent_time.png'));

subplot(2, 2, 4);
y=EM(14,:);x=EM(18,:);
scatter(x, y, 5,'filled'); 
title('TG4 expression');
xlabel('Latent time');
ylabel('TG4 expression');
saveas(gcf, fullfile(save_path, 'FourGenes_expression_with_latent_time.png'));

%% plot gene expression vs spatial location

x = Cell_xy(:,1);y = Cell_xy(:,2); 

figure,
set(gcf, 'PaperPosition', [0, 0, 35, 28]); 
set(gcf, 'PaperSize', [35, 28]); 

subplot(2, 2, 1);
values = EM(11,:); % TG1 expression vs spatial location
scatter(x, y, 10, values, 'filled');  
colorbar;
title('Spatial gene expression of TG1');
xlabel('X Axis');
ylabel('Y Axis');
%saveas(gcf, fullfile(save_path, 'Spatial_TG1_expression.png'));

subplot(2, 2, 2);
values = EM(12,:); % TG2 expression vs spatial location
scatter(x, y, 10, values, 'filled'); 
colorbar;
title('Spatial gene expression of TG2');
xlabel('X Axis');
ylabel('Y Axis');
%saveas(gcf, fullfile(save_path, 'Spatial_TG2_expression.png'));

subplot(2, 2, 3);
values = EM(13,:); % TG3 expression vs spatial location
scatter(x, y, 10, values, 'filled'); 
colorbar;
title('Spatial gene expression of TG3');
xlabel('X Axis');
ylabel('Y Axis');
%saveas(gcf, fullfile(save_path, 'Spatial_TG3_expression.png'));

subplot(2, 2, 4);
values = EM(14,:); % TG4 expression vs spatial location
scatter(x, y, 10, values, 'filled'); 
colorbar;
title('Spatial gene expression of TG4');
xlabel('X Axis');
ylabel('Y Axis');
saveas(gcf, fullfile(save_path, 'FourGenes_Spatial_expression.png'));
%% save simulation expression data

save(fullfile(save_path, 'EM.mat'), 'EM');
save(fullfile(save_path, 'paraM1.mat'), 'paraM1');
save(fullfile(save_path, 'paraM2.mat'), 'paraM2');

xlswrite(fullfile(save_path, 'ExpressionMatrix.xlsx'), EM);
xlswrite(fullfile(save_path, 'SpatialCoordinates.xlsx'), Cell_xy);
xlswrite(fullfile(save_path, 'LR-TF_parameters.xlsx'), paraM1);
xlswrite(fullfile(save_path, 'TF-TG_parameters.xlsx'), paraM2);
xlswrite(fullfile(save_path, 'GWnoise.xlsx'), GW);


    
