% Chuong trinh hieu chuan 6 vi tri mo rong cho iPhone
% Tac gia: Tran Trung Chuyen
function Calib6PosExt
clc, clear all  clf

%% Zup
dataFile = 'iPhone6Plus_20160109_101105 AM';
sensor = importdata(strcat(dataFile,'-acce.txt'),' ');
zUpAcce = sensor.data;
sensor = importdata(strcat(dataFile,'-gyro.txt'),' ');
zUpGyro = sensor.data;
sensor = importdata(strcat(dataFile,'-acce1.txt'),' ');
zUpAcce1 = sensor.data;
sensor = importdata(strcat(dataFile,'-gyro1.txt'),' ');
zUpGyro1 = sensor.data;
sensor = importdata(strcat(dataFile,'-gravity.txt'),' ');
zUpGravity = sensor.data;
%% Zdown
dataFile = 'iPhone6Plus_20160109_101443 AM';
sensor = importdata(strcat(dataFile,'-acce.txt'),' ');
zDownAcce = sensor.data;
sensor = importdata(strcat(dataFile,'-gyro.txt'),' ');
zDownGyro = sensor.data;
sensor = importdata(strcat(dataFile,'-acce1.txt'),' ');
zDownAcce1 = sensor.data;
sensor = importdata(strcat(dataFile,'-gyro1.txt'),' ');
zDownGyro1 = sensor.data;
sensor = importdata(strcat(dataFile,'-gravity.txt'),' ');
zDownGravity = sensor.data;

%% Xup
dataFile = 'iPhone6Plus_20160109_101756 AM';
sensor = importdata(strcat(dataFile,'-acce.txt'),' ');
xUpAcce = sensor.data;
sensor = importdata(strcat(dataFile,'-gyro.txt'),' ');
xUpGyro = sensor.data;
sensor = importdata(strcat(dataFile,'-acce1.txt'),' ');
xUpAcce1 = sensor.data;
sensor = importdata(strcat(dataFile,'-gyro1.txt'),' ');
xUpGyro1 = sensor.data;
sensor = importdata(strcat(dataFile,'-gravity.txt'),' ');
xUpGravity = sensor.data;
%% Xdown
dataFile = 'iPhone6Plus_20160109_102053 AM';
sensor = importdata(strcat(dataFile,'-acce.txt'),' ');
xDownAcce = sensor.data;
sensor = importdata(strcat(dataFile,'-gyro.txt'),' ');
xDownGyro = sensor.data;
sensor = importdata(strcat(dataFile,'-acce1.txt'),' ');
xDownAcce1 = sensor.data;
sensor = importdata(strcat(dataFile,'-gyro1.txt'),' ');
xDownGyro1 = sensor.data;
sensor = importdata(strcat(dataFile,'-gravity.txt'),' ');
xDownGravity = sensor.data;

%% Yup
dataFile = 'iPhone6Plus_20160109_102354 AM';
sensor = importdata(strcat(dataFile,'-acce.txt'),' ');
yUpAcce = sensor.data;
sensor = importdata(strcat(dataFile,'-gyro.txt'),' ');
yUpGyro = sensor.data;
sensor = importdata(strcat(dataFile,'-acce1.txt'),' ');
yUpAcce1 = sensor.data;
sensor = importdata(strcat(dataFile,'-gyro1.txt'),' ');
yUpGyro1 = sensor.data;
sensor = importdata(strcat(dataFile,'-gravity.txt'),' ');
yUpGravity = sensor.data;
%% Ydown
dataFile = 'iPhone6Plus_20160109_102655 AM';
sensor = importdata(strcat(dataFile,'-acce.txt'),' ');
yDownAcce = sensor.data;
sensor = importdata(strcat(dataFile,'-gyro.txt'),' ');
yDownGyro = sensor.data;
sensor = importdata(strcat(dataFile,'-acce1.txt'),' ');
yDownAcce1 = sensor.data;
sensor = importdata(strcat(dataFile,'-gyro1.txt'),' ');
yDownGyro1 = sensor.data;
sensor = importdata(strcat(dataFile,'-gravity.txt'),' ');
yDownGravity = sensor.data;

g = 9.80665;
f=100;          % tan so lay mau
i1 = 3*f;       % cat doan dau du lieu (bo di 3 giay dau)
i2 = 150*f+i1-1;  % cat doan cuoi du lieu (lay 2,5 phut)

Omega=deg2rad(15.0);    % Toc do quay trai dat rad
phi=deg2rad(21.072575); % Vi do khu vuc thuc nghiem (Truong DH Mo-Dia chat)

[bxa,Sxa,bxg,Sxg,axu,axd,gxu,gxd,grxu,grxd]=Cal(1,i1,i2,xUpAcce1,xDownAcce1,xUpGyro1,xDownGyro1,xUpGravity,xDownGravity);
[bya,Sya,byg,Syg,ayu,ayd,gyu,gyd,gryu,gryd]=Cal(2,i1,i2,yUpAcce1,yDownAcce1,yUpGyro1,yDownGyro1,yUpGravity,yDownGravity);
[bza,Sza,bzg,Szg,azu,azd,gzu,gzd,grzu,grzd]=Cal(3,i1,i2,zUpAcce1,zDownAcce1,zUpGyro1,zDownGyro1,zUpGravity,zDownGravity);

ylim1=-0.03; % limit plot 
ylim2=0.03;
figure('Name', 'Sensors: Raw Data');
%% Accelerometer Z Up
fig1(1) = subplot(4,3,1);
AccePlot(azu,azu,'(1) Z_{Up}');

%% Accelerometer Z Down
fig1(2) = subplot(4,3,4);
AccePlot(azd,azd,'(2) Z_{Down}');

%% Accelerometer X Up
fig1(3) = subplot(4,3,2);
AccePlot(axu,axu,'(3) X_{Up}');

%% Accelerometer X Down
fig1(4) = subplot(4,3,5);
AccePlot(axd,axd,'(4) X_{Down}');

%% Accelerometer Y Up
fig1(5) = subplot(4,3,3);
AccePlot(ayu,ayu,'(5) Y_{Up}');

%% Accelerometer Y Down
fig1(6) = subplot(4,3,6);
AccePlot(ayd,ayd,'(6) Y_{Down}');

ylim1=-0.05; % limit plot 
ylim2=0.05;
%% Gyroscope Z Up
fig2(1) = subplot(4,3,7);
GyroPlot(gzu,gzu,'(1) Z_{Up}');

%% Gyroscope Z Down
fig2(2) = subplot(4,3,10);
GyroPlot(gzd,gzd,'(2) Z_{Down}');

%% Gyroscope X Up
fig2(1) = subplot(4,3,8);
GyroPlot(gxu,gxu,'(3) X_{Up}');

%% Gyroscope X Down
fig2(2) = subplot(4,3,11);
GyroPlot(gxd,gxd,'(4) X_{Down}');

%% Gyroscope Y Up
fig2(1) = subplot(4,3,9);
GyroPlot(gyu,gyu,'(5) Y_{Up}');

%% Gyroscope Y Down
fig2(2) = subplot(4,3,12);
GyroPlot(gyd,gyd,'(6) Y_{Down}');

hold off


%% Calibrated
ylim1=-0.01; % limit plot 
ylim2=0.01;
for i=1:length(ayu)
    axu(i,:)=axu(i,:)-[0 bxa*abs(grxu(i,2)) bya*abs(grxu(i,3)) bza*abs(grxu(i,4))];
    axd(i,:)=axd(i,:)-[0 bxa*abs(grxd(i,2)) bya*abs(grxd(i,3)) bza*abs(grxd(i,4))];
    ayu(i,:)=ayu(i,:)-[0 bxa*abs(gryu(i,2)) bya*abs(gryu(i,3)) bza*abs(gryu(i,4))];
    ayd(i,:)=ayd(i,:)-[0 bxa*abs(gryd(i,2)) bya*abs(gryd(i,3)) bza*abs(gryd(i,4))];
    azu(i,:)=azu(i,:)-[0 bxa*abs(grzu(i,2)) bya*abs(grzu(i,3)) bza*abs(grzu(i,4))];
    azd(i,:)=azd(i,:)-[0 bxa*abs(grzd(i,2)) bya*abs(grzd(i,3)) bza*abs(grzd(i,4))];
    
    gxu(i,:)=gxu(i,:)-[0 bxg byg bzg];
    gxd(i,:)=gxd(i,:)-[0 bxg byg bzg];
    gyu(i,:)=gyu(i,:)-[0 bxg byg bzg];
    gyd(i,:)=gyd(i,:)-[0 bxg byg bzg];
    gzu(i,:)=gzu(i,:)-[0 bxg byg bzg];
    gzd(i,:)=gzd(i,:)-[0 bxg byg bzg];
end
figure('Name', 'Sensors: 6Pos Ext Calibrated Data');
%% Accelerometer Z Up
fig1(1) = subplot(4,3,1);
AccePlot(azu,azu,'(1) Z_{Up}');

%% Accelerometer Z Down
fig1(2) = subplot(4,3,4);
AccePlot(azd,azd,'(2) Z_{Down}');

%% Accelerometer X Up
fig1(3) = subplot(4,3,2);
AccePlot(axu,axu,'(3) X_{Up}');

%% Accelerometer X Down
fig1(4) = subplot(4,3,5);
AccePlot(axd,axd,'(4) X_{Down}');

%% Accelerometer Y Up
fig1(5) = subplot(4,3,3);
AccePlot(ayu,ayu,'(5) Y_{Up}');

%% Accelerometer Y Down
fig1(6) = subplot(4,3,6);
AccePlot(ayd,ayd,'(6) Y_{Down}');

ylim1=-0.01; % limit plot 
ylim2=0.01;
%% Gyroscope Z Up
fig2(1) = subplot(4,3,7);
GyroPlot(gzu,gzu,'(1) Z_{Up}');

%% Gyroscope Z Down
fig2(2) = subplot(4,3,10);
GyroPlot(gzd,gzd,'(2) Z_{Down}');

%% Gyroscope X Up
fig2(1) = subplot(4,3,8);
GyroPlot(gxu,gxu,'(3) X_{Up}');

%% Gyroscope X Down
fig2(2) = subplot(4,3,11);
GyroPlot(gxd,gxd,'(4) X_{Down}');

%% Gyroscope Y Up
fig2(1) = subplot(4,3,9);
GyroPlot(gyu,gyu,'(5) Y_{Up}');

%% Gyroscope Y Down
fig2(2) = subplot(4,3,12);
GyroPlot(gyd,gyd,'(6) Y_{Down}');

hold off


%% Apple Refined
ylim1=-0.03; % limit plot 
ylim2=0.03;

figure('Name', 'Sensors: Apple Refined Data');
%% Accelerometer Z Up
fig1(1) = subplot(4,3,1);
AccePlot(azu,zUpAcce(i1:i2,:),'(1) Z_{Up}');

%% Accelerometer Z Down
fig1(2) = subplot(4,3,4);
AccePlot(azd,zDownAcce(i1:i2,:),'(2) Z_{Down}');

%% Accelerometer X Up
fig1(3) = subplot(4,3,2);
AccePlot(axu,xUpAcce(i1:i2,:),'(3) X_{Up}');

%% Accelerometer X Down
fig1(4) = subplot(4,3,5);
AccePlot(axd,xDownAcce(i1:i2,:),'(4) X_{Down}');

%% Accelerometer Y Up
fig1(5) = subplot(4,3,3);
AccePlot(ayu,yUpAcce(i1:i2,:),'(5) Y_{Up}');

%% Accelerometer Y Down
fig1(6) = subplot(4,3,6);
AccePlot(ayd,yDownAcce(i1:i2,:),'(6) Y_{Down}');

ylim1=-0.01; % limit plot 
ylim2=0.01;
%% Gyroscope Z Up
fig2(1) = subplot(4,3,7);
GyroPlot(gzu,zUpGyro(i1:i2,:),'(1) Z_{Up}');

%% Gyroscope Z Down
fig2(2) = subplot(4,3,10);
GyroPlot(gzd,zDownGyro(i1:i2,:),'(2) Z_{Down}');

%% Gyroscope X Up
fig2(1) = subplot(4,3,8);
GyroPlot(gxu,xUpGyro(i1:i2,:),'(3) X_{Up}');

%% Gyroscope X Down
fig2(2) = subplot(4,3,11);
GyroPlot(gxd,xDownGyro(i1:i2,:),'(4) X_{Down}');

%% Gyroscope Y Up
fig2(1) = subplot(4,3,9);
GyroPlot(gyu,yUpGyro(i1:i2,:),'(5) Y_{Up}');

%% Gyroscope Y Down
fig2(2) = subplot(4,3,12);
GyroPlot(gyd,yDownGyro(i1:i2,:),'(6) Y_{Down}');

%% Calib
    function [ba,Sa,bg,Sg,au,ad,gu,gd,gru,grd]=Cal(id,i1,i2,aUp,aDown,gUp,gDown,gravityUp,gravityDown)
        mean_ba=[0 0 0];
        mean_Sa=[1 1 1];

        mean_bg=[0 0 0];
        mean_Sg=[1 1 1];
        index=0;
        for i=i1:i2
            index=index+1;    
            t=aUp(i,1)-aUp(i1,1);
            au(index,:)=aUp(i,:)-gravityUp(i,:);
            au(index,1)=t;
            gru(index,:)=gravityUp(i,:);
            gru(index,1)=t;

            t=aDown(i,1)-aDown(i1,1);
            ad(index,:)=aDown(i,:)-gravityDown(i,:);
            ad(index,1)=t;
            grd(index,:)=gravityDown(i,:);
            grd(index,1)=t;

            ba = (au(index,id+1)+ad(index,id+1))/2;
            Sa = (au(index,id+1)-ad(index,id+1)-2*g)/(2*g);
            mean_ba(id) = (mean_ba(id) * (index-1) + ba) / index;
            mean_Sa(id) = (mean_Sa(id) * (index-1) + Sa) / index;

            t=gUp(i,1)-gUp(i1,1);
            gu(index,:)=gUp(i,:);
            gu(index,1)=t;

            t=gDown(i,1)-gDown(i1,1);
            gd(index,:)=gDown(i,:);
            gd(index,1)=t;

            bg = (gu(index,id+1)+gd(index,id+1))/2;
            Sg = (gd(index,id+1)-gu(index,id+1)-2.0*Omega*sin(phi))/(2.0*Omega*sin(phi));
            mean_bg(id) = (mean_bg(id) * (index-1) + bg) / index;
            mean_Sg(id) = (mean_Sg(id) * (index-1) + Sg) / index;
        end
        ba=mean_ba(id);
        Sa=abs(mean_Sa(id));
        bg=mean_bg(id);
        Sg=abs(mean_Sg(id));
    end
%% Cac ham ve
    function AccePlot(time,y,t)
        ylim([ylim1 ylim2]);
        hold on
        grid on

        plot(time(:,1),y(:,2), 'r');
        plot(time(:,1),y(:,3), 'g');
        plot(time(:,1),y(:,4), 'b');

        legend('a_x','a_y','a_z','Location','NorthEastOutside');
        xlabel('Thoi gian (s)');
        ylabel('Gia toc (g)');
        title(strcat(t));
    end

    function GyroPlot(time,y,t)
        ylim([ylim1 ylim2]);
        hold on
        grid on

        plot(time(:,1),y(:,2), 'r');
        plot(time(:,1),y(:,3), 'g');
        plot(time(:,1),y(:,4), 'b');

        legend('g_x','g_y','g_z','Location','NorthEastOutside');
        xlabel('Thoi gian (s)');
        ylabel('Van toc goc (rad/s)');
        title(strcat(t));
    end

%linkaxes(fig2, 'y');
hold off
end