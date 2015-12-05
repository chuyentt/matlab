%% Nhap du lieu tu phan mem SDLogger
%% Update intervals 1Hz
Gyroscope_1Hz = importdata('iPhone4_20151204_142004-gyro.txt',' ');
Gyroscope1_1Hz = importdata('iPhone4_20151204_142004-gyro1.txt',' ');
Accelerometer_1Hz = importdata('iPhone4_20151204_142004-acce.txt',' ');
Accelerometer1_1Hz = importdata('iPhone4_20151204_142004-acce1.txt',' ');

%% Update intervals 10Hz
Gyroscope_10Hz = importdata('iPhone4_20151204_150504-gyro.txt',' ');
Gyroscope1_10Hz = importdata('iPhone4_20151204_150504-gyro1.txt',' ');
Accelerometer_10Hz = importdata('iPhone4_20151204_150504-acce.txt',' ');
Accelerometer1_10Hz = importdata('iPhone4_20151204_150504-acce1.txt',' ');

%% Update intervals 100Hz
Gyroscope_100Hz = importdata('iPhone4_20151204_155005-gyro.txt',' ');
Gyroscope1_100Hz = importdata('iPhone4_20151204_155005-gyro1.txt',' ');
Accelerometer_100Hz = importdata('iPhone4_20151204_155005-acce.txt',' ');
Accelerometer1_100Hz = importdata('iPhone4_20151204_155005-acce1.txt',' ');

%% Gan che do test voi tung tan so
Gyroscope = Gyroscope_100Hz;
Gyroscope1 = Gyroscope1_100Hz;
Accelerometer = Accelerometer_100Hz;
Accelerometer1 = Accelerometer1_100Hz;

%% Ve bieu dien du lieu Gyroscope
figure('Name', 'Sensors Data');
axis(1) = subplot(4,1,1);
hold on;
plot(Gyroscope1.data(:,2), 'r');
plot(Gyroscope1.data(:,3), 'g');
plot(Gyroscope1.data(:,4), 'b');
legend('X', 'Y', 'Z');
xlabel('Thoi gian (s)');
ylabel('Toc do goc (rad/s)');
line1 = strrep(Gyroscope1.textdata(2), '// Device Info:', '');
line2 = strrep(Gyroscope1.textdata(3), '// Update intervals: ', ' Frequencies');
title(strcat('Gyroscope: ', line1, line2));
hold off;

%% Ve bieu dien du lieu Gyroscope
axis(2) = subplot(4,1,2);
hold on;
plot(Gyroscope.data(:,2), 'r');
plot(Gyroscope.data(:,3), 'g');
plot(Gyroscope.data(:,4), 'b');
legend('X', 'Y', 'Z');
xlabel('Thoi gian (s)');
ylabel('Toc do goc (rad/s)');
line1 = strrep(Gyroscope.textdata(2), '// Device Info:', '');
line2 = strrep(Gyroscope.textdata(3), '// Update intervals: ', ' Frequencies');
title(strcat('Gyroscope (Calibrated): ', line1, line2));
hold off;

%% Ve bieu dien du lieu Accelerometer1
axis(3) = subplot(4,1,3);
hold on;
plot(Accelerometer1.data(:,2), 'r');
plot(Accelerometer1.data(:,3), 'g');
plot(Accelerometer1.data(:,4), 'b');
legend('X', 'Y', 'Z');
xlabel('Thoi gian (s)');
ylabel('Gia toc (g)');
line1 = strrep(Accelerometer1.textdata(2), '// Device Info:', '');
line2 = strrep(Accelerometer1.textdata(3), '// Update intervals: ', ' Frequencies');
title(strcat('Accelerometer: ', line1, line2));
hold off;

%% Ve bieu dien du lieu Accelerometer
axis(4) = subplot(4,1,4);
hold on;
plot(Accelerometer.data(:,2), 'r');
plot(Accelerometer.data(:,3), 'g');
plot(Accelerometer.data(:,4), 'b');
legend('X', 'Y', 'Z');
xlabel('Thoi gian (s)');
ylabel('Gia toc (g)');
line1 = strrep(Accelerometer.textdata(2), '// Device Info:', '');
line2 = strrep(Accelerometer.textdata(3), '// Update intervals: ', ' Frequencies');
title(strcat('Accelerometer (Calibrated): ', line1, line2));
hold off;
linkaxes(axis, 'x');
