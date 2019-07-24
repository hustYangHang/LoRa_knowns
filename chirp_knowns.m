% chirp信号产生
clear;clc;close all;

BW = 500e3;  % 信号带宽
SF = 12;  % LoRa信号SF
F0 = 0;  % 信号中心频率

simple_time = 2^SF/BW;  % 一个LoRa符号周期
t = 1/BW:1/BW:simple_time;  % 仿真时间

f0 = -BW/2 + F0;  % 下边带
f1 = BW/2 + F0;  % 上边带
%  upChirp 2路正交信号
up_chirpI = chirp(t,f0,simple_time,f1,'linear',90);
up_chirpQ = chirp(t,f0,simple_time,f1,'linear',0);
up_chirp = complex(up_chirpI, up_chirpQ);
%  downChirp 2路正交信号
down_chirpI = chirp(t,f1,simple_time,f0,'linear',90);
down_chirpQ = chirp(t,f1,simple_time,f0,'linear',0);
down_chirp = complex(down_chirpI, down_chirpQ);

%  信号调制
bits_map = [2^0, 2^1, 2^2, 2^3, 2^4, 2^5, 2^6, 2^7, 2^8, 2^9, 2^10, 2^ 11];
infor_bits = [0,0,0,0,0,0,0,0,0,0,0,1];
shift_num = dot(infor_bits, bits_map(1,1:SF));  % 计算起始频率的位置

x_signal_0 = circshift(down_chirp, shift_num);  % 符号1：循环移位，使信号从起始频率开始发送
x_signal_1 = down_chirp;   % 符号2： 发送信息为[0,0,0,0,0,0,0,0,0,0,0,0],起始频率不变
%  解调
y = [x_signal_0, x_signal_1];  % 接收信号
de_chirp = repmat(up_chirp, 1, ceil(length(y) / length(up_chirp)));  % 解调信号
signal = y.*de_chirp;
symbol_num = length(y)/2^SF;
for i = 1:symbol_num
    s(:,i) = fft(signal(1,((i - 1) * 2^SF) + 1:i * 2^SF));
end

s_matrix = abs(s);
[~, symbols_recovered] = max(abs(s));
symbols_recovered = symbols_recovered-1;  % 得到起始频率点
% 根据起始频率点恢复bit信息
bit_recovered = 0;
for i = 1 : length(symbols_recovered)
    bit_recovered = [ bit_recovered , bitget(symbols_recovered(i), 1:SF) ];
end
bit_recovered = bit_recovered(2:end);

% 查看chirp信号的频谱
[S,F,T,P] = spectrogram(y,32,30,64,BW);
surf(T,F,10*log10(P),'edgecolor','none'); axis tight;
view(0,90);
xlabel('Time (Seconds)'); ylabel('Hz');
% 查看解调信号的FFT结果
figure;
subplot(2,1,1);
plot(s_matrix(:,1));
title('symbol 1');
set(gca, 'Fontname', 'Times New Roman');
subplot(2,1,2);
plot(s_matrix(:,2));
title('symbol 2');
set(gca, 'Fontname', 'Times New Roman');
% 查看chirp信号的时域
figure;
plot(real(up_chirp(1,2049:2049+1000)))



