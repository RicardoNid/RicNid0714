%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: 主函数      
%  // ======================================================================
clear all
close all
clc
tic
BER_SNR = [];
BER_MC_Sim_total = [];
warning off all
OFDMParameters = InitOFDMParameters();
on = OFDMParameters.on;
bitNumber = OFDMParameters.bitNumber;
%20个子帧组成一个大的帧，一个子帧里有2训练序列，8个FTN符号

% [OFDMFrame,OFDMSymbols] = OFDMFrameGenerator(OFDMParameters);
OFDMFrame = OFDMBigFrameGenerator( OFDMParameters );
OFDMFrame = filter([1,0.8,0.1,0.05,0.01,0.005], 1, OFDMFrame);
OFDMFrameplot = 20*log10(abs(fft(OFDMFrame)));
figure;plot(OFDMFrameplot);
% scatterplot(fft(frame))
% p = rms(OFDMFrame);
%% EbNo 仿真
BER_Pre = [];
BER_post = [];
SNR = 20;
BER_MC_Sim = zeros(length(SNR),1);
snr=10^(SNR/10);
code_power=norm(OFDMFrame)^2/(length(OFDMFrame)); %信号的符号功率 =var(passchan_ofdm_symbol)
sigma=sqrt(code_power/(snr*2));%sigma如何计算，与当前SNR和信号平均能量有关系
[OFDMFrame_rec,awgn]=addnoise(OFDMFrame,sigma); % use randn 噪声只加了实部
%figure();plot(20*log10(abs(fft(OFDMFrame_rec))))
%% 接收端
%把1个大帧拆成20个子帧
OFDMFrame_total = reshape(OFDMFrame_rec, [],20);
number_of_error_all = 0;
pre_code_errors_all = 0;

for cir = 1:20
    OFDMFrame_rec_per = OFDMFrame_total(:,cir);
    [PerQAMtotal,PerQAMError,QAM_re_sum, pre_code_errors,number_of_error] = OFDMFrameReceiver(OFDMFrame_rec_per, OFDMParameters, cir);
    number_of_error_all = number_of_error_all + number_of_error;
    pre_code_errors_all = pre_code_errors_all + pre_code_errors;
    %收集20个子帧的点，画星座图, 对应的每个星座图的误码率是纠后的
    QAM_re_all(cir,:) = QAM_re_sum;
    PerQAMError_all(cir,:) = PerQAMError;
    PerQAMtotal_all(cir,:) = PerQAMtotal;
end

if on == 0
    BER_pre_woChow = pre_code_errors_all/(bitNumber*20) %bitNumber是一个小子帧的比特数，20个小子帧组成一个大帧
    BER_post_woChow = number_of_error_all/(bitNumber*20)
else    
  load('bitAllocSort.mat');
    BER_pre_wChow = pre_code_errors_all/(bitNumber*20) %bitNumber是一个小子帧的比特数，20个小子帧组成一个大帧
    BER_post_wChow = number_of_error_all/(bitNumber*20)
    PerQAMError_all = sum(PerQAMError_all);
    PerQAMtotal_all = sum(PerQAMtotal_all);
    for i = 1:length(bitAllocSort)
        QAM_re_per_total = QAM_re_all(:,i);
        QAM_re_per_total = cell2mat(QAM_re_per_total);
        scatterplot(QAM_re_per_total); 
        %取出每种QAM对应的错误数和比特总数
        errorNum = PerQAMError_all(1,i);
        totalNum = PerQAMtotal_all(1,i);
        BER_post = errorNum/totalNum;
        title(['BER-MC-post = ', num2str(BER_post)  '---NE= ', num2str(errorNum)]);
    end
end
toc