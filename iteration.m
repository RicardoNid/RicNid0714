%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: on=0,计算比特分配时的迭代，取最后一个子帧的最后一次迭代
%  //                结果进行比特分配计算      
%  // ======================================================================
function [pre_code_errors,pre_code_ber,BER_MC_HD,decodedMsg_HD,number_of_error_HD]=iteration(decodedMsg_HD,OFDMParameters,recvOFDMSignal, tblen,i,recoveredSymbols_FDE, cir)
FFTSize = OFDMParameters.FFTSize;
BitsPerSymbolQAM = OFDMParameters.BitsPerSymbolQAM;
CPLength = OFDMParameters.CPLength;
OFDMSymbolNumber = OFDMParameters.OFDMSymbolNumber;
DataCarrierPositions=OFDMParameters.DataCarrierPositions;
OFDMPositions = OFDMParameters.OFDMPositions;
SubcarriersNum = length(OFDMParameters.DataCarrierPositions);
iterationT = OFDMParameters.iteration;
SToPcol = OFDMParameters.SToPcol;
%% channel coding
% Code properties
codeRate = 1/2;
constlen = 7;
codegen  = [171 133];
trellis  = poly2trellis(constlen, codegen);
codedMsg  = convenc(decodedMsg_HD, trellis);

%% interleaving
depth=32;
codedMsg=reshape(codedMsg,depth,[]);
codedMsg1=[];
for k=1:depth
    codedMsg1=[codedMsg1,codedMsg(k,:)];
end
codedMsg1=codedMsg1.';

%% mapping
M=2^BitsPerSymbolQAM;
modObj   =modem.qammod('M', M, 'SymbolOrder', 'Gray', 'InputType', 'Bit');
QAMSymbols = modulate(modObj, codedMsg1);
QAMSymbols = QAMSymbols/rms(QAMSymbols)*sqrt(10);
QAMSymbols = reshape(QAMSymbols, length(DataCarrierPositions), SToPcol);

%% 结构简单，算ICI
%% 保存S'HD
S_HD = reshape(QAMSymbols, [], 1);
save S_HD S_HD
%% IFFT(重复发端操作）zero padding
ifftBlock = zeros(FFTSize, SToPcol);
ifftBlock(DataCarrierPositions, :) = QAMSymbols;
ifftBlock(FFTSize + 2 - DataCarrierPositions, :) = conj(ifftBlock(DataCarrierPositions,:));
OFDMSymbols = ifft(ifftBlock);
OFDMSymbols1 = OFDMSymbols(1:length(OFDMPositions),:);
OFDMSymbols1 = [OFDMSymbols1(end-CPLength/2+1:end, :); OFDMSymbols1; OFDMSymbols1(1:CPLength/2, :)];
OFDMSymbols = reshape(OFDMSymbols1, [], 1);
%% FFT
recovered = RecoverOFDMSymbols(OFDMSymbols, OFDMParameters);
recoveredSymbols = reshape(recovered,[],1);
recoveredSymbols = recoveredSymbols/rms(recoveredSymbols)*sqrt(10);%16QAM记得改
%% ICI
ICI = recoveredSymbols - S_HD;
recoveredSymbols = recoveredSymbols_FDE - ICI;
%%
% Code properties
codeRate = 1/2;
constlen = 7;
codegen  = [171 133];
trellis  = poly2trellis(constlen, codegen);

%% de-mapping
M=2^BitsPerSymbolQAM;
modObj   =modem.qammod('M', M, 'SymbolOrder', 'Gray', 'InputType', 'Bit');
demodObj =modem.qamdemod(modObj);
% Set up the demodulator object to perform hard decision demodulation
set(demodObj, 'DecisionType', 'Hard decision');
demodulatedMsg_HD = demodulate(demodObj, recoveredSymbols);
demodulatedMsg_HD=demodulatedMsg_HD';

load('codedMsg1.mat');
[pre_code_errors, pre_code_ber] = biterr(demodulatedMsg_HD', codedMsg1);
%% de-interleaving
depth=32;
len=length(demodulatedMsg_HD)/depth;
codedMsg1=[];
for k=1:depth
    codedMsg1=[codedMsg1;demodulatedMsg_HD(len*(k-1)+1:len*k)];
end
demodulatedMsg_HD=codedMsg1(:);
%% Use the Viterbi decoder in hard decision mode
decodedMsg_HD = vitdec(demodulatedMsg_HD, trellis, tblen, 'cont', 'hard');
% Compute the bit error rate
load('bits.mat');
bits=bits.';
decodedMsg_HD=[decodedMsg_HD(tblen+1:end);bits(length(bits)-tblen+1:length(bits))];

% output
[nErrors_HD, ber_HD]   = biterr(decodedMsg_HD, bits);
% figure;plot(xor(decodedMsg_HD,bits));
BER_MC_HD=ber_HD;
number_of_error_HD = nErrors_HD;
RecoveredSymbols = recoveredSymbols(:);
% statplot(RecoveredSymbols)
% title([num2str(i),'iteration','BER-MC = ', num2str(BER_MC_HD)])

%% Calculate the optimal bit allocation
if cir == 20
if i==iterationT
    QAMSymbols_trans = cell2mat(struct2cell(load('QAMSymbols_trans.mat')));
    SendSymbols = QAMSymbols_trans*sqrt(10);
    SendSymbols = reshape(SendSymbols,1,[]);
    recoveredQAMSymbols = reshape(recoveredSymbols,1,[]);
    ErrorLocation(recoveredQAMSymbols,SendSymbols,OFDMParameters);
    title([num2str(i),'iteration','ErrorPersubcarrier '])
    SNR = SNRLocation(recoveredQAMSymbols,SendSymbols,OFDMParameters);
    title([num2str(i),'iteration','SNRPersubcarrier '])
    
    BER=1e-3;%误码率
    %计算gap
%     gap=-log(5*BER)/1.5; %in dB  固定公式，某种调整参数
    SER = 1-(1-BER)^4;
    gap = 1/3*(qfuncinv(SER/4))^2;
    target_bit = 4;
    Rbit = round(target_bit*SubcarriersNum);
    miu = 1e-5;
    [bits_allo, power_allo,total_bits]=chow_algo_all(SNR,SubcarriersNum,gap,Rbit);
    bits_alloc = bits_allo;
    power_alloc = power_allo;
    
    if total_bits == 0
        flag = 1;
        while (BER >0) && (BER < 0.2) && (total_bits == 0)
            BER = BER + flag * miu;
            SER = 1-(1-BER)^4;
            % gap=-log(5*BER)/1.5;
            gap = 1/3*(qfuncinv(SER/4))^2;
            bits_alloc_record = bits_allo;
            power_alloc_record = power_allo;
            [bits_allo, power_allo, total_bits]=chow_algo_all(SNR, SubcarriersNum, gap, Rbit);
        end
    end
    if total_bits ~= 0
        flag = -1;
        while (BER > 0) && (BER < 0.2) && (total_bits > 0)
            BER = BER + flag * miu;
            SER = 1-(1-BER)^4;
            %  gap=-log(5*BER)/1.5;
            gap = 1/3*(qfuncinv(SER/4))^2;
            bits_alloc_record = bits_allo;
            power_alloc_record = power_allo;
            [bits_allo, power_allo, total_bits]=chow_algo_all(SNR, SubcarriersNum, gap, Rbit);
        end
    end
    BER = BER- flag * miu;
    disp(BER)
    
    bits_alloc = bits_alloc_record;
    power_alloc = power_alloc_record;
  
    figure; plot(bits_alloc);
    xlabel('subcarrier');
    ylabel('bits alloc');
    figure;plot(10*log10(abs(power_alloc)));%取dB
    figure; plot(power_alloc);%不取dB
    xlabel('subcarrier');
    ylabel('power alloc');
    bits_alloc = bits_alloc';
    power_alloc = power_alloc';
    sum(power_alloc);
    figure;plot(10*log10(abs(power_alloc)));
    % Bit allocation subcarrier summary
    [bitAllocSort,BitAllocSum] =bits_alloc_position_sum(bits_alloc,SubcarriersNum);
    save bitAllocSort bitAllocSort;
    save BitAllocSum BitAllocSum;
    save power_alloc power_alloc;
end
end