%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: ¹À¼ÆÐÅµÀ 
%  // ======================================================================

function [H] = ChannelEstimationByPreamble(recvPreamble, OFDMParameters)
%% parameters
PreambleBitsPerSymbolQAM = OFDMParameters.PreambleBitsPerSymbolQAM;
PreambleCarrierPositions = OFDMParameters.PreambleCarrierPositions;
FFTSize = OFDMParameters.FFTSize;
CPLength = OFDMParameters.CPLength;
numCarrier = length(OFDMParameters.PreambleCarrierPositions);
PreambleNumber = OFDMParameters.PreambleNumber;

bitsNumber = numCarrier*PreambleBitsPerSymbolQAM;
preambleBits = randint(bitsNumber, 1, 2, OFDMParameters.PreambleSeed);
GQAMSymbols = GrayQAMCoder(preambleBits, PreambleBitsPerSymbolQAM);

recvPreambleSignal = reshape(recvPreamble, FFTSize+CPLength, []);
recvPreambleSignal = recvPreambleSignal(CPLength/2+1:end-CPLength/2,:);

recvQAMSignal = fft(recvPreambleSignal);
recvQAMSignal = recvQAMSignal(PreambleCarrierPositions,:);

H_temp = zeros(numCarrier, PreambleNumber);

for i = 1:PreambleNumber
    H_temp(:, i) = recvQAMSignal(:,i)./GQAMSymbols*sqrt(10);
end
H = mean(H_temp, 2);
