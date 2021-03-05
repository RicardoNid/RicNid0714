%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: 做频域均衡
%  // ======================================================================
function [recoveredQAMSymbols, OFDMSymbols] = RecoverOFDMSymbolsWithPilot(recvOFDMSignal, OFDMParameters, H)
OFDMPositions = OFDMParameters.OFDMPositions;
OFDMSymbolNumber = OFDMParameters.OFDMSymbolNumber;
DataCarrierPositions = OFDMParameters.DataCarrierPositions;
FFTSize = OFDMParameters.FFTSize;
CPLength = OFDMParameters.CPLength;
recvOFDMSignal = reshape(recvOFDMSignal, length(OFDMPositions)+CPLength, []);
recvOFDMSignal = recvOFDMSignal(CPLength/2+1:end-CPLength/2, :);
% 插0，时域尾部去零
recvOFDMSignal_interp = zeros(FFTSize,OFDMSymbolNumber*2);
recvOFDMSignal_interp(1:length(OFDMPositions),:)= recvOFDMSignal; 
recvOFDMSignal = fft(recvOFDMSignal_interp);
dataQAMSymbols = recvOFDMSignal(DataCarrierPositions, :);
DataCarrierPositions_H = DataCarrierPositions;
    for i = 1:OFDMParameters.OFDMSymbolNumber*2
        recoveredQAMSymbols(:, i) = dataQAMSymbols(:, i) ./ H(DataCarrierPositions_H);
    end
 
%%     
ifftBlock = zeros(FFTSize, OFDMSymbolNumber*2); 
ifftBlock(DataCarrierPositions, :) = recoveredQAMSymbols; 
ifftBlock(FFTSize + 2 - DataCarrierPositions, :) = conj(ifftBlock(DataCarrierPositions,:));
OFDMSymbols = ifft(ifftBlock);
OFDMSymbols1 = OFDMSymbols(1:length(OFDMPositions),:);
OFDMSymbols1 = [OFDMSymbols1(end-CPLength/2+1:end, :); OFDMSymbols1; OFDMSymbols1(1:CPLength/2, :)];
OFDMSymbols = reshape(OFDMSymbols1, [], 1);
end
