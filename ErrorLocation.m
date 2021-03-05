%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: º∆À„¥ÌŒÛ∑÷≤º      
%  // ======================================================================
function ErrorPerSubcarrier = ErrorLocation(recoveredSymbols, transmittedSymbols, OFDMParameters)
SubcarriersNum = length(OFDMParameters.DataCarrierPositions);
recoveredSymbols = reshape(recoveredSymbols,SubcarriersNum,[]);
transmittedSymbols = reshape(transmittedSymbols,SubcarriersNum,[]);
BitsPerSymbolQAM = OFDMParameters.BitsPerSymbolQAM;
number_of_error = zeros(SubcarriersNum,1);
for i = 1:SubcarriersNum
    recvBits = GrayQAMDecoder(recoveredSymbols(i,:),BitsPerSymbolQAM);
    sendBits = GrayQAMDecoder(transmittedSymbols(i,:),BitsPerSymbolQAM);
    [number_of_error(i), BER_MC] = biterr(sendBits, recvBits);
end
ErrorPerSubcarrier = number_of_error;
figure; plot(ErrorPerSubcarrier);
end