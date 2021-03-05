%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: QAM”≥…‰     
%  // ======================================================================
function QAMSymbols = GrayQAMCoder(bits, BitsPerSymbol)

if length(bits) == 0
    QAMSymbols = [];
    return
end
if BitsPerSymbol < 1 || BitsPerSymbol ~= int8(BitsPerSymbol)
    error('BitsPerSymbol should be a positive integer');
end

bits = bits(:);
BitNumber = length(bits);
if rem(BitNumber, BitsPerSymbol) ~= 0 
    error('BitNumber must dividable by BitsPerSymbol')
end

ConstellationSize = 2^BitsPerSymbol;
symbols = bi2de(reshape(bits, BitsPerSymbol, []).','left-msb');
symbols = bin2gray(symbols, 'qam', ConstellationSize);
QAMSymbols = modulate(modem.qammod(ConstellationSize), symbols); 