%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: QAMΩ‚”≥…‰      
%  // ======================================================================
function bits = GrayQAMDecoder(QAMSymbols, BitsPerSymbol)

if BitsPerSymbol < 1 || BitsPerSymbol ~= int8(BitsPerSymbol)
    error('BitsPerSymbol should be a positive integer');
end

ConstellationSize = 2^BitsPerSymbol;
symbols = demodulate(modem.qamdemod(ConstellationSize), QAMSymbols);
symbols = gray2bin(symbols, 'qam', ConstellationSize);
bits = de2bi(symbols,BitsPerSymbol,'left-msb');
bits = reshape(bits.',[], 1);