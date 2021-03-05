function PAMSymbols = PAMCoder(bits, BitsPerSymbol)

if length(bits) == 0
    PAMSymbols = [];
    return
end
if BitsPerSymbol < 1 || BitsPerSymbol ~= int8(BitsPerSymbol) %int8 -127~128
    error('BitsPerSymbol should be a positive integer');
end

bits = bits(:);%one column
BitNumber = length(bits);
if rem(BitNumber, BitsPerSymbol) ~= 0 %mod
    error('BitNumber must dividable by BitsPerSymbol')
end

M = 2^BitsPerSymbol;
symbols = bi2de(reshape(bits, BitsPerSymbol, []).','left-msb');%binary->decimal
symbols = bin2gray(symbols, 'pam', M);%gray coding position
PAMSymbols = modulate(modem.pammod('M',M), symbols);
        