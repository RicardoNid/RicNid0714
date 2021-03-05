function bits = PAMDecoder(Symbols, BitsPerSymbol)

if BitsPerSymbol < 1 || BitsPerSymbol ~= int8(BitsPerSymbol)
    error('BitsPerSymbol should be a positive integer');
end

M = 2^BitsPerSymbol;
symbols = demodulate(modem.pamdemod('M',M), Symbols);
symbols = gray2bin(symbols, 'pam', M);
bits = de2bi(symbols,BitsPerSymbol,'left-msb');
bits = reshape(bits.',[], 1);