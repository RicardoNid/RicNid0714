%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: 生成训练序列，接收端对其进行处理，估计信道，进而做频域均衡   
%  // ======================================================================

function preamble = CreateOFDMPreamble(OFDMParameters)

bitsNumber = length(OFDMParameters.PreambleCarrierPositions) * OFDMParameters.PreambleBitsPerSymbolQAM;
preambleBits = randint(bitsNumber, 1, 2, OFDMParameters.PreambleSeed);
preambleQAMSymbols = GrayQAMCoder(preambleBits, OFDMParameters.PreambleBitsPerSymbolQAM);
preambleQAMSymbols = preambleQAMSymbols./rms(preambleQAMSymbols);

ifftBlock = zeros(OFDMParameters.FFTSize, 1);
ifftBlock(OFDMParameters.PreambleCarrierPositions) = preambleQAMSymbols;
ifftBlock(OFDMParameters.FFTSize+2 - OFDMParameters.PreambleCarrierPositions) = conj(preambleQAMSymbols);

sig = ifft(ifftBlock);
CPLength = OFDMParameters.CPLength;
sig = [sig(end-CPLength/2+1:end); sig; sig(1:CPLength/2)];
preamble = repmat(sig, OFDMParameters.PreambleNumber, 1); %将sig，重复4次，[sig;sig;sig;sig]

