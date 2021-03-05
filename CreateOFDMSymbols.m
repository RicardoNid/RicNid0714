%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: 生成FTN符号      
%  // ======================================================================
function OFDMSymbols = CreateOFDMSymbols(OFDMParameters)
on = OFDMParameters.on;
FFTSize = OFDMParameters.FFTSize;
OFDMPositions = OFDMParameters.OFDMPositions;
CPLength = OFDMParameters.CPLength;
OFDMSymbolNumber = OFDMParameters.OFDMSymbolNumber;
BitsPerSymbolQAM = OFDMParameters.BitsPerSymbolQAM;
DataCarrierPositions = OFDMParameters.DataCarrierPositions;
bitNumber = OFDMParameters.bitNumber;
SToPcol = OFDMParameters.SToPcol;

if on ==1
    %% bit loading %%
    load('bitAllocSort.mat');
    load('BitAllocSum.mat');
    load('power_alloc.mat');
    rmsAlloc = [];
    bits = [];
    QAMSymbols_trans = [];
    ifftBlock = zeros(FFTSize, SToPcol);
    for i=1:length(bitAllocSort)
        carrierPosition = BitAllocSum{i};
        bitNumber = OFDMSymbolNumber * length(carrierPosition) * bitAllocSort(i);
        bits_per = randint(bitNumber, 1, 2, OFDMParameters.Seed);
        bits = [bits;bits_per];
    end
    save('bits.mat','bits');
    
    % Code properties(channel coding)
    codeRate = 1/2;
    constlen = 7;
    codegen  = [171 133];
    trellis  = poly2trellis(constlen, codegen);
    codedMsg  = convenc(bits, trellis);
    
    % Conv(interleaving)
    depth = 32;
    codedMsg=reshape(codedMsg,depth,[]);
    codedMsg1=[];
    for k=1:depth
        codedMsg1=[codedMsg1,codedMsg(k,:)];
    end
    codedMsg1=codedMsg1.';
    save codedMsg1 codedMsg1;
    
    % mapping（自带的qammod)
    b=1;
    for i=1:length(bitAllocSort)
        if  bitAllocSort(i) == 0
            QAMSymbols=0;
            rmsAlloc=[0];
            %rmsAlloc=[rmsAlloc;0];
        else
            M=2^bitAllocSort(i);
            modObj   =modem.qammod('M', M, 'SymbolOrder', 'Gray', 'InputType', 'Bit');
            codeMsg1_per = OFDMSymbolNumber * bitAllocSort(i)*length(BitAllocSum{i})*2;
            codeMsg1_perloading = codedMsg1(b:b+codeMsg1_per-1,1);
            b = codeMsg1_per+b;
            if bitAllocSort(i) == 3 % QAM8
                QAM8 = [-1-sqrt(3), -1+1i, -1-1i, 1i*(1+sqrt(3)), -1i*(1+sqrt(3)), 1+1i, 1-1i, 1+sqrt(3)];
                qam8bit = reshape(codeMsg1_perloading, bitAllocSort(i), [])';
                qam8dec = bi2de(qam8bit,'left-msb');
                QAMSymbols = QAM8(qam8dec + 1);
                QAMSymbols = QAMSymbols';
                % scatterplot(QAMSymbols);
            else
                QAMSymbols = modulate(modObj, codeMsg1_perloading);
                % scatterplot(QAMSymbols);
            end
            rms_alloc=rms(QAMSymbols);
            rmsAlloc=[rmsAlloc;rms_alloc];
            QAMSymbols =QAMSymbols/rms_alloc;
            QAMSymbols_trans = [QAMSymbols_trans;QAMSymbols];
            QAMSymbols = reshape(QAMSymbols, length(BitAllocSum{i}), SToPcol);
        end
        carrierPosition = BitAllocSum{i};
        carrierPosition = carrierPosition+2;
        ifftBlock(carrierPosition, :) = QAMSymbols;
    end
    save rmsAlloc  rmsAlloc;
    save QAMSymbols_trans  QAMSymbols_trans;
    
    % 功率加载
    for i = 1:SToPcol
        ifftBlock(DataCarrierPositions,i) = ifftBlock(DataCarrierPositions,i).*sqrt(power_alloc');
    end
    %%
else
   %% cal %%
    bits = randint(1, bitNumber, 2, OFDMParameters.Seed);
    save('bits.mat','bits');
    
    % Code properties(channel coding)
    codeRate = 1/2;
    constlen = 7;
    codegen  = [171 133];
    trellis  = poly2trellis(constlen, codegen);
    codedMsg  = convenc(bits, trellis);
    
    % Conv(interleaving)
    depth = 32;
    codedMsg=reshape(codedMsg,depth,[]);
    codedMsg1=[];
    for k=1:depth
        codedMsg1=[codedMsg1,codedMsg(k,:)];
    end
    codedMsg1=codedMsg1.';
    save codedMsg1 codedMsg1;
    
    % mapping（自带的qammod)
    M=2^BitsPerSymbolQAM;
    modObj   =modem.qammod('M', M, 'SymbolOrder', 'Gray', 'InputType', 'Bit');
    QAMSymbols = modulate(modObj, codedMsg1);
    QAMSymbols = QAMSymbols/rms(QAMSymbols);
    save QAMSymbols_trans QAMSymbols;
    QAMSymbols = reshape(QAMSymbols, length(DataCarrierPositions), SToPcol);
    
    ifftBlock = zeros(FFTSize, SToPcol);
    ifftBlock(DataCarrierPositions, :) = QAMSymbols;
    
end
%% 共轭对称
ifftBlock(FFTSize + 2 - DataCarrierPositions, :) = conj(ifftBlock(DataCarrierPositions,:));
% N-ifft
OFDMSymbols = ifft(ifftBlock);
OFDMSymbols1 = OFDMSymbols(1:length(OFDMPositions),:);
OFDMSymbols1 = [OFDMSymbols1(end-CPLength/2+1:end, :); OFDMSymbols1; OFDMSymbols1(1:CPLength/2, :)];
OFDMSymbols = reshape(OFDMSymbols1, [], 1);


