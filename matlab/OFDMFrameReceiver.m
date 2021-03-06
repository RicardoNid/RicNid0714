%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: 接收端DSP
%  // ======================================================================
function  [PerQAMtotal,PerQAMError,QAM_re_sum,pre_code_errors,number_of_error] = OFDMFrameReceiver(recvOFDMFrame, OFDMParameters, cir)
on = OFDMParameters.on;
iterationT = OFDMParameters.iteration;
CPLength = OFDMParameters.CPLength;
BitsPerSymbolQAM = OFDMParameters.BitsPerSymbolQAM;
OFDMSymbolNumber = OFDMParameters.OFDMSymbolNumber;
DataCarrierPositions = OFDMParameters.DataCarrierPositions;
preambleNumber = OFDMParameters.PreambleNumber;
SToPcol = OFDMParameters.SToPcol;
FFTSize = OFDMParameters.FFTSize;
%% FDE
preamble = recvOFDMFrame(1:preambleNumber*(FFTSize+CPLength));
% load rmsOFDMSymbols
% preamble = preamble*rmsOFDMSymbols;
% figure;
% plot(20*log10(abs(fftshift(fft(reshape(recvOFDMFrame,[],1))))));
H = ChannelEstimationByPreamble(preamble, OFDMParameters);
tap = 20;
H = smooth(H,tap);
% figure; 
% plot(20*log10(abs(H)));
recvOFDMSignal = recvOFDMFrame(preambleNumber*(FFTSize+CPLength)+1:end);
[recovered, recvOFDMSignal]=  RecoverOFDMSymbolsWithPilot(recvOFDMSignal, OFDMParameters, H);

% 除对应功率
if on==1
    load('power_alloc.mat');
    for i = 1:SToPcol
        recovered(DataCarrierPositions-2,i) = recovered(DataCarrierPositions-2,i)./sqrt(power_alloc');
    end
else
    recoveredSymbols = reshape(recovered,[],1);
    % scatterplot(recoveredSymbols)
end

if on ==1
    %% bit loading %%
    load('bitAllocSort.mat');
    load('BitAllocSum.mat');
    load('rmsAlloc.mat');
    bitNumber_total =0;
    number_of_error_total = 0;
    demodulated_HD =[];
    for i=1:length(bitAllocSort)
        if bitAllocSort(i)~=0
            carrierPosition = BitAllocSum{i};
            QAM = reshape(recovered(carrierPosition,:),[],1);
            QAM_re = QAM/rms(QAM)*rmsAlloc(i);
            QAM_re_sum{i}= QAM_re;
            % Code properties(channel coding)
            codeRate = 1/2;
            constlen = 7;
            codegen  = [171 133];
            tblen    = 90;
            trellis  = poly2trellis(constlen, codegen);
            % de-mapping
            M=2^bitAllocSort(i);
            if bitAllocSort(i) == 3  %QAM8
                carrierPosition = BitAllocSum{i};
                outPut = reshape(recovered(carrierPosition,:),[],1);
                outPut = outPut';
                QAM8 = [-1-sqrt(3), -1+1i, -1-1i, 1i*(1+sqrt(3)), -1i*(1+sqrt(3)), 1+1i, 1-1i, 1+sqrt(3)] ./ sqrt(3 + sqrt(3));
                [~, index] = min(abs(repmat(outPut, 8, 1) - repmat(transpose(QAM8), 1, length(outPut))));
                temp = de2bi(index - 1, 3, 'left-msb');
                demodulatedMsg_HD = reshape(temp', 1, []);
            else
                modObj   =modem.qammod('M', M, 'SymbolOrder', 'Gray', 'InputType', 'Bit');
                demodObj =modem.qamdemod(modObj);
                set(demodObj, 'DecisionType', 'Hard decision');
                demodulatedMsg_HD = demodulate(demodObj, QAM_re);
                demodulatedMsg_HD=demodulatedMsg_HD';
            end
            % decisiong
            % Set up the demodulator object to perform hard decision demodulation
            demodulated_HD = [demodulated_HD, demodulatedMsg_HD];
        end
    end
    
    save QAM_re_sum QAM_re_sum;
    % de-interleaving
    depth=32;
    len=length(demodulated_HD)/depth;
    codedMsg1=[];
    for k=1:depth
        codedMsg1=[codedMsg1;demodulated_HD(len*(k-1)+1:len*k)];
    end
    demodulatedMsg_HD=codedMsg1(:);
    % viterbi decoder
    % Use the Viterbi decoder in hard decision mode(recvbits)
    load('bits.mat');
    decodedMsg_HD = vitdec(demodulatedMsg_HD, trellis, tblen, 'cont', 'hard');
    decodedMsg_HD=[decodedMsg_HD(tblen+1:end);bits(length(bits)-tblen+1:length(bits))];
    % sendbits
    c=1;
    for i=1:length(bitAllocSort)
        if bitAllocSort(i)~=0
            carrierPosition = BitAllocSum{i};
            bitNumber = OFDMSymbolNumber * length(carrierPosition) * bitAllocSort(i);
            bits_per = randint(bitNumber, 1, 2, OFDMParameters.Seed);
            bitNumber_total =bitNumber_total+bitNumber;
            decodedMsg_HD_per = decodedMsg_HD(c:c+bitNumber-1,1);
            c = bitNumber+c;
            [nErrors_HD, ber_HD] = biterr(decodedMsg_HD_per, bits_per);
            BER_MC=ber_HD;
            number_of_error=nErrors_HD;
            load('QAM_re_sum.mat');
            %               scatterplot(QAM_re_sum{i});
            %               title(['BER-MC = ', num2str(ber_HD)  '----NE= ', num2str(nErrors_HD)]);
            number_of_error_total = number_of_error_total+nErrors_HD;
        end
    end
%     BER_MC =number_of_error_total/bitNumber_total
    % iteration 
    for iter = 1:iterationT
       [PerQAMtotal,PerQAMError,QAM_re_sum,pre_code_errors,decodedMsg_HD,number_of_error_HD]=iteration_alloc(decodedMsg_HD,OFDMParameters,tblen,iter,recovered);
    end
   
%     BER_MC = BER_MC_HD;
%     BER_pre = pre_code_ber;
    number_of_error = number_of_error_HD;  
   %% BER_MC =number_of_error_total/bitNumber_total;
    
    %%
else
    %% cal %%
    recoveredSymbols = recoveredSymbols/rms(recoveredSymbols)*sqrt(10);
    recoveredSymbols_FDE = recoveredSymbols;
    % Code properties(channel coding)
    codeRate = 1/2;
    constlen = 7;
    codegen  = [171 133];
    tblen    = 90;
    trellis  = poly2trellis(constlen, codegen);
    
    % de-mapping
    M=2^BitsPerSymbolQAM;
    modObj   =modem.qammod('M', M, 'SymbolOrder', 'Gray', 'InputType', 'Bit');
    demodObj =modem.qamdemod(modObj);
    
    % decisiong
    % Set up the demodulator object to perform hard decision demodulation
    set(demodObj, 'DecisionType', 'Hard decision');
    demodulatedMsg_HD = demodulate(demodObj, recoveredSymbols);
    demodulatedMsg_HD=demodulatedMsg_HD';
    
    % de-interleaving
    depth=32;
    len=length(demodulatedMsg_HD)/depth;
    codedMsg1=[];
    for k=1:depth
        codedMsg1=[codedMsg1;demodulatedMsg_HD(len*(k-1)+1:len*k)];
    end
    demodulatedMsg_HD=codedMsg1(:);
    
    % viterbi decoder
    % Use the Viterbi decoder in hard decision mode
    decodedMsg_HD = vitdec(demodulatedMsg_HD, trellis, tblen, 'cont', 'hard');
    % Compute the bit error rate
    load('bits.mat');
    bits=bits.';
    decodedMsg_HD=[decodedMsg_HD(tblen+1:end);bits(length(bits)-tblen+1:length(bits))];
    [nErrors_HD, ber_HD] = biterr(decodedMsg_HD, bits);
    
    number_of_error=nErrors_HD;
    BER_MC=ber_HD;
    
    %iteration
    for i = 1:iterationT
        [pre_code_errors,pre_code_ber,BER_MC_HD,decodedMsg_HD,number_of_error_HD]= iteration(decodedMsg_HD,OFDMParameters,recvOFDMSignal,tblen,i,recoveredSymbols_FDE, cir);
    end
    BER_MC = BER_MC_HD;
    BER_pre = pre_code_ber;
    number_of_error = number_of_error_HD;
    
    %因为on=1的时候，需要输出这几个变量，所以on=0的时候，也要添加这几个输出,为了计算on=1时的星座图和不同调制格式的误码率
    QAM_re_sum=0;
    PerQAMError=0;
    PerQAMtotal=0;
    %%
end
