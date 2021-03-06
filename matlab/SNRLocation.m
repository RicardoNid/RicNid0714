%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: ¼ÆËã¸÷×ÓÔØ²¨SNR
%  // ======================================================================
function SNRPerSubcarrier = SNRLocation(recoveredSymbols, transmittedSymbols, OFDMParameters)
SubcarriersNum = length(OFDMParameters.DataCarrierPositions);
recoveredSymbols = reshape(recoveredSymbols,SubcarriersNum,[]);
transmittedSymbols = reshape(transmittedSymbols,SubcarriersNum,[]);
SNR = zeros(SubcarriersNum,1);
SNRdB = zeros(SubcarriersNum,1);
 
for i = 1:SubcarriersNum
    SNR(i) = sum(abs(transmittedSymbols(i,:)).^2) / sum(abs(recoveredSymbols(i,:) - transmittedSymbols(i,:)).^2);
    SNRdB(i) = 10*log10(SNR(i));
end
 
SNRPerSubcarrier = SNR;
save SNRdB SNRdB
figure; plot(SNRdB);
xlabel('subcarrier');
ylabel('SNRdB');

% fid = fopen('E:\MATLABCode\matlabcode8-26\8-26ofdm\SNRPerSubcarrier.txt','wt');
% fprintf(fid,'%g\t',SNR);
% fclose(fid);


end