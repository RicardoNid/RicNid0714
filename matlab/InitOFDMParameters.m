%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: 参数设置     
%  // ======================================================================
function OFDMParameters = InitOFDMParameters()
% on=0算比特分配，on=1比特加载
OFDMParameters.on = 0;
OFDMParameters.CPLength = 20; 
OFDMParameters.PreambleNumber = 2;
OFDMParameters.FFTSize = 512;
OFDMParameters.OFDMSymbolNumber = 8; %OFDM符号数量
OFDMParameters.DataSymbolAmplitude = 1; %ASK信号数据幅度
OFDMParameters.DCCarrierPosition = OFDMParameters.FFTSize/2+1;
OFDMParameters.iteration = 5;
%复数
% OFDMParameters.DataCarrierPositions = [1:93,108:200];
%实数
OFDMParameters.DataCarrierPositions = [3:226];
OFDMParameters.OFDMPositions = sort([1 OFDMParameters.DataCarrierPositions  OFDMParameters.FFTSize/2+1 OFDMParameters.FFTSize + 2 - OFDMParameters.DataCarrierPositions ]);
OFDMParameters.BitsPerSymbolQAM = 4;
OFDMParameters.Seed = 21;
% OFDMParameters.Seed = 10;

OFDMParameters.PreambleCarrierPositions = [-OFDMParameters.FFTSize/2+1:-1] + OFDMParameters.DCCarrierPosition;
OFDMParameters.PreambleBitsPerSymbolQAM = 4;
OFDMParameters.PreambleSeed = 20;
OFDMParameters.PreambleSymbolAmplitude = 1;

% OFDMSymbolNumber*2  S/P的列数 编码效率1/2,k=1,n=2
OFDMParameters.codeRate = 1/2;
OFDMParameters.bitNumber = OFDMParameters.OFDMSymbolNumber * length(OFDMParameters.DataCarrierPositions) * OFDMParameters.BitsPerSymbolQAM;
OFDMParameters.SToPcol = ((OFDMParameters.bitNumber*(1/OFDMParameters.codeRate))/OFDMParameters.BitsPerSymbolQAM)/length(OFDMParameters.DataCarrierPositions);