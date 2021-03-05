%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: 1¸ö×ÓÖ¡
%  // ======================================================================
function [OFDMSmallFrame, OFDMSymbols] = OFDMFrameGenerator(OFDMParameters)
OFDMSymbols = CreateOFDMSymbols(OFDMParameters);
preamble = CreateOFDMPreamble(OFDMParameters);
OFDMSmallFrame = [preamble; OFDMSymbols];


