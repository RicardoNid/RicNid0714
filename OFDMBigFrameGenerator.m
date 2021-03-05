%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: 20����֡���һ����֡
%  // ======================================================================

function OFDMFrame = OFDMBigFrameGenerator( OFDMParameters )
OFDMBigFrame = [];
for cir = 1:20
    [OFDMSmallFrame,OFDMSymbols] = OFDMFrameGenerator(OFDMParameters);
    OFDMBigFrame = [OFDMBigFrame,OFDMSmallFrame];
end
OFDMFrame = reshape(OFDMBigFrame, [], 1);

end

