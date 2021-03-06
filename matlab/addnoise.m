%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  //  @description: º”‘Î…˘      
%  // ======================================================================
function [out,noise]=addnoise(input,attn)
%%%%this function is used to generate AWGN nosie
iout=randn(1,length(input)).*attn;
noise=iout';
out=input+noise;
end