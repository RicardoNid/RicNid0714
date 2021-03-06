%  // ======================================================================
%  //  Jinan University
%  //  @Author: JiZhou CanyangXiong
%  //  @Last Modified time: 2021-03-05      
%  // ======================================================================
function statplot(signal_plot)
rr = real(signal_plot);
ii = imag(signal_plot);
figure
test = [ii rr];
C = hist3(test,[101,101]);
XMax=max(rr);
XMin=min(rr);
YMax=max(ii);
YMin=min(ii);
X = [XMin XMax];
Y = [YMin YMax];
imagesc(X,Y,C);
% imagesc(C);
h = colormap(jet);
h(1,:) = 1;
colormap(h);
colorbar;
set(gca,'YDir','normal');

% title('统计分布星座图','fontweight','bold','fontsize',20);
end