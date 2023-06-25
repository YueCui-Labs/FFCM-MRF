function [Trans_Image] = Hist_match3D(Target_Hist,Image)
% This function is used to do histogram matching
%
% [Trans_Image] = Hist_match3D(Target_Hist,Image)
%
% inputs,
%   Target_Hist : The target histogram 
%   Image : 3D volume data whose histogram will be transformed to the
%   target histogram
% outputs,
%   Trans_Image : The 3D volume data after histogram matching 
%
% Written by  Baochang Zhang
% E-mail：bc.zhang@siat.ac.cn or 1748879448@qq.com

img = Image((Image>=1));%大于1的像素值集合
LengthIMG = numel(img);%求图像中元素个数
Max_intensity = max(img(:));%最大像素点
Test_Hist=zeros(floor(Max_intensity)+1,1);%预设直方图（列向量）
Hist_Cell=cell(floor(Max_intensity)+1,1);%空元胞数组
tic;%秒表开始计时
for i=1:length(Test_Hist)
    Hist_Cell{i}=img(img>=(i-1) & img<i);%小于i像素值的所有像素值（有小数点）集合
    Test_Hist(i)=numel(Hist_Cell{i});%统计每个像素值处的累计数
end
Test_Hist = Test_Hist/ LengthIMG;%归一化（0~1）
t=toc;%显示时间
disp(['calculate the hist of the Image---runtime = ' num2str(t)]);pause(0.1);
% figure
% subplot(2,2,1)
% X=0:1:length(Target_Hist)-1;
% plot(X,Target_Hist,'-k','LineWidth',2);%目标直方图显示（归一化）
% subplot(2,2,2)
% X=0:1:length(Test_Hist)-1;
% plot(X,Test_Hist,'-k','LineWidth',2);%输入图片直方图显示


% calculate  the cumulative histogram计算累计直方图

Target_lens=length(Target_Hist);%目标直方图数组长度
C_target_hist=[];
for i=1:Target_lens
   C_target_hist=[C_target_hist sum(Target_Hist(1:i))];%按序在矩阵后追加目标直方图的像素累计值 
end
% subplot(2,2,3)
% plot(0:Target_lens-1,C_target_hist,'-b','LineWidth',2);%目标直方图的累计图（从0到1）
% the Cumulative distribution of data that will be transformed
Test_lens=length(Test_Hist);
C_test_hist=[];               
for i=1:Test_lens
   C_test_hist=[C_test_hist sum(Test_Hist(1:i))]; 
end
% subplot(2,2,4)
% plot(0:Test_lens-1,C_test_hist,'-r','LineWidth',2);%输入直方图的累计图（从0到1）
% sgtitle('Left: Target-hist,Right:Pre trans-hist');pause(0.1);
% transfer the hist of the image to specific hist
tic;
Map_table=zeros([Target_lens,3]);%第一列：输入图的映射下界；第二列：输入图的映射上界；第三列：输入图的映射点
tick=0;
for i = 1:Target_lens% 
    Map_table(i,1)=tick;% the intensity of the lower limit of the nearest point
    temp=C_test_hist-C_target_hist(i);%输入累计直方图-目标第i累计值
%     [~,index1]=max(temp(find(temp<=0))); 
    [~,index1]=max(temp(temp<=0)); % Lower limit of the nearest point 找到临界值——组映射
    if length(index1)>1%如果有多个符合条件的临界值，选择最大索引（末尾那个）
        index1=index1(end);
    end
    index2=min(index1+1,Test_lens);% Upper limit of the nearest point 找到上界
    new_V=C_test_hist(index1);%下界对应的累计值
    while(Test_Hist(index2)==0)%找到在输入直方图中临界值之上最近的一个非零个像素直方图点
        index2=index2+1;
        if index2==Test_lens
            break;
        end
    end
    A=Test_Hist(index2);% the intensity of the upper limit of the nearest point 上界对应的直方图值
    B=C_target_hist(i)-C_test_hist(index1);
    Ratio=B/A*100;%映射误差
    if Ratio==0
        tick=index2-1-1e-10;%减去一个细微值，作为下一次的下界，也作为本次上界:
    else%在下界到上界中间，找一个能补全映射误差的值作为新的上界（上界下移）
        tick = prctile(Hist_Cell{index2},B/A*100);% 根据区间中的百分比 p 返回小于index2像素值的所有像素点集合 中元素的百分位数。
        new_V=new_V+numel(find(Hist_Cell{index2}>=tick))/LengthIMG;%下界对应的累计分布概率更新：加上
    end
    Map_table(i,2)=tick;%上界
    % Cumulative distribution probability corresponding to transformation
    Map_table(i,3)=new_V; %累计分布概率
end
t=toc;
disp(['transfer the hist of the image to specific hist---runtime = ' num2str(t)]);pause(0.1);
% Complete the transformation table

% Image transformation based on transformation table
tic;
Trans_Image=zeros(size(Image));%映射后图片
for i=1:Target_lens
    inf_V=Map_table(i,1);%本级像素下界
    sup_V=Map_table(i,2);%上界
    Trans_Image(Image>inf_V & Image<=sup_V)=i-1;%满足条件的像素点的值为第i-1灰度级
end
t=toc;
disp(['Image transformation---runtime = ' num2str(t)]);pause(0.1);
img = Trans_Image((Trans_Image>0));
LengthIMG = numel(img);
Max_img = max(img(:));
[~,~,c] = size(Trans_Image);
[N,~] = hist(img(:),0:Max_img); 
Trans_Hist=N'/LengthIMG;

% figure;
% subplot(2,2,1);imshow(imrotate(Image(:,:,fix(c/2)),-90),[]);
% subplot(2,2,2);plot(0:length(Test_Hist)-1,Test_Hist,'-k','LineWidth',2);
% subplot(2,2,3);imshow(imrotate(Trans_Image(:,:,fix(c/2)),-90),[0,400]);
% subplot(2,2,4);plot(0:length(Trans_Hist)-1,Trans_Hist,'-k','LineWidth',2);
% sgtitle('Up: Original Image and hist| Down: Tranfered Image and hist');pause(0.1)

end

