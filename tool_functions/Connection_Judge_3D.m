function [Connect_Asked, C_number,list] = Connection_Judge_3D(Binary_Img, maxconnectnum,connectIndex,minC_number,mode)
%
% This function is used to get the connection region in the 3D binarydata, the
% connection regions have been resorted in descending order
%
% [Connect_Asked, C_number,list] = Connection_Judge_3D(Binary_Img, maxconnectnum,connectIndex,minC_number,mode)
%
% inputs,
%   Binary_Img : The input image volume (binary Image)
%   maxconnectnum :A constant, which determines the number of biggest
%                  connection region.
%   connectIndex: A list, which determines some special index of
%                 connection region
%   minC_number: A constant, which requires the voxels of connection
%                region must bigger than the given contanst.
%   mode: 
%         if mode:1 =>  Get the top 'maxconnectnum' maximum connection region
%         if mode:2 =>  Get the connection regionof the specified list(connectIndex)
%         if mode:3 =>  Get the connection regions whose the voxels are bigger 
%                       than the given contanst minC_number.
% outputs,
%   Connect_Asked : A volume containing only the required connection
%                   regions
%   C_number : A count of all connection regions
%   list: A list that recodes the voxels in each connection region 
%
% Example,
%	[Connections, ~,~] = Connection_Judge_3D(Binary_Img, 5,[1 3,5 7 9],300,1)
%   [Connections, ~,~] = Connection_Judge_3D(Binary_Img, 5,[1 3,5 7 9],300,2)
%   [Connections, ~,~] = Connection_Judge_3D(Binary_Img, 5,[1 3,5 7 9],300,3)
%


    Connect_Asked =zeros(size(Binary_Img));
    % find components
    CC = bwconncomp(Binary_Img);
    % component number
    C_number = CC.NumObjects;
    list = [];
    for k = 1:C_number
        cPixel_N = length( CC.PixelIdxList{k} );
        list = [list;cPixel_N,k];
    end
    list = sortrows(list,-1);
    if mode==1
         for k=1:maxconnectnum
             index = list(k,2);
             Connect_Asked(CC.PixelIdxList{index}) = 1;
         end
    end
    if mode ==2
        for k = 1:length(connectIndex)
            index = list(connectIndex(k),2);
            Connect_Asked(CC.PixelIdxList{index}) = 1;
        end
    end
    if mode ==3
        k=1;
        while list(k,1)>minC_number
            index = list(k,2);
            Connect_Asked(CC.PixelIdxList{index}) = 1;
            k=k+1;
        end
    end
end
