function [dImg] = dilate(img, mask)
%     [ai,bi]=size(img);
%     [am,bm]=size(mask);
%     deltaa = floor(am/2);
%     deltab = floor(bm/2);
%     out = img;
% 
%     for i = ceil(am/2):ai-deltaa
%         for j = ceil(bm/2):bi-deltab
%             fenster = img(i-deltaa:i+deltaa,j-deltab:j+deltab);
%             fenster = fenster(logical(mask));
%             out(i,j) = max(fenster);
%         end
%     end
%     dImg = out;

    dImg = ordfilt2(img, length(find(mask)), mask);
end