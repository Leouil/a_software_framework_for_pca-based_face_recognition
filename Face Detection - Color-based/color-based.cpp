function BinImg = GetFaceBin(rgb)
Ycbcr = rgb2ycbcr(rgb);

fThreshold = 0.22;

[M,N,D]=size(Ycbcr);
FaceProbImg = zeros(M,N,1);
BinImg = uint8(zeros(M,N,1));
Mean = [117.4316 148.5599]';
C = [97.0946 24.4700;
     24.4700 141.9966];
cbcr = zeros(2,1);
for i=1:M
for j=1:N
cbcr(1) = Ycbcr(i,j,2);
cbcr(2) = Ycbcr(i,j,3);
FaceProbImg(i,j)=exp(-0.5*(cbcr-Mean)'*inv(C)*(cbcr-Mean));
                     if FaceProbImg(i,j)>fThreshold
                     BinImg(i,j) = 1;
                     end
                     end
                     end
                     se=strel('disk',3');
                              BinImg = imopen(BinImg,se);
                              imdilate(BinImg,se);
                              % subplot(122);imshow(BinImg*255);title('');
                              CC = bwconncomp(BinImg);
                              numPixels = cellfun(@numel,CC.PixelIdxList);。
                              [biggest,idx] = max(numPixels);
                              for i=1:CC.NumObjects
                              if i~=idx
                              BinImg(CC.PixelIdxList{i}) = 0;
                              end
                              end
                              % figure(2);imshow(uint8(BinImg*255));
                              end
