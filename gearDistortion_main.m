clc;
% -----------------
% 2024.10.02 仇韶强
% -----------------
%% 对齿轮图像去畸变并绘制
% 读取齿轮图像
gearImage = imread('12.bmp');

% 相机标定
[distance, cameraParams] = cameraCheck();

% 对齿轮图像进行畸变校正
undistortedGearImage = undistortImage(gearImage, cameraParams);

% 灰度图像处理
% 计算差异
differenceImage = imabsdiff(gearImage, undistortedGearImage);

% 设置阈值，找出发生变化的区域
threshold1 = 20;
binaryMask = differenceImage > threshold1;

markerImage = cat(3, gearImage, gearImage, gearImage); % 将灰度图复制为三通道
markerImage(:,:,1) = gearImage + uint8(binaryMask) * 255; % 红色通道增强
markerImage(:,:,2) = gearImage .* uint8(~binaryMask); % 抑制绿色通道
markerImage(:,:,3) = gearImage .* uint8(~binaryMask); % 抑制蓝色通道

% 显示结果
figure;
imshowpair(gearImage, markerImage, 'montage');
title('原图与标注去畸变区域的图像');

%% 对图像进行裁剪，只保留齿轮部分

% 获取图像的大小
[height, width, ~] = size(undistortedGearImage);

% 定义裁剪区域 (保留中间部分)
cropWidth = round(width * 0.4); % 裁剪区域宽度为图像宽度的 50%
cropHeight = round(height * 0.5); % 裁剪区域高度为图像高度的 50%
centerX = round(width / 2);
centerY = round(height / 2);

% 偏移量：使裁剪区域偏上
verticalOffset = round(cropHeight * 0.25);  % 向上偏移 25%

% 计算裁剪区域的起始坐标
xStart = centerX - round(cropWidth / 2);
yStart = centerY - round(cropHeight / 2) - verticalOffset; % 上移裁剪框

% 裁剪图像
croppedImage = imcrop(undistortedGearImage, [xStart, yStart, cropWidth, cropHeight]);

% 显示裁剪后的图像
figure;
imshow(croppedImage);
title('裁剪后的齿轮图像');

%% 检测齿轮的边缘并测量齿根圆、齿顶圆和孔径
% 把第一节课的作业封装成函数，直接调用
[hole_radius, addendum_radius, dedendum_radius] = process_gear_image(croppedImage,100);

%% 计算比例
% 假设标定板方格的实际边长为 6 毫米
squareSize = 6; % 单位：毫米

% 计算像素到毫米的转换比例
mmPerPixel = distance;

% 显示像素到毫米的比例
fprintf('每像素对应 %.2f 毫米\n', mmPerPixel);

%% 实际齿轮尺寸
radii = [hole_radius, addendum_radius, dedendum_radius];
% 使用 mmPerPixel 将像素转换为毫米
actualRadii = radii * mmPerPixel;

% 输出实际尺寸
fprintf('齿根圆直径: %.2f mm\n', actualRadii(3) * 2);
fprintf('齿顶圆直径: %.2f mm\n', actualRadii(2) * 2);
fprintf('孔径: %.2f mm\n', actualRadii(1) * 2);