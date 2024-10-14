threshold = 130; 
rectangleImage = imread('11.bmp'); % 读取图像

% 相机标定
[distance, cameraParams] = cameraCheck();

% 对冲压件图像进行畸变校正
undistortedrectangleImage = undistortImage(rectangleImage, cameraParams);

% 灰度图像处理
% 计算差异
differenceImage2 = imabsdiff(rectangleImage, undistortedrectangleImage);

% 设置阈值，找出发生变化的区域
threshold1 = 20;
binaryMask = differenceImage2 > threshold1;

markerImage2 = cat(3, rectangleImage, rectangleImage, rectangleImage); % 将灰度图复制为三通道
markerImage2(:,:,1) = rectangleImage + uint8(binaryMask) * 255; % 红色通道增强
markerImage2(:,:,2) = rectangleImage .* uint8(~binaryMask); % 抑制绿色通道
markerImage2(:,:,3) = rectangleImage .* uint8(~binaryMask); % 抑制蓝色通道

% 显示结果
figure;
imshowpair(rectangleImage, markerImage2, 'montage');
title('原图与标注去畸变区域的图像');
    
%% 对图像进行裁剪，只保留冲压件部分

% 获取图像的大小
[height2, width2, ~] = size(undistortedrectangleImage);

% 定义裁剪区域 (保留中间部分)
cropWidth2 = round(width2 * 0.5); % 裁剪区域宽度为图像宽度的 50%
cropHeight2 = round(height2 * 0.5); % 裁剪区域高度为图像高度的 50%
centerX2 = round(width2 / 2);
centerY2 = round(height2 / 2);

% 偏移量：使裁剪区域偏上
verticalOffset2 = round(cropHeight2 * 0.25);  % 向上偏移 25%

% 计算裁剪区域的起始坐标
xStart2 = centerX2 - round(cropWidth2 / 2);
yStart2 = centerY2 - round(cropHeight2 / 2) - verticalOffset2; % 上移裁剪框

% 裁剪图像
croppedImage2 = imcrop(undistortedrectangleImage, [xStart2, yStart2, cropWidth2, cropHeight2]);

% 显示裁剪后的图像
% figure;
% imshow(croppedImage2);
% title('裁剪后的冲压件图像');

%% 二值化
% 查找颜色比偏白颜色更白的像素（像素值大于阈值）
white_mask = croppedImage2 > threshold;

% 创建一个全黑的图像
gray_img = zeros(size(croppedImage2, 1), size(croppedImage2, 2), 'uint8');

% 将符合条件的像素设为白色
gray_img(white_mask) = 255;

% % 显示二值化图片
 figure;
 imshow(gray_img);

%% 找圆
[centers,radii] = detectCircularHoles(gray_img);

%% 拟合直线
grayImg = im2gray(croppedImage2);

    % 边缘检测
    edges = edge(grayImg, 'Canny');
    
    % % 形态学操作：膨胀和腐蚀去掉圆角
    se = strel('line', 5, 0);
    dilatedEdges = imdilate(edges, se);
    cleanedEdges = imerode(dilatedEdges, se);

    % 霍夫变换提取直线
    [H, T, R] = hough(edges);
    % plot(H);
    peaks = houghpeaks(H, 4, 'threshold',ceil(0.15*max(H(:)))); % 获取前4个峰值
    lines = houghlines(cleanedEdges, T, R, peaks);

    % 提取直线参数
    lineParams = zeros(length(lines), 4);
    for k = 1:length(lines)
        % 计算直线方程参数
        x1 = lines(k).point1(1);
        y1 = lines(k).point1(2);
        x2 = lines(k).point2(1);
        y2 = lines(k).point2(2);
        
        % 线的斜率和截距
        if x1 ~= x2
            slope = (y2 - y1) / (x2 - x1);
            intercept = y1 - slope * x1;
            lineParams(k, :) = [slope, intercept, 1, k]; % 线的参数 (斜率, 截距, 类型, 索引)
        else
            % 垂直线
            lineParams(k, :) = [Inf, x1, 0, k]; % x = const
        end
    end

    deleteIndex = [];

    % 判断是否有斜率截距相近的直线，将其合并
    for i = 1:length(lines)
        for j = i+1:length(lines)
            % 计算两条直线的斜率和截距之差
            diff = abs(lineParams(i, 1) - lineParams(j, 1)) + abs(lineParams(i, 2) - lineParams(j, 2));
            if diff < 10
                % 合并两条直线，取平均值
                lineParams(i, 1) = (lineParams(i, 1) + lineParams(j, 1)) / 2;
                lineParams(i, 2) = (lineParams(i, 2) + lineParams(j, 2)) / 2;
                deleteIndex = [deleteIndex, j];
            end
        end
    end

    % 删除合并的直线
    lineParams(deleteIndex, :) = [];
    lines(deleteIndex) = [];

    for i = 1:length(lineParams)
        lineParams(i, 4) = i;
    end

    % 计算平行边的距离
    distances = zeros(length(lineParams)/2, 1);

    % 将直线按斜率排序
    lineParams = sortrows(lineParams, 1);

    % 计算两两直线之间的距离
    for i = 1:2
        slope1 = lineParams(i*2-1, 1);
        % intercept1 = lineParams(i, 2);
        
        slope2 = lineParams(i*2, 1);
        % intercept2 = lineParams(i*2, 2);
        
        % 根据线的类型计算距离
        if slope1 ~= Inf && slope2 ~= Inf
            % 两条非垂直直线
            % d = abs(intercept2 - intercept1) / sqrt(1 + slope1^2);
            x1 = (lines(lineParams(i*2-1,4)).point1(1)+lines(lineParams(i*2-1,4)).point2(1))/2;
            y1 = (lines(lineParams(i*2-1,4)).point1(2)+lines(lineParams(i*2-1,4)).point2(2))/2;
            x2 = (lines(lineParams(i*2,4)).point1(1)+lines(lineParams(i*2,4)).point2(1))/2;
            y2 = (lines(lineParams(i*2,4)).point1(2)+lines(lineParams(i*2,4)).point2(2))/2;
            k = (slope1+slope2)/2;
            A1 = (slope1+slope2)/2;
            B1 = -1;
            C1 = -k*x1+y1;
            % A2 = A1;
            % B2 = -1;
            C2 = -k*x2+y2;
            d=abs(C2-C1)/sqrt(A1^2+B1^2);
        else
            % 垂直线之间的距离
            d = abs(lineParams(i, 2) - lineParams(i+2, 2));
        end
        
        distances(i) = d;
    end

    % 输出距离
    disp('平行边的距离:');
    for i = 1:length(distances)
        fprintf('边 %d 和 边 %d 之间的距离: %.2f mm\n', i, i + 2, distances(i)*distance);
    end


%绘图
figure;
imshow(croppedImage2);
hold on;

% 绘制直边
for k = 1:length(lines)
    xy = [lines(k).point1; lines(k).point2];
    plot(xy(:, 1), xy(:, 2), 'LineWidth', 2, 'Color', 'red');
end

title('提取的矩形直边');
hold off;

%% 孔边距
% 初始化一个存储点到直线距离的数组
pointToLineDistances = zeros(size(centers, 1), length(lines));

% 计算每个点到每条线的距离
for i = 1:size(centers, 1) % 对于每个中心点
    x0 = centers(i, 1);
    y0 = centers(i, 2);
    
    for j = 1:length(lines) % 对于每条直线
        x1 = lines(j).point1(1);
        y1 = lines(j).point1(2);
        x2 = lines(j).point2(1);
        y2 = lines(j).point2(2);
        
        % 计算距离
        distance2 = abs((y2 - y1) * x0 - (x2 - x1) * y0 + x2 * y1 - y2 * x1) / ...
                   sqrt((y2 - y1)^2 + (x2 - x1)^2);
        
        % 存储距离
        pointToLineDistances(i, j) = distance2;
    end
end

% 输出每个点到每条线的距离
for i = 1:size(centers, 1)
    %fprintf('点 %d 到四条线的距离: ', i);
    for j = 1:length(lines)
        %fprintf('线 %d: %.2f ', j, pointToLineDistances(i, j));
    end
    fprintf('\n');
end
% 初始化变量
w1 = 0;
h1 = 0;

% 遍历每个点
for i = 1:size(centers, 1)
    % 获取当前点到四条线的距离
    distances = pointToLineDistances(i, :);
    
    % 找到最小值和第二小值
    sortedDistances = sort(distances); % 对距离进行排序
    minDist = sortedDistances(1); % 最小值
    secondMinDist = sortedDistances(2); % 第二小值
    
    % 更新 w1 和 h1
    w1 = w1 + minDist;
    h1 = h1 + secondMinDist;
end

% 计算平均值
w1 = w1 / size(centers, 1);
h1 = h1 / size(centers, 1);

% 输出结果
%fprintf('最小距离的平均值 w1: %.2f\n', w1);
%fprintf('第二小距离的平均值 h1: %.2f\n', h1);



%% 计算比例
% 计算像素到毫米的转换比例
mmPerPixel = distance;

% 显示像素到毫米的比例
fprintf('每像素对应 %.2f 毫米\n', mmPerPixel);

%% 实际尺寸


% 使用 mmPerPixel 将像素转换为毫米
%actual_rectangle_length = rectangle_length * mmPerPixel;
%actual_rectangle_width = rectangle_width * mmPerPixel;

actual_w1 = w1 * mmPerPixel;
actual_h1 = h1 * mmPerPixel;

% 输出实际尺寸
% fprintf('冲压件长W: %.2f mm\n', actual_rectangle_length);
% fprintf('冲压件宽H: %.2f mm\n', actual_rectangle_width);
fprintf('孔边距w1: %.2f mm\n', actual_w1);
fprintf('孔边距h1: %.2f mm\n', actual_h1);

%% 读取距离
figure;
imshow(croppedImage2);
title('请点击图像上的两个点');

hold on;

% 初始化变量以存储之前的点和距离
previousPoints = [];
previousDistance = [];

while true
    % 使用 ginput 获取用户点击的两个点
    [x, y] = ginput(2); % 选择两个点

    % 计算两点之间的距离
    distance3 = sqrt((x(2) - x(1))^2 + (y(2) - y(1))^2);

    % 在图像上标记点击的点
    plot(x, y, 'ro', 'MarkerSize', 5); % 在点上绘制红色圆点

    % 绘制连接两个点的线
    line([x(1), x(2)], [y(1), y(2)], 'Color', 'b', 'LineWidth', 2); % 绘制蓝色线

    % 在图像上显示距离
    midpoint = [(x(1) + x(2)) / 2, (y(1) + y(2)) / 2]; % 计算中点
    text(midpoint(1), midpoint(2), sprintf('%.2f mm', distance3 * mmPerPixel), 'Color', 'r', 'FontSize', 12, 'HorizontalAlignment', 'center');

    % 存储当前点和距离
    previousPoints = [previousPoints; x, y];
    previousDistance = [previousDistance; distance3];

    % 等待用户按下键
    waitforbuttonpress; % 等待按键输入
    key = get(gcf, 'CurrentKey'); % 获取当前按键

    if strcmp(key, 'escape') % 如果按下 Esc 键
        % 清除图像上的所有内容
        clf;
        imshow(croppedImage2); % 重新显示图像
        title('请点击图像上的两个点'); % 重新设置标题
        hold on; % 继续保持当前图像

    elseif strcmp(key, 'return') % 如果按下 Enter 键
        % 保留上次点的两个点和显示图像
        for i = 1:length(previousDistance)
            % 在图像上标记之前的点
            plot(previousPoints(i, 1), previousPoints(i, 2), 'ro', 'MarkerSize', 5); % 在点上绘制红色圆点

            % 绘制连接之前的点的线
            if mod(i, 2) == 0 % 每两个点一组
                line([previousPoints(i-1, 1), previousPoints(i, 1)], ...
                     [previousPoints(i-1, 2), previousPoints(i, 2)], 'Color', 'b', 'LineWidth', 2); % 绘制蓝色线
                % 在图像上显示距离
                midpoint = [(previousPoints(i-1, 1) + previousPoints(i, 1)) / 2, ...
                             (previousPoints(i-1, 2) + previousPoints(i, 2)) / 2]; % 计算中点
                text(midpoint(1), midpoint(2), sprintf('%.2f mm', previousDistance(i/2) * mmPerPixel), ...
                     'Color', 'r', 'FontSize', 12, 'HorizontalAlignment', 'center');
            end
        end
    end
end

hold off;