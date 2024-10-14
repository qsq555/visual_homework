function [centers,radii] = detectCircularHoles(image)
    
    % 相机标定
    [distance, cameraParams] = cameraCheck();

    % 读取图像
    rectangleImage = image; 

    % 使用形态学操作去除噪声
    binaryImage = imopen(rectangleImage, strel('disk', 5));

    % 霍夫变换检测圆形
    [centers, radii] = imfindcircles(binaryImage, [25 50], 'ObjectPolarity', 'bright', 'Sensitivity', 0.94);

    % 在原图上绘制检测到的圆
    imshow(rectangleImage);
    hold on;

    % 绘制圆的边界
    viscircles(centers, radii, 'EdgeColor', 'blue', 'LineWidth', 1);
    
    actual_radii = radii * distance;
     % 在圆心位置标记
    for i = 1:length(actual_radii)
        % 用蓝色标注半径
        text(centers(i, 1), centers(i, 2) + actual_radii(i) + 5, sprintf('R=%.2f', actual_radii(i)), 'Color', 'r', 'FontSize', 8);
        
        % 用蓝色标注圆心坐标，偏移圆心位置
        text(centers(i, 1) + 5, centers(i, 2) - 5, sprintf('(%d, %d)', round(centers(i, 1)), round(centers(i, 2))), 'Color', 'blue', 'FontSize', 8);
    end


    %% test
    


    title('检测到的圆孔');
    hold off;
end