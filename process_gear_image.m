function [hole_radius, addendum_radius, dedendum_radius] = process_gear_image(img, threshold)
    % 如果图像是多通道，将其转换为灰度图像
    if size(img, 3) > 1
        img = rgb2gray(img);
    end
    
    % 查找颜色比偏白颜色更白的像素（像素值大于阈值）
    white_mask = img > threshold;

    % 创建一个全黑的图像
    gray_img = zeros(size(img, 1), size(img, 2), 'uint8');
    
    % 将符合条件的像素设为白色
    gray_img(white_mask) = 255;

    % 查找齿轮的边缘轮廓
    boundaries = bwboundaries(gray_img);

    boundary_sizes = [];
    
    % 创建一个数组来存储符合条件的轮廓
    valid_boundaries = {};

    % 遍历每个轮廓，过滤掉点数少于100的轮廓
    for k = 1:length(boundaries)
        boundary = boundaries{k};
        if size(boundary, 1) >= 100
            valid_boundaries{end+1} = boundary;
            boundary_sizes(end+1) = size(boundary, 1);
        end
    end

    % 去掉最大的轮廓
    [~, max_idx] = max(boundary_sizes);
    valid_boundaries(max_idx) = [];
    boundary_sizes(max_idx) = [];

    % 使用最大轮廓点数的轮廓
    [~, max_idx] = max(boundary_sizes);
    selected_boundary = valid_boundaries{max_idx};
    
    % 拟合圆
    [center, ~] = fit_circle_to_points(selected_boundary);
    center_real = center;
    
    % 创建一个数组来存储所有轮廓点到圆心的距离
    all_distances = [];

    % 遍历每个轮廓
    for k = 1:length(valid_boundaries)
        boundary = valid_boundaries{k};
    
        % 计算当前轮廓所有点到圆心的距离
        distances = sqrt((boundary(:,1) - center_real(2)).^2 + (boundary(:,2) - center_real(1)).^2);
    
        % 将这些距离添加到总数组中
        all_distances = [all_distances; distances];
    end

    % 对所有距离进行排序
    sorted_distances = sort(all_distances);

    % 获取最小和最大距离
    min_distance = sorted_distances(1);
    max_distance = sorted_distances(end);

    % 绘制以最小和最大距离为半径的圆
    viscircles(center_real(1, :), min_distance, 'EdgeColor', 'g', 'LineWidth', 2); % 绿色为最小圆
    viscircles(center_real(1, :), max_distance, 'EdgeColor', 'b', 'LineWidth', 2); % 蓝色为最大圆
    %fprintf('齿轮孔半径：%f\n', min_distance);
    %fprintf('齿顶圆半径：%f\n', max_distance);

    % 找出相邻点之间差值大于2的距离
    significant_distances = [];

    for i = 1:length(sorted_distances)-1
        % 如果相邻两个距离的差值大于2
        if sorted_distances(i+1) - sorted_distances(i) > 2
            % 记录较大的那个距离
            significant_distances = [significant_distances; sorted_distances(i+1)];
        end
    end

    % 如果找到符合条件的距离，绘制对应的圆
    if ~isempty(significant_distances)
        for i = 1:length(significant_distances)
            viscircles(center_real(1, :), significant_distances(i) , 'EdgeColor', 'b', 'LineWidth', 2); % 蓝色圆
            %fprintf('齿根圆半径：%f\n', significant_distances(i) );
        end
   
    else
        disp('没有找到符合条件的距离');
    end
    hole_radius = min_distance;
    addendum_radius = max_distance;
    dedendum_radius = significant_distances;  
end


function [center, radius] = fit_circle_to_points(points)
    % 适用于点集的圆拟合
    % 使用最小二乘法拟合圆
    % 输入: points - N x 2 矩阵，其中每一行是一个点 (x, y)
    
    % 转置点集
    x = points(:,1);
    y = points(:,2);
    
    % 使用最小二乘法拟合圆
    A = [x, y, ones(size(x))];
    B = -(x.^2 + y.^2);
    
    % 计算圆心和半径
    coeff = A \ B;
    x0 = -coeff(1)/2;
    y0 = -coeff(2)/2;
    radius = sqrt((coeff(1)^2 + coeff(2)^2)/4 - coeff(3));
    
    center = [y0, x0];
end







