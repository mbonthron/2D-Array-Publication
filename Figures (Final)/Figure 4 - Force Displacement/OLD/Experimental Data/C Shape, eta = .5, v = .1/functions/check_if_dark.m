function [invert] = check_if_dark(image)
    % Check if white or black background
    if size(image,3) == 3
        grayImg = rgb2gray(image);
    else
        grayImg = image;
    end
    avgBrightness = mean(grayImg(:)) / 255;
    if avgBrightness < 0.5
        invert = true; % black for bright images
    else
        invert = false; % white for dark images
    end
end

