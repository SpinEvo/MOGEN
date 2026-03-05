function lumen_pixel_indices = get_luminal_pixel_indices(radius_in_pixel, img_size, center_pixel)

    [xx, yy] = meshgrid(-radius_in_pixel : radius_in_pixel);
    local_center_pixel = [radius_in_pixel + 1, radius_in_pixel + 1];
    temp_mask = sqrt(xx.^2 + yy.^2) <= radius_in_pixel;
    [lumen_pixel_row, lumen_pixel_column] = find(temp_mask);

    lumen_pixel_row = lumen_pixel_row - local_center_pixel(1) + center_pixel(1);
    lumen_pixel_column = lumen_pixel_column - local_center_pixel(2) + center_pixel(2);

    lumen_pixel_indices = sub2ind(img_size, lumen_pixel_row, lumen_pixel_column);

end