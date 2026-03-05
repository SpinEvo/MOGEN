function lumen_xy_loc = get_luminal_xy_loc(radius, center_xy_loc, resolution)

    [xx, yy] = meshgrid(-radius : resolution : radius);
    temp_mask = sqrt(xx.^2 + yy.^2) <= radius;
    lumen_x_loc = xx(temp_mask) + center_xy_loc(1);
    lumen_y_loc = yy(temp_mask) + center_xy_loc(2);
    lumen_xy_loc = [lumen_x_loc, lumen_y_loc];

end