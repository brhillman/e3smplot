# Georeferencing raster images
For the blue marble plots, we need to georeference the NASA blue marble PNG images (or similar). I used GDAL for this (`conda install -c conda-forge gdal`). I kind of just got lucky by not needing to add the coordinate information and getting the right result by running `gdal_translate` with just the input/output file/format arguments. E.g.,
```
gdal_translate -of netcdf -if png img00082.png img00082.nc
```

The resulting files can then be plotted using my `plot_image_gdal.py` code. I am sure there is an intelligent way of adding the coordinate information to this command, and should look into that in the future.
