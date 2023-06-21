# Simple cropping script for slices.  You may have to modify to meet your
# needs.

for i in $*
do
    echo cropping ${i}
    convert -crop 841x53+55+17 ${i}_x.jpg ${i}_legend_horizontal.jpg
    convert -rotate 270 ${i}_legend_horizontal.jpg ${i}_legend_vertical.jpg 
    convert -crop 841x841+55+70 ${i}_x.jpg ${i}_x_cropped.jpg
    convert -crop 841x841+55+70 ${i}_y.jpg ${i}_y_cropped.jpg
    convert -crop 841x841+55+70 ${i}_z.jpg ${i}_z_cropped.jpg
    montage -tile 4x1 -geometry "-0-0" -gravity South ${i}_legend.jpg \
            ${i}_x_cropped.jpg ${i}_y_cropped.jpg ${i}_z_cropped.jpg  \
            ${i}_tiled_legend.jpg
    montage -tile 3x1 -geometry "-0-0" -gravity South \
            ${i}_x_cropped.jpg ${i}_y_cropped.jpg ${i}_z_cropped.jpg  \
            ${i}_tiled_nolegend.jpg
done
