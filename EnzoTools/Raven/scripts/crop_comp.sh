# Simple cropping script for slices.  You may have to modify to meet your
# needs.

for i in `seq -w 0 199`
do
    echo cropping ${i}
    convert -quality 100 -crop 837x837+55+70 max_0${i}_Density_x.jpg max_0${i}_x_cropped.jpg
    convert -quality 100 -crop 837x837+55+70 sum_0${i}_Density_x.jpg sum_0${i}_x_cropped.jpg
    montage -quality 100 -geometry "+0+0" max_0${i}_x_cropped.jpg sum_0${i}_x_cropped.jpg comb_0${i}_x.jpg
done
