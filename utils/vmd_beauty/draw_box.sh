#! /bin/bash

box=$1

minx=$(awk 'NR==1{print $1}' $box)
miny=$(awk 'NR==2{print $1}' $box)
minz=$(awk 'NR==3{print $1}' $box)
maxx=$(awk 'NR==1{print $2}' $box)
maxy=$(awk 'NR==2{print $2}' $box)
maxz=$(awk 'NR==3{print $2}' $box)

echo "draw materials off"
echo "     draw line \"$minx $miny $minz\" \"$maxx $miny $minz\""
echo "     draw line \"$minx $miny $minz\" \"$minx $maxy $minz\""
echo "     draw line \"$minx $miny $minz\" \"$minx $miny $maxz\""
echo "     draw line \"$maxx $miny $minz\" \"$maxx $maxy $minz\""
echo "     draw line \"$maxx $miny $minz\" \"$maxx $miny $maxz\""
echo "     draw line \"$minx $maxy $minz\" \"$maxx $maxy $minz\""
echo "     draw line \"$minx $maxy $minz\" \"$minx $maxy $maxz\""
echo "     draw line \"$minx $miny $maxz\" \"$maxx $miny $maxz\""
echo "     draw line \"$minx $miny $maxz\" \"$minx $maxy $maxz\""
echo "     draw line \"$maxx $maxy $maxz\" \"$maxx $maxy $minz\""
echo "     draw line \"$maxx $maxy $maxz\" \"$minx $maxy $maxz\""
echo "     draw line \"$maxx $maxy $maxz\" \"$maxx $miny $maxz\""
