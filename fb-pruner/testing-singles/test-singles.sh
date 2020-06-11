#!/bin/bash

FB_PRUNER=../build/fb-pruner
MAKE_LAYERED_CLOUD=../data-vis/scripts/make_layered_cloud.sh

for FILE in $(ls)
do
	cd $FILE
	bash ../$MAKE_LAYERED_CLOUD target/ query/
	cd ../
done
