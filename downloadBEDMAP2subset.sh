

# This is a fairly dumb script to download a subset of the BEDMAP2 data set for use in making pretty plots of Antarctica.
# Because we only need a subset of the data, and the compressed data set is fairly large and only used in plotting,
# this installation is optional and not forced by anitaBuildTool.

# You can install the files manually by making sure the files, e.g. bedmap_bed.flt end up in $ANITA_UTIL_INSTALL_DIR/share/anitaCalib/bedmap2_bin/
# This is the recommended way of getting them on the ice... let's hope someone brought them with them!

calibDir=$ANITA_UTIL_INSTALL_DIR/share/anitaCalib
if [[ -z "${ANITA_UTIL_INSTALL_DIR}" ]]; then
    echo "Please set ANITA_UTIL_INSTALL_DIR and try again"
    exit 1;
else

    bedmap2Dir=$calibDir/bedmap2_bin
    if [ -f $bedmap2Dir/bedmap2_bed.flt ] && [ -f $bedmap2Dir/bedmap2_bed.hdr ] && \
	   [ -f $bedmap2Dir/bedmap2_bed.flt ] && [ -f $bedmap2Dir/bedmap2_bed.hdr ] && \
	   [ -f $bedmap2Dir/bedmap2_icemask_grounded_and_shelves.flt ] && [ -f $bedmap2Dir/bedmap2_icemask_grounded_and_shelves.hdr ] && \
	   [ -f $bedmap2Dir/bedmap2_surface.flt ] && [ -f $bedmap2Dir/bedmap2_surface.hdr ] && \
	   [ -f $bedmap2Dir/bedmap2_thickness.flt ] && [ -f $bedmap2Dir/bedmap2_thickness.hdr ]; then
	echo "It seems you've already downloaded the relevent BEDMAP2 data."
	exit 1;
    fi

    echo "*******************************************************************"
    echo "* You are about to download a subset of the BADMAP2 data set.     *"
    echo "* The zip file containing the data set is ~50 MB.                 *"
    echo "* This may take some time depending on your internet speed.       *"
    echo "* Go have a cup of tea in honour of the British Antarctic Survey  *"
    echo "*******************************************************************"

    # Just in case you're doing a complete fresh install and haven't done the make install step
    mkdir -p $calibDir
    echo "Installing to" $calibDir
fi

# We want the binary data
whichData=bedmap2_bin
repo=subsetOfBedmap2Data
releaseNum=1.0
releaseTag=v$releaseNum

curl -L https://github.com/anitaNeutrino/$repo/archive/$releaseTag.tar.gz > $calibDir/$repo-$releaseTag.tar.gz
sleep 3; # Doesn't like trying to unzip straight away

# mv to calibDir to extract
cd $calibDir;

# gzipped tarball of a gzipped tarball...
tar -xvf $repo-$releaseTag.tar.gz
tar -xvf $repo-$releaseNum/$whichData.tar.gz2

# Tidy up
rm $repo-$releaseTag.tar.gz
rm -r $repo-$releaseNum

# return to initial directory
cd -


# Well done
echo "***************************************************"
echo "* Installation complete!                          *"
echo "* See https://www.bas.ac.uk/project/bedmap-2/     *"
echo "* for more info on the data sets.                 *"
echo "***************************************************"
