#!/bin/bash

module purge
module load openmolcas/23.06/gcc12.3.0

export WorkDir=`echo $(mktemp -d)`
mkdir -p $WorkDir

pymolcas input > output

exit 0
