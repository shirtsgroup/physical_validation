#!/bin/sh

case $OPERATING_SYSTEM in
    windows*)
        eval "$(${CONDA}/condabin/conda.bat shell.bash hook)";;
    *)
        eval "$(${CONDA}/condabin/conda shell.bash hook)";;
esac