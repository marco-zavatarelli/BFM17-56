#!/bin/bash

#-include-directory=./class -include-directory=./chaps --include-directory=./img 

NAME=bfmProfilingReport
OUTDIR=./out
export TEXINPUTS=./sty:./img:./tex:./bib:

#clean the output files generated before
CLEAN=
if [[ $1 && $1 == 'clean' ]]; then CLEAN=1; fi
if [ $CLEAN ]; then
    rm -rf $OUTDIR/*
fi


ls $OUTDIR/*| grep -v $NAME.vrs | xargs rm -rf

echo "------------LATEX----------------"
latex -output-directory=$OUTDIR $NAME.tex
if [ $? != 0 ]; then echo "$?"; echo "ERROR STEP 1"; exit; fi
echo "------------BIBTEX--------------"
bibtex ./out/$NAME
latex -output-directory=$OUTDIR $NAME.tex
if [ $? != 0 ]; then echo "$?"; echo "ERROR STEP 2"; exit; fi
echo ""
echo "------------RERUN--------------"
echo ""
latex -output-directory=$OUTDIR $NAME.tex
if [ $? != 0 ]; then echo "$?"; echo "ERROR STEP 3"; exit; fi
echo "------------PDF-----------------"
dvips -t a4 -Ppdf -o $OUTDIR/$NAME.ps $OUTDIR/$NAME.dvi
ps2pdf $OUTDIR/$NAME.ps $OUTDIR/$NAME.pdf
