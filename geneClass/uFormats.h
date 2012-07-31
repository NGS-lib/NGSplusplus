#ifndef UFORMATS_H_INCLUDED
#define UFORMATS_H_INCLUDED
//Our generic genomic entry format. Everything derives from this.
#include "uGeneException.h"
#include "uFormatBase.h"
#include "uFormatChrom.h"
#include "uFormatExperiment.h"




typedef uGenericNGSExperiment<uGenericNGSChrom<uGenericNGS>,uGenericNGS> GenExp;


/*
Note: We need to seriously refactor most of this code for names.
Standardize our member variables, our parameters names, etc.
Before all the chr, ourchr, pchr get out of hand.

// TODO Most of our Set functions should work by const ref, would optimise stuff
*/
#endif // UFORMATS_H_INCLUDED
