# GeneStructurePlot

This script is design to draw gene **exon-intron structure** for a few of genes based on GFF3 format. And output high quality vector image in SVG format

> please see examples/sample.svg for an example

## Requirements:

### Perl Modules:

   * [FuhaoPerl5Lib](https://github.com/lufuhao/FuhaoPerl5Lib) (MiscKit and GffKit)

   * Getopt::Long

   * SVG


## Options

> perl $0 -i input.gff3 -c in.config -o out.svg [Options]
> 
> Version: v20181128
> 
>    --help|-h
> 
>        Print this help/usage;
> 
>    --input|-i  <in.GFF3>
> 
>        Input GFF3 file
> 
>    --config|-c <in.config>
> 
>        Configure file, see examples/sample.config
> 
>    --output|-o <out.svg>
> 
>        Output SVG image in vector format
> 
>    --debug
> 
>        Try in to locate code problems in dubug mode
> 
>    --verbose
> 
>        Detailed output for trouble-shooting;
> 
>    --version|v!
> 
>        Print current SCRIPT version;

## Example:

> perl $0 -i examples/sample.gff3 -c examples/sample.config -o examples/sample.svg

## Author:

>    Fu-Hao Lu
> 
>    Post-Doctoral Scientist in Micheal Bevan laboratory
> 
>    Cell and Developmental Department, John Innes Centre
> 
>    Norwich NR4 7UH, United Kingdom
> 
>    E-mail: Fu-Hao.Lu\@jic.ac.uk
