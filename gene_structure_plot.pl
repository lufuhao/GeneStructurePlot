#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use FuhaoPerl5Lib::MiscKit qw/ReadConfig MrnaSort/;
use FuhaoPerl5Lib::GffKit qw/ReadGff3 GffAddUTR/;
use List::Util qw/min max/;
use SVG;
use constant USAGE=><<EOH;

SYNOPSIS:

perl $0 -i input.gff3 -o out.svg [Options]

Version: v20181128

Requirements:
    Modules: FuhaoPerl5Lib::MiscKit;
             FuhaoPerl5Lib::GffKit;
             Getopt::Long;
             SVG;

Descriptions:
    Draw exon-intron structure based on GFF3
    Output as SVG format

Options:
    --help|-h
        Print this help/usage;
    --input|-i  <in.GFF3>
        Input GFF3 file
    --config|-c <in.config>
        Configure file, see examples/sample.config
    --output|-o <out.svg>
        Output SVG image in vector format
    --debug
        Try in to locate code problems in dubug mode
    --verbose
        Detailed output for trouble-shooting;
    --version|v!
        Print current SCRIPT version;

Example:
    perl $0 -i examples/sample.gff3 -c examples/sample.config -o examples/sample.svg

Author:
    Fu-Hao Lu
    Post-Doctoral Scientist in Micheal Bevan laboratory
    Cell and Developmental Department, John Innes Centre
    Norwich NR4 7UH, United Kingdom
    E-mail: Fu-Hao.Lu\@jic.ac.uk

EOH
###HELP ends#########################################################
die USAGE unless @ARGV;



###Receving parameter################################################
my ($help, $verbose, $debug, $version);
my ($fil_gff3, $output, $file_config);

GetOptions(
	"help|h!" => \$help,
	"input|i=s" => \$fil_gff3,
	"config|c=s" => \$file_config,
	"output|o:s" => \$output,
	"debug!" => \$debug,
	"verbose!" => \$verbose,
	"version|v!" => \$version) or die USAGE;
($help or $version) and die USAGE;



### Defaults ########################################################
$debug=0 unless (defined $debug);
$verbose=0 unless (defined $verbose);



### input and output ################################################



### Main ############################################################
my $config=ReadConfig($file_config);
my $test_utr;
my ($success, $referenceids, $gene, $gene2mrna, $mrnas, $exons, $cds, $utr)=ReadGff3($fil_gff3);
if (scalar(keys %{$utr})==0) {
	($test_utr, $utr)=GffAddUTR($mrnas, $exons, $cds);
	if (scalar(keys %{$utr})==0) {
		$test_utr=0;
	}
}
else {
	$test_utr=1;
}
unless ($success) {
	die "Error: read GFF3 error\n";
}
my $plot_x_width=${$config}{'main'}{'plot_width'}-${$config}{'main'}{'plot_margin_left_size'}-${$config}{'main'}{'plot_margin_right_size'};


my $num_tracks=0;
my %num_mRNAs=();
my $total_num_mrna=0;
my %chrom_start=();
my %chrom_end=();
my %track_start=();
my $max_x=0;

### setting plot height

my $ind_y_start=${$config}{'main'}{'plot_margin_top_size'};
foreach my $ind_ref (sort keys %{$referenceids}) {
	my @draw_chr_length=();
	$track_start{$ind_ref}=$ind_y_start;
	unless (${$config}{'chrom'}{'chrom_cut_end'}) {
		unless (exists ${$config}{'chrom'}{$ind_ref} and ${$config}{'chrom'}{$ind_ref}=~/^\d+$/) {
			die "Error: specified [chrom]chrom_cut_end but not provide chromosome length\n";
		}
		push (@draw_chr_length, 1);
		push (@draw_chr_length, ${$config}{'chrom'}{$ind_ref});
	}
	$num_tracks++;
	my $largest_num_mRNA=0;
	foreach my $indpos (sort {$a<=>$b} keys %{${$referenceids}{$ind_ref}}) {
		foreach my $ind_gene (sort keys %{${$referenceids}{$ind_ref}{$indpos}}) {
			if (${$config}{'chrom'}{'chrom_cut_end'}) {
				unless (exists ${$gene}{$ind_gene} and exists ${$gene}{$ind_gene}{'end'} and exists ${$gene}{$ind_gene}{'start'} and ${$gene}{$ind_gene}{'end'}=~/^\d+$/ and ${$gene}{$ind_gene}{'start'}=~/^\d+$/ and ${$gene}{$ind_gene}{'start'}<${$gene}{$ind_gene}{'end'}) {
					die "Error: invalid gene coordinates: $ind_gene\n";
				}
				push (@draw_chr_length, ${$gene}{$ind_gene}{'start'});
				push (@draw_chr_length, ${$gene}{$ind_gene}{'end'});
			}
			if (scalar(keys %{${$gene2mrna}{$ind_gene}})>$largest_num_mRNA) {
				$largest_num_mRNA=scalar(keys %{${$gene2mrna}{$ind_gene}});
			}
		}
	}
	@draw_chr_length=sort {$a<=>$b} @draw_chr_length;
	$chrom_start{$ind_ref}= (($draw_chr_length[0]-${$config}{'chrom'}{'chrom_flanking_length'})>0) ? ($draw_chr_length[0]-${$config}{'chrom'}{'chrom_flanking_length'}) : 1;
	if (exists ${$config}{'chrom'} and exists ${$config}{'chrom'}{lc($ind_ref)}) {
		$chrom_end{$ind_ref}=(($draw_chr_length[-1]+${$config}{'chrom'}{'chrom_flanking_length'})<=${$config}{'chrom'}{lc($ind_ref)}) ? ($draw_chr_length[-1]+${$config}{'chrom'}{'chrom_flanking_length'}) : ${$config}{'chrom'}{lc($ind_ref)};
	}
	else {
		$chrom_end{$ind_ref}=$draw_chr_length[-1]+${$config}{'chrom'}{'chrom_flanking_length'};
	}
	if ($verbose) {
		print "Test: chr: start - end $ind_ref: ", $chrom_start{$ind_ref}, " - ", $chrom_end{$ind_ref}, "\n";
	}
	my $idv_chr_len=$chrom_end{$ind_ref}-$chrom_start{$ind_ref}+1;
	if ($idv_chr_len>0) {
		if ($idv_chr_len>$max_x) {
			$max_x=$idv_chr_len;
		}
	}
	$num_mRNAs{$ind_ref}=$largest_num_mRNA;
	$total_num_mrna+=$largest_num_mRNA;
	$ind_y_start=$ind_y_start+${$config}{'gene'}{'gene_track_height'}+${$config}{'chrom'}{'chrom_track_space'}+$largest_num_mRNA * (${$config}{'exon'}{'exon_track_space'} + ${$config}{'exon'}{'exon_track_height'});
	if (exists ${$config}{'intron'} and exists ${$config}{'intron'}{'intron_joiner_draw'} and ${$config}{'intron'}{'intron_joiner_draw'}=~/^true$/i) {
		$ind_y_start=$ind_y_start+$largest_num_mRNA * ${$config}{'intron'}{'intron_joiner_height'};
	}
}
my $total_y=${$config}{'main'}{'plot_margin_top_size'}+${$config}{'main'}{'plot_margin_bottom_size'}+($num_tracks-1)*${$config}{'chrom'}{'chrom_track_space'} + $num_tracks * ${$config}{'gene'}{'gene_track_height'}+ $total_num_mrna*(${$config}{'exon'}{'exon_track_space'} + ${$config}{'exon'}{'exon_track_height'});

if (exists ${$config}{'intron'} and exists ${$config}{'intron'}{'intron_joiner_draw'} and ${$config}{'intron'}{'intron_joiner_draw'}=~/^true$/i) {
	$total_y=$total_y+$total_num_mrna*${$config}{'intron'}{'intron_joiner_height'};
}

print "INFO: plot width : ", ${$config}{'main'}{'plot_width'}, " height: $total_y\n";
my $x_factor=$plot_x_width/$max_x;
if ($verbose) {
	print "Test: max chromosome length: $max_x\n";
}


### start plot
my $vectorout=SVG->new(width=>${$config}{'main'}{'plot_width'}, height=>$total_y);

if (${$config}{'main'}{'plot_background_draw'}=~/true$/i) {
	$vectorout->rectangle(x => 0, 
				y => 0, 
				width  	=> ${$config}{'main'}{'plot_width'}, 
				height => $total_y,
				id=> "background",
				style => {'fill' => ${$config}{'main'}{'plot_fill_color'},
				'fill-opacity'   => ${$config}{'main'}{'plot_fill_opaque'},
				'stroke'         => ${$config}{'main'}{'plot_stroke_color'},
				'stroke-width'   => ${$config}{'main'}{'plot_stroke_size'},
				'stroke-opacity' => ${$config}{'main'}{'plot_stroke_opaque'},
						},
				);
}
my ($xx1,$xx2,$yy1,$yy2);
foreach my $ind_ref (sort keys %{$referenceids}) {
	my $x_start=${$config}{'main'}{'plot_margin_left_size'};
	my $y_start=$track_start{$ind_ref};
	$xx1=$x_start;$xx2=$x_start + ($chrom_end{$ind_ref}-$chrom_start{$ind_ref}+1)*$x_factor;
	$yy1=$y_start;$yy2=$y_start;
	if ($verbose) {
		print "SEQ : $ind_ref;  COORDINATES : X1 $xx1 X2 $xx2 Y1 $yy1 Y2 $yy2\n";
	}
	$vectorout->line ( id=> "$ind_ref",
						x1 => $xx1,
						y1 => $yy1+${$config}{'gene'}{'gene_track_height'}/2, 
						x2 => $xx2, 
						y2 => $yy2+${$config}{'gene'}{'gene_track_height'}/2, 
						stroke => ${$config}{'chrom'}{'chrom_stroke_color'}, 
						"stroke-width" => ${$config}{'chrom'}{'chrom_stroke_size'}
					);
	$vectorout->text(	x => $xx1, 
						y => $yy1, 
						width => ${$config}{'main'}{'plot_font_size'}, 
						height => ${$config}{'main'}{'plot_font_size'}, 
						"font-family"=>${$config}{'main'}{'plot_font_family'}, 
						"text-anchor"=>"start",
						"font-size"=>${$config}{'main'}{'plot_font_size'}, 
						"-cdata" => "$ind_ref");
	foreach my $indpos (sort {$a<=>$b} keys %{${$referenceids}{$ind_ref}}) {
		foreach my $ind_gene (sort keys %{${$referenceids}{$ind_ref}{$indpos}}) {
			print "  Gene: $ind_gene\n";
			$vectorout->rectangle(x => $x_start+(${$gene}{$ind_gene}{'start'}-$chrom_start{$ind_ref}+1)*$x_factor, 
						y => $y_start, 
						width  	=> (${$gene}{$ind_gene}{'end'}-${$gene}{$ind_gene}{'start'}+1)*$x_factor, 
						height => ${$config}{'gene'}{'gene_track_height'},
						id=> "$ind_gene",
						style => {'fill' => ${$config}{'gene'}{'gene_fill_color'},
								'fill-opacity'   =>  ${$config}{'gene'}{'gene_fill_opaque'},
								'stroke'         => ${$config}{'gene'}{'gene_stroke_color'},
								'stroke-width'   =>  ${$config}{'gene'}{'gene_stroke_opaque'},
								'stroke-opacity' =>  1,
							},
						);
			$vectorout->text(	x => $x_start+(${$gene}{$ind_gene}{'start'}+${$gene}{$ind_gene}{'end'}-2*$chrom_start{$ind_ref}+2)*$x_factor/2, 
						y => $y_start+${$config}{'gene'}{'gene_track_height'}+${$config}{'main'}{'plot_font_size'}, 
						width => ${$config}{'main'}{'plot_font_size'}, 
						height => ${$config}{'main'}{'plot_font_size'}, 
						"font-family"=>${$config}{'main'}{'plot_font_family'}, 
						"text-anchor"=>"middle",
						"font-size"=>${$config}{'main'}{'plot_font_size'}, 
						"-cdata" => "$ind_gene");
			if (${$config}{'gene'}{'gene_orientation_draw'}=~/^true$/i) {
				my $gene_x1=$x_start+(${$gene}{$ind_gene}{'start'}-$chrom_start{$ind_ref}+1)*$x_factor;
				my $gene_x2=$x_start+(${$gene}{$ind_gene}{'end'}-$chrom_start{$ind_ref}+1)*$x_factor;
				my @x_coords=&GetMarkerCoords($gene_x1, $gene_x2, ${$config}{'gene'}{'gene_orientation_width'}, ${$config}{'gene'}{'gene_orientation_space'});
				foreach my $gene_idv_x (@x_coords) {
						my $xv = [$gene_idv_x,$gene_idv_x,$gene_idv_x];
						if (${$gene}{$ind_gene}{'strand'} eq '+') {
							$xv = [$gene_idv_x-${$config}{'gene'}{'gene_orientation_width'}/2,$gene_idv_x+${$config}{'gene'}{'gene_orientation_width'}/2,$gene_idv_x-${$config}{'gene'}{'gene_orientation_width'}/2];
						}
						elsif (${$gene}{$ind_gene}{'strand'} eq '-') {
							$xv = [$gene_idv_x+${$config}{'gene'}{'gene_orientation_width'}/2,$gene_idv_x-${$config}{'gene'}{'gene_orientation_width'}/2,$gene_idv_x+${$config}{'gene'}{'gene_orientation_width'}/2];
						}
						my $yv = [$y_start+${$config}{'gene'}{'gene_track_height'}/2-${$config}{'gene'}{'gene_orientation_height'}/2, $y_start+${$config}{'gene'}{'gene_track_height'}/2, $y_start+${$config}{'gene'}{'gene_track_height'}/2+${$config}{'gene'}{'gene_orientation_height'}/2];
						my $points = $vectorout->get_path(
						    x       => $xv,
						    y       => $yv,
						    -type   => 'polyline',
						    -closed => 'false' #specify that the polyline is closed.
						);
						$vectorout->polyline (
							%$points,
							id    =>"gene_orientation-$gene_idv_x",
							style => {
								'fill-opacity' => 0,
								'stroke'       => ${$config}{'gene'}{'gene_orientation_color'}, 
								"stroke-width" => ${$config}{'gene'}{'gene_orientation_size'}
							}
						);
				}
			}
			
			my $y_mRNAstart=$track_start{$ind_ref}+${$config}{'gene'}{'gene_track_height'}+${$config}{'gene'}{'gene_track_height'}+${$config}{'exon'}{'exon_track_space'};
			my $num_mrna=0;
			my @mrna_ids=MrnaSort(keys %{${$gene2mrna}{$ind_gene}});
			foreach my $ind_mrna (@mrna_ids) {
				if ($verbose) {
					print "    mRNA ID: $ind_mrna\n";
				}
				$num_mrna++;
				my $this_mrna_y=$y_mRNAstart+($num_mrna-1)*(${$config}{'exon'}{'exon_track_space'}+${$config}{'exon'}{'exon_track_height'});
				if (${$config}{'intron'}{'intron_joiner_draw'}=~/^true$/i) {
					$this_mrna_y=$this_mrna_y+$num_mrna*${$config}{'intron'}{'intron_joiner_height'};
				}
				### joiner
				if (${$config}{'intron'}{'intron_joiner_draw'}=~/^true$/i or ${$config}{'intron'}{'intron_draw'}=~/^true$/i) {
					my $lastexonend=0;
					foreach my $ind_exon_start (sort {$a<=>$b} keys %{${$exons}{$ind_mrna}{'exon'}}) {
						my @exon_ends=keys %{${$exons}{$ind_mrna}{'exon'}{$ind_exon_start}};
						my $ind_exon_end=shift @exon_ends;
						if ($lastexonend==0) {
							$lastexonend=$ind_exon_end;
							next;
						}
						if (${$config}{'intron'}{'intron_joiner_draw'}=~/^true$/i) {
							my $xv = [$x_start+($lastexonend-$chrom_start{$ind_ref}+2)*$x_factor, $x_start+($lastexonend+$ind_exon_start-2*$chrom_start{$ind_ref})*$x_factor/2, $x_start+($ind_exon_start-$chrom_start{$ind_ref})*$x_factor];
							my $yv = [$this_mrna_y, $this_mrna_y-${$config}{'intron'}{'intron_joiner_height'}, $this_mrna_y];
							my $points = $vectorout->get_path(
							    x       => $xv,
							    y       => $yv,
							    -type   => 'polyline',
							    -closed => 'false' #specify that the polyline is closed.
							);
							$vectorout->polyline (
								%$points,
								id    =>"$ind_mrna-intron-$lastexonend-$ind_exon_start-joiner",
								style => {
									'fill-opacity' => 0,
									'stroke'       => ${$config}{'intron'}{'intron_joiner_stroke_color'}, 
									"stroke-width" => ${$config}{'intron'}{'intron_joiner_stroke_size'}
								}
							);
=pod
							$vectorout->line ( id=> "$ind_mrna-intron-$lastexonend-$ind_exon_start-left",
								x1 => $x_start+($lastexonend-$chrom_start{$ind_ref}+2)*$x_factor,
								y1 => $this_mrna_y, 
								x2 => $x_start+($lastexonend+$ind_exon_start-2*$chrom_start{$ind_ref})*$x_factor/2, 
								y2 => $this_mrna_y-${$config}{'intron'}{'intron_joiner_height'}, 
								stroke => ${$config}{'intron'}{'intron_joiner_stroke_color'}, 
								"stroke-width" => ${$config}{'intron'}{'intron_joiner_stroke_size'}
								);
							$vectorout->line ( id=> "$ind_mrna-intron-$lastexonend-$ind_exon_start-right",
								x1 => $x_start+($lastexonend+$ind_exon_start-2*$chrom_start{$ind_ref})*$x_factor/2,
								y1 => $this_mrna_y-${$config}{'intron'}{'intron_joiner_height'}, 
								x2 => $x_start+($ind_exon_start-$chrom_start{$ind_ref})*$x_factor,
								y2 => $this_mrna_y, 
								stroke => ${$config}{'intron'}{'intron_joiner_stroke_color'}, 
								"stroke-width" => ${$config}{'intron'}{'intron_joiner_stroke_size'}
								);
=cut
						}
						if (${$config}{'intron'}{'intron_draw'}=~/^true$/i) {
							$vectorout->line ( id=> "$ind_mrna-intron-$lastexonend-$ind_exon_start",
								x1 => $x_start+($lastexonend-$chrom_start{$ind_ref}+1)*$x_factor,
								y1 => $this_mrna_y+${$config}{'exon'}{'exon_track_height'}/2, 
								x2 => $x_start+($ind_exon_start-$chrom_start{$ind_ref}+1)*$x_factor, 
								y2 => $this_mrna_y+${$config}{'exon'}{'exon_track_height'}/2, 
								stroke => ${$config}{'intron'}{'intron_stroke_color'}, 
								"stroke-width" => ${$config}{'intron'}{'intron_stroke_size'}
							);
						}
						$lastexonend=$ind_exon_end;
					}
				}
				##UTR
				my $draw_utr5=0;
				my $draw_utr3=0;
				if (${$config}{'utr'}{'utr_draw'}=~/^true$/i) {
					if (exists ${$utr}{$ind_mrna} and exists ${$utr}{$ind_mrna}{'utr5'}) {
						$draw_utr5=1;
					}
					if (exists ${$utr}{$ind_mrna} and exists ${$utr}{$ind_mrna}{'utr3'}) {
						$draw_utr3=1;
					}
				}
				if ($test_utr and $draw_utr5) {
					foreach my $ind_utr_start (sort {$a<=>$b} keys %{${$utr}{$ind_mrna}{'utr5'}}) {
						my @utr_ends=keys %{${$utr}{$ind_mrna}{'utr5'}{$ind_utr_start}};
						my $ind_utr_end=shift @utr_ends;
						$vectorout->rectangle(x => $x_start+($ind_utr_start-$chrom_start{$ind_ref}+1)*$x_factor, 
						y => $this_mrna_y+${$config}{'exon'}{'exon_track_height'}/2-${$config}{'utr'}{'utr_track_height'}/2, 
						width  	=> ($ind_utr_end-$ind_utr_start+1)*$x_factor, 
						height => ${$config}{'utr'}{'utr_track_height'},
						id=> "$ind_mrna-utr5-$ind_utr_start-$ind_utr_end",
						style => {'fill' => ${$config}{'utr'}{'utr_fill_color'},
								'fill-opacity'   =>  ${$config}{'utr'}{'utr_fill_opaque'},
								'stroke'         => ${$config}{'utr'}{'utr_stroke_color'},
								'stroke-width'   => ${$config}{'utr'}{'utr_stroke_size'},
								'stroke-opacity' => ${$config}{'utr'}{'utr_stroke_opaque'},
							},
						);
					}
				}
				if ($test_utr and $draw_utr3) {
					foreach my $ind_utr_start (sort {$a<=>$b} keys %{${$utr}{$ind_mrna}{'utr3'}}) {
						my @utr_ends=keys %{${$utr}{$ind_mrna}{'utr3'}{$ind_utr_start}};
						my $ind_utr_end=shift @utr_ends;
						$vectorout->rectangle(x => $x_start+($ind_utr_start-$chrom_start{$ind_ref}+1)*$x_factor, 
						y => $this_mrna_y+${$config}{'exon'}{'exon_track_height'}/2-${$config}{'utr'}{'utr_track_height'}/2, 
						width  	=> ($ind_utr_end-$ind_utr_start+1)*$x_factor, 
						height => ${$config}{'utr'}{'utr_track_height'},
						id=> "$ind_mrna-utr3-$ind_utr_start-$ind_utr_end",
						style => {'fill' => ${$config}{'utr'}{'utr_fill_color'},
								'fill-opacity'   =>  ${$config}{'utr'}{'utr_fill_opaque'},
								'stroke'         => ${$config}{'utr'}{'utr_stroke_color'},
								'stroke-width'   => ${$config}{'utr'}{'utr_stroke_size'},
								'stroke-opacity' => ${$config}{'utr'}{'utr_stroke_opaque'},
							},
						);
					}
				}
				##Exon
				if (${$config}{'exon'}{'exon_draw'}=~/^true$/i and exists ${$cds}{$ind_mrna} and exists ${$cds}{$ind_mrna}{'cds'}) {
					foreach my $ind_cds_start (sort {$a<=>$b} keys %{${$cds}{$ind_mrna}{'cds'}}) {
						my @cds_ends=keys %{${$cds}{$ind_mrna}{'cds'}{$ind_cds_start}};
						my $ind_cds_end=shift @cds_ends;
						$vectorout->rectangle(x => $x_start+($ind_cds_start-$chrom_start{$ind_ref}+1)*$x_factor, 
						y => $this_mrna_y, 
						width  	=> ($ind_cds_end-$ind_cds_start+1)*$x_factor, 
						height => ${$config}{'exon'}{'exon_track_height'},
						id=> "$ind_mrna-CDS-$ind_cds_start-$ind_cds_end",
						style => {'fill' => ${$config}{'exon'}{'exon_fill_color'},
								'fill-opacity'   =>  ${$config}{'exon'}{'exon_fill_opaque'},
								'stroke'         => ${$config}{'exon'}{'exon_stroke_color'},
								'stroke-width'   => ${$config}{'exon'}{'exon_stroke_size'},
								'stroke-opacity' => ${$config}{'exon'}{'exon_stroke_opaque'},
							},
						);
					}
				}
				### Exon orientation
				if (${$config}{'exon'}{'exon_orientation_draw'}=~/^true$/i and exists ${$cds}{$ind_mrna} and exists ${$cds}{$ind_mrna}{'cds'}) {
					foreach my $ind_cds_start (sort {$a<=>$b} keys %{${$cds}{$ind_mrna}{'cds'}}) {
						my @cds_ends=keys %{${$cds}{$ind_mrna}{'cds'}{$ind_cds_start}};
						my $ind_cds_end=shift @cds_ends;
						my $cds_x1=$x_start+($ind_cds_start-$chrom_start{$ind_ref}+1)*$x_factor;
						my $cds_x2=$x_start+($ind_cds_end-$chrom_start{$ind_ref}+1)*$x_factor;
						my @x_coords=&GetMarkerCoords($cds_x1, $cds_x2, ${$config}{'exon'}{'exon_orientation_width'}, ${$config}{'exon'}{'exon_orientation_space'});
						foreach my $cds_idv_x (@x_coords) {
								my $xv = [$cds_idv_x,$cds_idv_x,$cds_idv_x];
								if (${$cds}{$ind_mrna}{'strand'} eq '+') {
									$xv = [$cds_idv_x-${$config}{'exon'}{'exon_orientation_width'}/2,$cds_idv_x+${$config}{'exon'}{'exon_orientation_width'}/2,$cds_idv_x-${$config}{'exon'}{'exon_orientation_width'}/2];
								}
								elsif (${$cds}{$ind_mrna}{'strand'} eq '-') {
									$xv = [$cds_idv_x+${$config}{'exon'}{'exon_orientation_width'}/2,$cds_idv_x-${$config}{'exon'}{'exon_orientation_width'}/2,$cds_idv_x+${$config}{'exon'}{'exon_orientation_width'}/2];
								}
								my $yv = [$this_mrna_y+${$config}{'exon'}{'exon_track_height'}/2-${$config}{'exon'}{'exon_orientation_height'}/2, $this_mrna_y+${$config}{'exon'}{'exon_track_height'}/2, $this_mrna_y+${$config}{'exon'}{'exon_track_height'}/2+${$config}{'exon'}{'exon_orientation_height'}/2];
								my $points = $vectorout->get_path(
									x       => $xv,
									y       => $yv,
									-type   => 'polyline',
									-closed => 'false' #specify that the polyline is closed.
								);
								$vectorout->polyline (
									%$points,
									id    =>"$ind_mrna-CDS_orientation-$cds_idv_x",
									style => {
										'fill-opacity' => 0,
										'stroke'       => ${$config}{'exon'}{'exon_orientation_color'}, 
										"stroke-width" => ${$config}{'exon'}{'exon_orientation_size'}
									}
								);
						}
					}
				}
				###
				$vectorout->text(x => $x_start+(${$mrnas}{$ind_mrna}{'start'}+${$mrnas}{$ind_mrna}{'end'}-2*$chrom_start{$ind_ref}+2)*$x_factor/2,
						y => $this_mrna_y+${$config}{'exon'}{'exon_track_height'}+${$config}{'main'}{'plot_font_size'}, 
						width => ${$config}{'main'}{'plot_font_size'}, 
						height => ${$config}{'main'}{'plot_font_size'}, 
						"font-family"=>${$config}{'main'}{'plot_font_family'}, 
						"text-anchor"=>"middle",
						"font-size"=>${$config}{'main'}{'plot_font_size'}, 
						"-cdata" => "$ind_mrna"
						);
			}
		}
	}
}
if (${$config}{'ruler'}{'ruler_draw'}=~/^true$/i) {
	my $ruler_length=substr ($max_x, 0, 1) * 10**(length($max_x)-2);
	if (exists ${$config}{'ruler'}{'ruler_stroke_seq_length'} and ${$config}{'ruler'}{'ruler_stroke_seq_length'}=~/^\d+$/) {
		$ruler_length=${$config}{'ruler'}{'ruler_stroke_seq_length'};
	}
	if ($verbose) {
		print "Max_x $max_x ruler $ruler_length\n";
	}
	
	$vectorout->line (id=> "ruler",
					x1 => ${$config}{'main'}{'plot_width'}-${$config}{'main'}{'plot_margin_right_size'}-$ruler_length*$x_factor,
					y1 => $total_y-${$config}{'main'}{'plot_margin_bottom_size'}/2, 
					x2 => ${$config}{'main'}{'plot_width'}-${$config}{'main'}{'plot_margin_right_size'}, 
					y2 => $total_y-${$config}{'main'}{'plot_margin_bottom_size'}/2, 
					stroke => ${$config}{'ruler'}{'ruler_stroke_color'}, 
					"stroke-width" => ${$config}{'ruler'}{'ruler_stroke_size'}
				);
	$vectorout->text(	x => ${$config}{'main'}{'plot_width'}-${$config}{'main'}{'plot_margin_right_size'}-$ruler_length*$x_factor/2, 
					y => $total_y-${$config}{'main'}{'plot_margin_bottom_size'}/2+${$config}{'main'}{'plot_font_size'}, 
					width => ${$config}{'main'}{'plot_font_size'}, 
					height => ${$config}{'main'}{'plot_font_size'}, 
					"font-family"=>${$config}{'main'}{'plot_font_family'}, 
					"text-anchor"=>"middle",
					"font-size"=>${$config}{'main'}{'plot_font_size'}, 
					"-cdata" => "$ruler_length bp");
}
my $finalout = $vectorout->xmlify;
open SVGFILE, "> $output";
print SVGFILE $finalout;
close SVGFILE;


#####################################################################
###                         sub functions                         ###
#####################################################################

sub GetMarkerCoords {
	my ($GMCx1, $GMCx2, $GMClen, $GMCspace)=@_;
	
	if ($debug) {
		print "    Test: $GMCx1, $GMCx2, $GMClen, $GMCspace\n";
	}
	my $GMCmiddle=($GMCx1+$GMCx2)/2;
	my @GMCret_arr=();
	
	my $GMAright_border=$GMCx2-max($GMClen, $GMCspace)/2;
	my $GMAleft_border=$GMCx1+max($GMClen, $GMCspace)/2;
	
	if (($GMCx2-$GMCx1)>$GMClen) {
		push (@GMCret_arr, $GMCmiddle);
	}
	if ($debug) {
		print "      Coords: ",join(",", @GMCret_arr),"\n";
	}
	my $GMCborder=$GMCmiddle-$GMCspace-$GMClen-$GMClen/2;
	while ($GMCborder>=$GMAleft_border) {
		unshift (@GMCret_arr, $GMCborder+$GMClen/2);
		$GMCborder=$GMCborder-$GMCspace-$GMClen;
	}
	if ($debug) {
		print "      Coords: ",join(",", @GMCret_arr),"\n";
	}
	my $GMCborder=$GMCmiddle+$GMCspace+$GMClen+$GMClen/2;
	while ($GMCborder<=$GMAright_border) {
		push (@GMCret_arr, $GMCborder-$GMClen/2);
		$GMCborder=$GMCborder+$GMCspace+$GMClen;
	}
	
	if ($debug) {
		print "      Coords: ",join(",", @GMCret_arr),"\n";
	}
	return @GMCret_arr;
}

exit 0;
