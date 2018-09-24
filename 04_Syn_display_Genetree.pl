#!/bin/perl
use strict;
use warnings;
use Bio::Graphics;
use Bio::SeqFeature::Generic;
use Bio::Graphics::Glyph::segments;
use Bio::Graphics::Glyph::heat_map;
use Bio::Graphics::Glyph::generic;

#usage: perl 04_Syn_display_Genetree.pl SEQ-file.synt3.GFF.unify SEQ-file.pairs SEQ-file.synt.longest_prot.fa.bp-clustI14.out GENETREE

my $inputGFF=$ARGV[0];
my $PAIRS_file=$ARGV[1];
my $CLUSTER_file=$ARGV[2];
my $GENETREE=$ARGV[3];

my $totalNumClusters=`grep -P '\\t' $CLUSTER_file | wc -l`;
print "number of clusters: $totalNumClusters\n";

my $Scaffolds=`cut -f1 -d' ' $inputGFF | sort | uniq`;
chomp $Scaffolds;

my @scaffolds = split(/\n/, $Scaffolds);
print "scaffolds:\n@scaffolds\n";
####################################################
my $GENENAMES_Tree=`sed -r 's/\\)|\\(|:|,/\\n/g' $GENETREE| sed "s/'//g" |  grep -P "^\\w\\w\\w\\w_" | xargs -I % grep '%' $PAIRS_file | cut -f3`;

#my $GENENAMES_Tree=`sed -r 's/\)|\(|:|,/\n/g' $GENETREE | grep -P "^\w\w\w\w_" | xargs -I % grep '%' $PAIRS_file | cut -f3` ;

my @newScaffolds= split(/\n/, $GENENAMES_Tree);

print "new Scaffolds:\n@newScaffolds\n\n";

sleep 2;
########################################

my $output=$inputGFF.'_SynStack.svg';
open (my $SVGOUT, '>', $output) or die "$! could not open filehandle to file $output\n\n";


#go through new gff collection and for each scaffold. calculate the maximum span
#fürs erste nehmen wir 300.000 

#check max span in each scaffold
my %lows;
my %dis2start_collect;
my $totalDist=0;
my $topDist2startGOI=0;
my $topDist2EndGOI=0;
my $topGOIsize=0;

for(@newScaffolds){

my $GeneOfInterest = `grep "$_" $PAIRS_file | cut -f2  `;
chomp $GeneOfInterest;

#total span
my $lowest_start = `grep "$_" $inputGFF | sort -n -k4 | head -1 | cut -d' ' -f4`;
chomp $lowest_start;
$lows{$_}=$lowest_start;

my $highest_stop=`grep "$_" $inputGFF | sort -nr -k5 | head -1 | cut -d' ' -f5`;
chomp $highest_stop;

#my $span = $highest_stop - $lowest_start;
#if ($span > $totalDist ){$totalDist = $span;}

#GOI position
my $start_GOI= `grep "$_" $inputGFF | grep "ID=$GeneOfInterest" | cut -d' ' -f4`;
chomp $start_GOI;
my $stop_GOI = `grep "$_" $inputGFF | grep "ID=$GeneOfInterest" | cut -d' ' -f5`;
chomp $stop_GOI;
my $sizeGOI=$stop_GOI - $start_GOI;
if ($sizeGOI>$topGOIsize){print "we have a higher GOIsize\n";$topGOIsize=$sizeGOI ;}

$dis2start_collect{$_}=$start_GOI - $lowest_start;
if ($dis2start_collect{$_} > $topDist2startGOI){$topDist2startGOI= $dis2start_collect{$_}; }
   
my $dis2stopGOI=$highest_stop - $stop_GOI ;
if ($dis2stopGOI > $topDist2EndGOI ){ print "yes\n";$topDist2EndGOI = $dis2stopGOI ; }

###
print "goi: $GeneOfInterest\tsize: $sizeGOI\ttopsize: $topGOIsize\tdistance to start: $dis2start_collect{$_}\ttop distance to start: $topDist2startGOI\tdistance to stop: $dis2stopGOI\t top dis to stop: $topDist2EndGOI\n";

}
$totalDist = $topGOIsize + $topDist2startGOI + $topDist2EndGOI;


#build the first panel
my $panel=Bio::Graphics::Panel->new(-length     =>  $totalDist,
                                    -width      =>  2500,
                                    -pad_left   =>  200,
                                    -pad_right  =>  50,
                                    -key_style  =>  "left",
                                    -image_class=>  'GD::SVG',
                                    -key_font=>"gdLargeFont",
                                    -spacing=>20, #default 5
                                    
                                    );

                            
#start building tracks
my $full_length=Bio::SeqFeature::Generic->new(-start=>0,-end=>$totalDist);
$panel->add_track(
                $full_length,
		  -glyph=>'arrow',
		  -tick=>2,
		  -fgcolor=>'black',
		  -double=>1,
		  -font=>"gdLargeFont",
		  -label_intervals=>0, ####?
		  
		 
		  );
		  
  
my @track;
#for each scaffold make a new panel

for (my $i=0;$i<scalar(@newScaffolds);$i++){
my $scaf= $newScaffolds[$i];
print "$scaf $i\n";
my $factor2center=$topDist2startGOI-$dis2start_collect{$scaf};

$track[$i]=$panel->add_track(
			  -glyph =>'heat_map',
			  -label => 1,
			  -strand_arrow=>1,
			  -key=>"$scaf",
			  -min_score=>0,
			  -max_score=>$totalNumClusters +30 , #scale is dependent on number of clusters with more than 1 sequence....
			  -start_color=>"grey",
			  -end_color=>"red",
			  -font=>"gdLargeFont", #gdGiantFont
			   -height=>18,
			   -label_position=>'top', ###########?
			  );

			  
my $PART = `grep "$scaf" $inputGFF `;
my @lines = split(/\n/, $PART);
my $count=0;



for my $LINE (@lines){#

    $count++;
    my $score=0;
    print "$count $LINE\n";
    my @Elements = split(/ /, $LINE);
    #my $ID = `echo $Elements[8] | tr ";" "\n" | grep '^ID=' -  `;
    my @subElements = split(/;/, $Elements[8]);
    my $ID=$subElements[0];
    $ID=~s/ID=//;
    chomp $ID;
    $ID =~ s/^ID=//;
    my $name = $ID;
    #my $original = `grep "\t$ID" $PAIRS_file`;
     my @original = split(/__/, $scaf);

    my $cluster = `grep -P "\t" $CLUSTER_file | grep -n "$ID" | cut -f1 -d':'`;
    chomp $cluster;
   #if($ID=~m/dmel/i){print "$ID\t$cluster\n";}
    if ($original[0] eq $ID ){$name = "$ID";$score=$totalNumClusters +30;}
    elsif ($cluster ){$name = "C$cluster"; $score=$cluster+20;}
    else{$name ='';}
    
    #print "seq_id $ID\tstart $Elements[3]\tend $Elements[4]\tstrand $Elements[6]\tscore $score\n";
    
    my $feature = Bio::SeqFeature::Generic->new(
                                                    -seq_id=>$name, #no sequence ID---- this should be optional....
						    -start=>$Elements[3]-$lows{$scaf}+$factor2center,
						    -end=>$Elements[4]-$lows{$scaf}+$factor2center,
						    -strand=>$Elements[6],
						    -score=>$score
						    );

    $track[$i]->add_feature($feature);
}



}

print $SVGOUT $panel->svg;