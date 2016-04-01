#!/usr/bin/perl

open(FILEIN,"laplace.aei");
while(<FILEIN>)
{
    my @vars = split;
    while($vars[1]<-180.)
    {
	$vars[1] = $vars[1] + 360.;
    }
    while($vars[1]>180.)
    {
	$vars[1] = $vars[1] - 360.;
    }
    print "$vars[0] $vars[1] \n";
}
