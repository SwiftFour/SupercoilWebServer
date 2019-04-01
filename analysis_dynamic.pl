#! /usr/bin/perl -w
use strict;

if (!$ARGV[0]){
	print "perl analysis_dynamic.pl expression-fast-0.6.dynamic\n\n";
	exit(1);
}
open (fp1, "$ARGV[0]") || die "$!";
chomp ( my @dynamic=<fp1>);
close fp1;

my $startNum =0;
my $boundayNum = 0;
my $hegNum = 0;
my $napNum = 0;
my $hegNum1 = 0;
my $napNum1 = 0;
my $hegNapNum =0;
my $hegNapNum1 =0;
my $obValue = 0;
my $perturbation = 0;
my $hegNumR= 0;
my $napNumR = 0;
my $hegNapNumR =0;
my $flagStart = 0;
my $flagR = 0;
my $parameter = 0;
my $parameter1 = 0;
my $parameter2 = 0;
my $parameter3 = 0;
my $parameter4 = 0;


print "startNum\tp1\tp2\tp3\tp4\tp5\tobValue\tperturbation\tboundayNum\thegNum\tnapNum\thegNapNum\tPhegNum\tPnapNum\tPhegNapNum\thegNumR\tnapNumR\thegNapNumR\n";

for (my $i=0; $i<@dynamic; $i++){
	if ($dynamic[$i] =~ /Start/ && !$flagStart){
		chomp(my @a=split(/\t/, $dynamic[$i]));
		chomp(my @aa=split(/ /, $dynamic[$i]));
		$obValue = $a[1];
		$parameter = $a[2];
		$parameter1 = $a[3];	
		$parameter2 = $a[4];	
		$parameter3 = $a[5];	
		$parameter4 = $a[6];	
		$startNum = $aa[3];
		$flagStart = 1;
		next;
	}
	if ($dynamic[$i] =~ /Random/){
		$flagR = 1;
		next;
	}
	if ($dynamic[$i] =~ /Perturbation/){
		chomp(my @a=split(/ /, $dynamic[$i]));
		$perturbation = $a[5];
		next;
	}
	chomp(my @b=split(/\t/,$dynamic[$i]));
	if (!$flagR){
		if ($b[6]){
			$hegNum++;
		}
		if ($b[7]){
			$napNum++;
		}
		if ($b[6] && $b[7]){
			$hegNapNum++;
		}
		if ($b[8]){
			$hegNum1++;
		}
		if ($b[9]){
			$napNum1++;
		}
		if ($b[8] && $b[9]){
			$hegNapNum1++;
		}
		$boundayNum++;
	}else{
		if ($b[6]){
			$hegNumR++;
		}
		if ($b[7]){
			$napNumR++;
		}
		if ($b[6] && $b[7]){
			$hegNapNumR++;
		}
	}
	if ($dynamic[$i] =~ /Start/ && $flagStart){
		print "$startNum\t$parameter\t$parameter1\t$parameter2\t$parameter3\t$parameter4\t$obValue\t$perturbation\t$boundayNum\t$hegNum\t$napNum\t$hegNapNum\t$hegNum1\t$napNum1\t$hegNapNum1\t$hegNumR\t$napNumR\t$hegNapNumR\n";
		$flagR = 0;
		$hegNum =0;
		$napNum = 0;
		$hegNum1 =0;
		$napNum1 = 0;
		$hegNumR =0;
		$napNumR = 0;
		$boundayNum = 0;
		$hegNapNum = 0;
		$hegNapNum1 = 0;
		$hegNapNumR =0;
		chomp(my @a=split(/\t/, $dynamic[$i]));
		chomp(my @aa=split(/ /, $dynamic[$i]));
		$obValue = $a[1];
		$parameter = $a[2];
		$parameter1 = $a[3];	
		$parameter2 = $a[4];	
		$parameter3 = $a[5];	
		$parameter4 = $a[6];	
		$startNum = $aa[3];
		next;
	}
}
print "$startNum\t$parameter\t$parameter1\t$parameter2\t$parameter3\t$parameter4\t$obValue\t$perturbation\t$boundayNum\t$hegNum\t$napNum\t$hegNapNum\t$hegNum1\t$napNum1\t$hegNapNum1\t$hegNumR\t$napNumR\t$hegNapNumR\n";
