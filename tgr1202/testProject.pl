#!/usr/bin/perl 
use strict;
use warnings;
use Getopt::Long;
use JSON;

my  $dirroot = 'E:/nodejs/tgr1202';

print "dirroot = $dirroot \n";

my $config = $dirroot . "/app/common/config.js";

open FILE, "$config" or die "Couldn't open file: $config";
my $string = join("", <FILE>);
close FILE;

my ($server) = $string =~ /cqsdbUrl:\s'(\S+)'/;
my ($applicationid) = $string =~ /applicationId:\s'(\S+)'/;

print "server = $server \n";
print "applicationid = $applicationid \n";

my $ret = `curl -d username=foo -d password=bar $server/auth/login`;
my ($token) = $ret =~ /"key":"(\S+)"/;

print "token = $token \n";

my $project = `curl -H "Authorization: $token" $server/projects/$applicationid/collections/StudyEnrollments/54fdcc99eb3b1f600f32cb9c`;

my $json_print = JSON->new;
print $json_print->pretty->encode($json_print->decode($project));

