
use strict;
use warnings;

require "DownloadGEO.pl";

#my @cellline=("GSE7561", "GSE7848","GSE8565",GSE8597,GSE11324,GSE11506,GSE11352,GSE9187,GSE9586,GSE6800,GSE8087,GSE8562,GSE8096,GSE7382,GSE4006,GSE3188,GSE7765,GSE6462,GSE6883,GSE4025,GSE4668 "GSE16200");
#GSE4917,GSE3542,GSE2225,GSE4292,GSE1647,GSE3529,GSE2222, GSE2251, GSE1400,GSE1045,GSE848,GSE763

#my @datasets = (
#  "GSE1378", "GSE1379", "GSE1456", "GSE1561", "GSE1992", "GSE2034", "GSE2109",
#  "GSE2429", "GSE2603", "GSE2607", "GSE2740", "GSE2741", "GSE2990", "GSE3013",
#  "GSE3143", "GSE3165", "GSE3494", "GSE3521", "GSE3744", "GSE3893", "GSE4779",
#  "GSE4913",
#  "GSE4922", "GSE5327", "GSE5364", "GSE5460", "GSE5462", "GSE5847", "GSE5764",
#  "GSE5847", "GSE6128", "GSE6130", "GSE6434", "GSE6532", "GSE6577", "GSE6772",
#  "GSE6861", "GSE7390", "GSE7849", "GSE7904", "GSE8465", "GSE9195",
#  "GSE10270", "GSE10797", "GSE10886", "GSE10893", "GSE11121",
#  "GSE12093", "GSE12276", "GSE12763", "GSE16446", "GSE17705", "GSE18229",
#  "GSE18864",
#  "GSE19615", "GSE20194", "GSE20437", "GSE20711", "GSE21653", "GSE21997",
#  "GSE22226", "GSE22358", "GSE23720", "GSE31448"
#);

my @datasets = (
  "GSE25066"
);

DownloadGeoDatasets(@datasets);

print "Finished!\n";
