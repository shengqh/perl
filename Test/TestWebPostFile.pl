use strict;
use WWW::Mechanize;

my $url = "http://rbpmap.technion.ac.il/";

my $agent = WWW::Mechanize->new(
  autocheck => 1,
  agent     => "Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/47.0.2526.73 Safari/537.36"
);
$agent->get($url);

my $form = $agent->form_name("RBPmap_form");

$form->find_input('genome')->value("human");
$form->find_input('db')->value("hg19");
$form->find_input('seq_method')->value("file");
$form->find_input('seq_file')->disabled(0);
$form->find_input('seq_file')->value("E:/PUM2pos_test.fa");
$form->find_input('sequence')->disabled("disabled");
$form->find_input('db_selection')->value("on");
$form->find_input('motifs_selection')->value("None");
$form->find_input('is_user_pssm')->value("on");
$form->find_input('pssm_file')->disabled(0);
$form->find_input('pssm_file')->value("E:/PUM2.pwm.txt");
$form->find_input('stringency')->value("medium");
$form->find_input('is_conservation')->value("on");
$form->find_input('is_email')->value("TRUE");
$form->find_input('email')->disabled(0);
$form->find_input('email')->value("quanhu.sheng\@vanderbilt.edu");
$form->find_input('job_name')->value("JobByP");

$agent->submit_form();

die unless ( $agent->success );

if($agent->content() =~ /The calculation time depends on the number of sequences/){
  print "Submit job succeed, the result page = ", $agent->uri();
}else{
  print "Submit job failed, the whole result = ", print $agent->content();
}
