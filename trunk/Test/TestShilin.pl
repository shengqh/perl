use WWW::Mechanize
$w = WWW::Mechanize->new;

#$w->get("https://mail.yahoo.com");
$w->get("https://apps.tn.gov/dlappts-app/appt.jsp");

print $w->uri; 
print $w->content;