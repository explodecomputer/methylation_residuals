mkdir splitfiles

split -a 3 -l 1000 -d 15upnorm.phen splitfiles/15upnorm.phen.
split -a 3 -l 1000 -d antenatalnorm.phen splitfiles/antenatalnorm.phen.
split -a 3 -l 1000 -d cordnorm.phen splitfiles/cordnorm.phen.
split -a 3 -l 1000 -d F7norm.phen splitfiles/F7norm.phen.
split -a 3 -l 1000 -d FOMnorm.phen splitfiles/FOMnorm.phen.
