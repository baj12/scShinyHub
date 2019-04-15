system("if [ -L contributions/bjContributions ] ; then 
echo 'it was here' ; 
rm contributions/bjContributions; else echo 'it was not here'; 
ln -s ../bjContributions contributions/bjContributions
       fi")

