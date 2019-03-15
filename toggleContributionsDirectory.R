system("if [ -L contributions/bjContributions ] ; then 
echo 'it is here' ; 
rm contributions/bjContributions else echo 'it is not here'; 
ln -s ../bjContributions contributions/bjContributions
       fi")
