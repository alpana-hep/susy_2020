#!/bin/bash
# a=$1 #"Hello i am pass"; h_Sbins_LL_MET_200
# if [ `echo $hist1 | grep -c "MET200" ` -gt 0 ]
# then
#     echo "Success"
# elif [`echo $hist1 | grep -c "MET100" ` -gt 0]
# then
#      echo "Success"

# elif [`echo $hist1 | grep -c "MET300" ` -gt 0]
# then
#      echo "Success"

# else
#   echo "Fail";
# fi

#!/bin/bash

a=$1
mass_p=100
if [ $(echo $a | grep -c "MET200") -gt 0 ]
then
    echo "Success"
    mass_p=200
    echo $mass_p
elif [ $(echo $a | grep -c "MET100") -gt 0 ]
then
    echo "Success"
elif [ $(echo $a | grep -c "MET300") -gt 0 ]
then
    echo "Success"
else
    echo "Fail"
fi
