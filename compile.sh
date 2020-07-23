# This source file is part of the Global Atmospheric Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

echo "Starting to compile ndvar."
mpicc src/* -leccodes -lnetcdf -lm -lgeos95 -latmostracers -Wl,-rpath=/lib -Wall -o ndvar
if [ $? -ne 0 ]
then
echo -e ${RED}Ndvar compilation failed.$NC
else
echo "ndvar compiled sucessfully."
fi
