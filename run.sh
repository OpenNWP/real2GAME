# This source file is part of NDVAR (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/ndvar

echo "Starting the assimilation process ..."
mpirun -np 1 ./bin/ndvar $year $month $day $hour
# valgrind ./test_generator $test_id
if [ $? -ne 0 ]
then
echo -e ${RED}Data assimilation failed.$NC
else
echo "Model input file created sucessfully."
fi
