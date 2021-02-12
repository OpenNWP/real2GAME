#!/bin/bash

# This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
# Github repository: https://github.com/MHBalsmeier/game

# This script can be used to replace $old_string by $new_string in all files of the project.

old_string="https:\/\/github.com\/MHBalsmeier\/ndvar"
new_string="https:\/\/github.com\/AUN4GFD\/ndvar"
grep -rl $old_string . | xargs sed -i 's/'$old_string'/'$new_string'/g'
