#!/bin/sh
#Script to extract the unique gene IDs from the module files
#Run a loop through all the module files merging the columns with gene ids into one single column and keep
#only unique ids.
for MODULE in input/*_modules/*; do
   cat ${MODULE} | awk '{print $1}' > tmp/tmp_ids
   cat ${MODULE} | awk '{print $2}' >> tmp/tmp_ids
   MODULE_NAME=$(echo $MODULE | cut -d '/' -f 3 | cut -d '.' -f 1)
   sort -u tmp/tmp_ids > tmp/${MODULE_NAME}_unique_names.ids
   rm -f tmp/tmp_ids
done
