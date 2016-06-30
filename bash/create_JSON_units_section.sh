#!/bin/bash
for i in *; do cont1=$(echo $(ls ${i})|cut -f 1 -d " "); cont2=$(echo $(ls ${i})|cut -f 2 -d " "); echo "\""$i"\" : [\""$cont1"\", \""${cont2}"\"]," ; done
