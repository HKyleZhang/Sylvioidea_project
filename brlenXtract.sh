#!/bin/bash

read -p "Total Branch Length is extracted with current folder? (y/n) " answer
if [ "${answer}" = "y" ]; then

  read -p "
Specify a taxon name or type 'all': " taxa
  if [ "${taxa}" != " " ]; then

    if [ "${taxa}" != "all" ]; then

      for file in *; do
        filetype=$(echo "${file}" | cut -d "." -f 2)
        if [ "${filetype}" != "sh" ]; then
          name=$(echo "${file}" | cut -d "." -f 1)
          brlen=$(nw_distance -n ${file} | grep "${taxa}" | cat | sed "s/${taxa}\t//")
          echo -e "${brlen}\t${name}" >>${taxa}_brlen
        fi
      done

      echo "
Job done!"

    elif [ "${taxa}" = "all" ]; then

      for file in *; do
        filetype=$(echo "${file}" | cut -d "." -f 2)
        if [ "${filetype}" != "sh" ]; then
          name=$(echo "${file}" | cut -d "." -f 1)
          brlensum=$(nw_distance -n ${file}) >${name}_brlensum
          brlen1=$(nw_distance -n ${file} | grep "ZF" | cat | sed "s/ZF\t//")
          brlen2=$(nw_distance -n ${file} | grep "BC" | cat | sed "s/BC\t//")
          brlen3=$(nw_distance -n ${file} | grep "CW" | cat | sed "s/CW\t//")
          brlen4=$(nw_distance -n ${file} | grep "MW" | cat | sed "s/MW\t//")
          brlen5=$(nw_distance -n ${file} | grep "GRW" | cat | sed "s/GRW\t//")
          brlen6=$(nw_distance -n ${file} | grep "CRW" | cat | sed "s/CRW\t//")

          nw_distance -n ${file} >${name}_brlensum
          echo -e "${brlen1}\t${name}" >>ZF_brlen
          echo -e "${brlen2}\t${name}" >>BC_brlen
          echo -e "${brlen3}\t${name}" >>CW_brlen
          echo -e "${brlen4}\t${name}" >>MW_brlen
          echo -e "${brlen5}\t${name}" >>GRW_brlen
          echo -e "${brlen6}\t${name}" >>CRW_brlen
        fi
      done

      mkdir branch_length
      mkdir brlen_sum
      mv ./*brlensum ./brlen_sum/
      mv ./brlen_sum ./*brlen ./branch_length/

      echo "
Job done!"

    fi

  fi

fi
