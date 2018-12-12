#!/bin/bash

curdir=$(pwd)

rm -rf ${curdir}/KC_Distance
mkdir ${curdir}/KC_Distance

for file in ${curdir}/Cluster_Assignment/*assign; do
  name=$(basename ${file} | cut -d "." -f 1 | cut -d "-" -f 1)

  lnum=$(cat ${file} | wc -l)
  l=1
  while [ "${l}" -le "${lnum}" ]; do
    cl=$(cat ${file} | sed -n "${l}p")
    if [ "${cl}" -lt 10 ]; then clmd=$(echo "cluster0${cl}"); else clmd=$(echo "cluster${cl}"); fi
    if [ ! -d ${curdir}/KC_Distance/${clmd} ]; then mkdir ${curdir}/KC_Distance/${clmd}; fi
    dis=$(cat ${curdir}/KCmetric_distance/${name}/${name}-${clmd}.kcd | sed -n "${l}p")
    echo "${name} ${dis}" >>${curdir}/KC_Distance/${clmd}/${name}.dis
    ((l++))
  done
done

clnum=$(ls ${curdir}/KC_Distance | wc -l)
rm -rf ${curdir}/KC_Distance_summary
mkdir ${curdir}/KC_Distance_summary
lnum=$(cat ${curdir}/order | wc -l)
l=1
while [ "${l}" -le "${lnum}" ]; do
  gene=$(cat ${curdir}/order | sed -n "${l}p")
  cl=1
  while [ "${cl}" -le "${clnum}" ]; do
    if [ "${cl}" -lt 10 ]; then clmd=$(echo "cluster0${cl}"); else clmd=$(echo "cluster${cl}"); fi
    touch ${curdir}/KC_Distance_summary/${clmd}.dsum
    fe=$(ls ${curdir}/KC_Distance/${clmd} | grep "${gene}")
    if [ -n "${fe}" ]; then
      sed -i "s/$/ ${l}/" ${curdir}/KC_Distance/${clmd}/${gene}.dis
      cat ${curdir}/KC_Distance_summary/${clmd}.dsum ${curdir}/KC_Distance/${clmd}/${gene}.dis >>${curdir}/KC_Distance_summary/${clmd}-temp.dsum
      mv ${curdir}/KC_Distance_summary/${clmd}-temp.dsum ${curdir}/KC_Distance_summary/${clmd}.dsum
    else
      echo "${gene} NA ${l}" >>${curdir}/KC_Distance_summary/${clmd}.dsum
    fi
    ((cl++))
  done
  ((l++))
done
