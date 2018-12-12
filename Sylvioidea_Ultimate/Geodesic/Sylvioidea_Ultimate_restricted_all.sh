#!/bin/bash

#Functions
usage()
{
 echo "Usage: [-a] <Alignment folder>"
}

IQTREE()
{
#Execute the analysis
 iternum=2
 iterindex=1
 while [ "${iterindex}" -le "${iternum}" ]
 do
   rm -rf ${1}/.run${iterindex}
   mkdir ${1}/.run${iterindex}
   cp -p ${2} ${1}/.run${iterindex}
   
   cd ${1}/.run${iterindex}
   for iterfile in *
   do   
     seqlen=$(sed -n "1p" ${iterfile} | cut -d " " -f 3)
     if [ "${seqlen}" -lt 10000 ]; then 
         nth="1"
     elif [ "${seqlen}" -lt 1000000 ]; then 
         nth="2"
     elif [ "${seqlen}" -lt 10000000 ]; then
         nth="3"
     else
         nth="4"
     fi
     if [ "${3}" == "-bb" ]; then
        iqtree -s ${iterfile} -m TEST -AICc -ninit 200 -ntop 50 -nt ${nth} -keep-ident ${3} ${4} ${5} ${6} > /dev/null
     else
        iqtree -s ${iterfile} -m TEST -AICc -ninit 200 -ntop 50 -nt ${nth} -keep-ident > /dev/null
     fi
    done
    ((iterindex++))
  done 

#Identify the run with bigger log likelihood
  cd ${1}  
  iterindex=1
  while [ "${iterindex}" -le "${iternum}" ]
  do
    loglkh[${iterindex}]=$(cat ${1}/.run${iterindex}/*iqtree | grep "Log-likelihood of the tree:" | cut -d " " -f 5)
    loglkhmd[${iterindex}]=$(echo "${loglkh[${iterindex}]} * 10000" | bc | cut -d "." -f 1 | cut -d "-" -f 2)
    ((iterindex++))
  done
  
  iterindex=2; loglkht=${loglkhmd[1]}; tnum=1
  while [ "${iterindex}" -le "${iternum}" ]
  do
    if [ "${loglkht}" -gt "${loglkhmd[${iterindex}]}" ]; then
       loglkht=${loglkhmd[${iterindex}]}
       tnum=${iterindex}
    fi
    ((iterindex++))
  done
  
#File sorting
  if [ "${3}" == "-bb" ]; then
     mv ${1}/.run${tnum}/*ufboot ${1}/${7}/bs-file
     mv ${1}/.run${tnum}/*iqtree ${1}/${7}/report/
     mv ${1}/.run${tnum}/*treefile ${1}/${7}/
     mv ${1}/.run${tnum}/*log ${1}/.run${tnum}/*gz ${1}/.run${tnum}/*bionj ${1}/.run${tnum}/*mldist ${1}/${7}/trivia
  else
     mv ${1}/.run${tnum}/*iqtree ${1}/${4}/report/
     mv ${1}/.run${tnum}/*treefile ${1}/${4}/
     mv ${1}/.run${tnum}/*log ${1}/.run${tnum}/*gz ${1}/.run${tnum}/*bionj ${1}/.run${tnum}/*mldist ${1}/${4}/trivia
  fi
  rm -rf ${1}/.run*
}

aftertreebuild()
{
 rm -rf ${1}/${2}_rooted
 mkdir ${1}/${2}_rooted
 #mkdir ${1}/treeXtract
 #mkdir ${1}/treeXtract/${2}
 #mkdir ${1}/treeXtract/${2}_rooted

 for file in ${1}/${2}/*treefile
 do
   namepre=$(echo "${file}")
   name=$(basename ${namepre})
   nw_reroot ${file} ZF > ${1}/${2}_rooted/${name}
 done
 
 if [ "${3}" -gt 1 ]; then
    #for file in ${1}/${2}/*treefile
    #do
    #   namepre=$(echo "${file}")
    #   name=$(basename ${namepre} | cut -d "." -f 1 | cut -d "-" -f 2)
    #   seq=$(sed -n "1p" ${file})
    #  
    #   echo "$name $seq" >> ${1}/treeXtract/${2}/gnTrees_collection.tre
    #   echo "$seq" >> ${1}/treeXtract/${2}/gnTrees_noname.tre
    #done

    for file in ${1}/${2}_rooted/*treefile
    do
      namepre=$(echo "${file}")
      if [ "${3}" -gt 5 ]; then
         name=$(basename ${namepre} | cut -d "." -f 1 | cut -d "_" -f 3 | cut -d "-" -f 1)
      else    
         name=$(basename ${namepre} | cut -d "." -f 1 | cut -d "-" -f 2)
      fi 
      seq=$(sed -n "1p" ${file}) 
 
      echo "$name $seq" >> ${1}/gnTrees_collection.tre
   #  echo "$name $seq" >> ${1}/treeXtract/${2}_rooted/gnTrees_collection.tre
   #  echo "$seq" >> ${1}/treeXtract/${2}_rooted/gnTrees_noname.tre
   done
 fi
}

XtractTree()
{
 for file in ${1}/${2}/*treefile
 do
    namepre=$(echo "${file}")
    name=$(basename ${namepre} | cut -d "." -f 1 | cut -d "-" -f 2)
    seq=$(sed -n "1p" ${file})
  
    #if [ ${index} -eq 1 ]; then
    #   echo "$name $seq" >> ${1}/gnTrees_collection.tre
    #elif [ ${index} -eq 2 ]; then
    #   echo "$seq" >> ${1}/gnTrees_noname.tre
    #elif [ ${index} -eq 3 ]; then
       echo "$name $seq" >> ${1}/gnTrees_collection.tre
    #   echo "$seq" >> ${1}/gnTrees_noname.tre
    #fi
 done
}

scrutinize()           #Function for picking out the trees wiht at least one boostrapping value < threshold value.
{
 stlen=$(echo "${#2}")
 if [ "${3}" -lt 100 ]; then
    fiend=5
 elif [ "${3}" -lt 1000 ]; then
    fiend=6
 elif [ "${3}" -lt 10000 ]; then
    fiend=7
 fi  
 step=$((fiend-1))
 while [ "${fiend}" -le "${stlen}" ]
 do
   fist=$((fiend - step))
   unit=$(echo "${2}" | cut -c ${fist}-${fiend})
   left=$(echo "${unit}" | grep ")")
   right=$(echo "${unit}" | grep ":")
   if [ -n "${left}" ] && [ -n "${right}" ]; then
      unit=$(echo "${unit}" | cut -d ")" -f 2 | cut -d ":" -f 1)
      if [ -n "${unit}" ]; then
         if [ "${unit}" -eq "${unit}" ] 2>/dev/null; then
            if [ "${unit}" -lt "${3}" ]; then
               res="FAIL"
            fi
         fi
      fi
   fi
   ((fiend++))
 done
}

tree_scrutinizer()
{
 thva=$(echo "${2}")
 if [ -z "${thva}" ]; then
    thva=70
 else
    if [ "${thva}" -eq "${thva}" ] 2>/dev/null; then   
       if [ "${thva}" -eq 1 ]; then
          thva=70
       elif [ "${thva}" -eq 2 ]; then
            read -p $'\nSpecify the value: ' thva
            if [ -z "${thva}"  ]; then echo -e "\nERROR!" && exit; fi 
       fi
    else
       echo -e "\nERROR!" && exit
    fi      
 fi
      
 lnum=$(cat ${1}/gnTrees_collection.tre | wc -l)
 i=0
 while [ "${i}" -lt "${lnum}" ]
 do
   ((i++))
  
   #if [ $((i%7)) -eq 0 ]; then echo -e "Processing..."; fi
  
   sentence=$(cat ${1}/gnTrees_collection.tre | sed -n "${i}p")
   name=$(echo "${sentence}" | cut -d " " -f 1)
   res="PASS"
   scrutinize ${sentence} ${thva}
   if [ "${res}" = "FAIL" ]; then
      echo "${sentence}" >> ${1}/output_list
   fi
 done

 if [ -e output_list ]; then
    echo "
     ********************************************************
         Msg:The results has been written to output_list.    
     ********************************************************"
 else
    echo "
     ********************************************************
               Msg:No untrustworthy tree is found.    
     ********************************************************"
 fi
}

makedatasets()
{
i=$((RANDOM/37+1))

echo -e "\nThe seed for the generation of the 1000-replicate permutation data is: ${3}.\nThe generated datasets(100 datasets in total) start from No.${i} dataset." > impMsg.txt
j=$((i+100))
filenum=1
percent=5
rm -rf ${1}/permutation_data
mkdir ${1}/permutation_data  
 
until [ $i -eq $j ]
do
    k=$((i+1))
    startnum=$(cat outfile | grep -n "${2}" | cut -d ":" -f 1 | sed -n "${i}p")
    end=$(cat outfile | grep -n "${2}" | cut -d ":" -f 1 | sed -n "${k}p")
    endnum=$((end-1))

    cat ${1}/outfile | sed -n "${startnum},${endnum}p" > ${1}/permutation_data/d${filenum}.phy
    i=$((i+1))
    filenum=$((filenum+1))
done
}

pl_partition()
{
 seqpos1=1
 seqt=0
 for file in ${1}/${2}/*phy
 do 
   name=$(echo "${file}" | cut -d "." -f 1 | cut -d "-" -f 2)
   seq=$(cat ${file} | cut -d " " -f 3)
   seqt=$(echo "$((seqt+seq))")
   echo "${name} = $seqpos1-$seqt;" >> ${1}/Analysis/permutationGTRI_analysis/partition
   seqpos1=$(echo "$((seqt +1))")
 done
}

partSort()
{
 op=$(echo "${1}_output")
   
 lnum=$(cat partition | wc -l)
 pi=1
 while [ "${pi}" -le "${lnum}" ]
 do
   gnname[$pi]=$(cat partition | sed -n "${pi}p" | cut -d " " -f 1)
   part[$pi]=$(cat partition | sed -n "${pi}p" | cut -d " " -f 3 | cut -d ";" -f 1)
   ((pi++))
 done

 cli=1
 while [ "${cli}" -le "${lnum}" ]
 do
   cl[$cli]=$(cat ${1} | grep "${gnname[$cli]}" | cut -d " " -f 2)
   ((cli++))
 done

 rm -rf ${op}
 i=1
 k=1
 while [ "${i}" -le "${lnum}" ]
 do
   if [ "${i}" -eq "${lnum}" ]; then
      j=${i}
   else
      j=$((i+1))
   fi

   if [ "${cl[$i]}" != "NA" ]; then
      partcon=$(echo "${part[$i]}")
      while [ "${j}" -le "${lnum}" ]
      do
        if [ "${cl[$j]}" != "NA" ]; then
           if [ "${cl[$j]}" -eq "${cl[$i]}" ]; then
              if [ "${i}" -lt "${lnum}" ]; then
                 partcon=$(echo "${partcon},${part[$j]}")
                 cl[$j]="NA"
              fi
           fi
        fi

        ((j++))
      done  
      if [ "${k}" -lt 10  ]; then
         echo "cluster0${k} = ${partcon}" >> ${op}
      else
         echo "cluster${k} = ${partcon}" >> ${op}
      fi
      ((k++))
   fi
   ((i++))
 done
}

pickbestrun()
{
  iterindex=1
  while [ "${iterindex}" -le "${2}" ]
  do
    loglkh[${iterindex}]=$(cat ${1}_run${2}/*iqtree | grep "Log-likelihood of the tree:" | cut -d " " -f 5)
    loglkhmd[${iterindex}]=$(echo "${loglkh[${iterindex}]} * 10000" | bc | cut -d "." -f 1 | cut -d "-" -f 2)
    ((iterindex++))
  done
  
  iterindex=2; loglkht=${loglkhmd[1]}; tnum=1
  while [ "${iterindex}" -le "${2}" ]
  do
    if [ "${loglkht}" -gt "${loglkhmd[${iterindex}]}" ]; then
       loglkht=${loglkhmd[${iterindex}]}
       tnum=${iterindex}
    fi
    ((iterindex++))
  done
}

logLXtract()
{
     seqname=${1}
     seqnamemd=$(echo "${1}_finish")

  clfiup=0 #number of cluster file                         #Obtain the total cluster number
  for clfile in ${2}/clinfo/*
  do
    fex=$(echo "${clfile}" | cut -d "." -f 2)
    if [ "${fex}" = "cl" ]; then
       ((clfiup++))
       kclnum[$clfiup]=$(basename ${clfile} | cut -d "." -f 1)
    fi
  done                                                     #Finish obtaining the total cluster number

  ki=0
  while [ "${ki}" -lt "${clfiup}" ]
  do
    ((ki++))
    clnumup=$(echo "${kclnum[$ki]}" | cut -d "k" -f 2)
    clnum=0
    ttllogLhmd=0
    while [ "${clnum}" -lt "${clnumup}" ]
    do
      ((clnum++))
      if [ "${clnum}" -lt 10 ]; then
         clnumfm=$(echo "0${clnum}")
      else
         clnumfm=${clnum}
      fi

      cd ${2}/${seqnamemd}/${kclnum[$ki]}/${seqname}_cluster${clnumfm}-out_bestrun
      logLh=$(cat *iqtree | grep "Log-likelihood of the tree:" | cut -d " " -f 5)
      logLhmd=$(echo "${logLh} * 10000" | bc | cut -d "." -f 1 | cut -d "-" -f 2)
      ttllogLhmd=$((ttllogLhmd + logLhmd))
    done
    echo -e "${seqname} k=${clnumup} -loglikelihood*10^4:\t${ttllogLhmd}" >> ${2}/${seqname}_logLhcollection
  done
  if [ ! -d "${2}/logLikelihood" ]; then mkdir ${2}/logLikelihood; fi
  mv ${2}/*logLhcollection ${2}/logLikelihood/   
}





#--------------------------------------------------------------------------------#
#Loading the flags
while getopts ":h:a:u:" opt
do
  case $opt in
               h) usage
                  exit;;
               a) alnf=$(echo "${OPTARG}" | cut -d "/" -f 1);;
   esac
done

if [ $OPTIND -eq 1 ]; then
   usage
   exit
fi


curdir=$(pwd)
resize -s 28 110 > /dev/null

read -p "
 Files checklist: 1.Alignment folder [-a]
                  2.get_cluster.R
                  3.logLkhforR_make.R
                  4.graphs.R
                  5.find_potent_knum.R
                  6.clusterAssign.R

    Dependencies: 1.IQ-Tree           www.iqtree.org
                  2.Newick Utilities  cegg.unige.ch/newick_utils
                  3.AMAS              github.com/marekborowiec/AMAS
                  4.fseqboot          emboss.open-bio.org/rel/rel6/embassy/phylipnew/fseqboot.html
                    (from Embassy-phylip,require Emboss).
                  5.GTP
                  6.R-base            www.r-project.org

      R packages: 1.distory           cran.r-project.org/web/packages/distory/index.html
                  2.phangorn          cran.r-project.org/web/packages/phangorn/index.html
                    (require rgl, igraph>=1.0,ape)
                  3.heatmap3          cran.r-project.org/web/packages/heatmap3/index.html
                  4.gplots            cran.r-project.org/web/packages/gplots/index.html
                  5.ggplot2           cran.r-project.org/web/packages/ggplot2/index.html
                  6.tidyr             cran.r-project.org/web/packages/tidyr/index.html
                  7.dplyr             cran.r-project.org/web/packages/dplyr/index.html
                  8.svglite           cran.r-project.org/web/packages/svglite/index.html

Continue?(Y/n): " answer

if [ "${answer}" == "Y" -o "${answer}" == "y" ]; then
   if [ -z "${alnf}" ] || [ ! -d "${curdir}/${alnf}" ] || [ ! -e get_cluster.R ] || [ ! -e logLkhforR_make.R ] || [ ! -e graphs.R ] || [ ! -e find_potent_knum.R ] || [ ! -e clusterAssign.R ]; then echo -e "\nMsg: Some files are missing.\n" && usage && exit; fi  #Check if all the files exist 

#read -p $'\nSpecify the directory of AMAS.py (Default is $HOME/Software/amas) ' amasdir
   amasdir="$HOME/Software/amas"
#read -p $'\nSpecify the directory of gtp.jar (Default is $HOME/Software/gtp) ' gtpdir
   gtpdir="$HOME/Software/gtp"

tstamp=$(date)
echo "
   Analysis started at ${tstamp}
   ---------------------------------------------------------------"

all="y"

#1.Building gene trees with bootstrapping values
#*************************************************************
echo "
         Part 1. Gene trees construction on the orignal data"

fdname="iqMF-autoWZgn_1000UFB"
if [ -d "${curdir}/${fdname}" ]; then rm -rf ${curdir}/${fdname}; fi
mkdir ${curdir}/${fdname} ${curdir}/${fdname}/report ${curdir}/${fdname}/trivia ${curdir}/${fdname}/bs-file

bsset="-bb 1000 -bnni -wbtl" #Perform Ultra Fast Bootstrapping with 1000 replicates.
for alnfile in ${curdir}/${alnf}/*phy
do
  IQTREE ${curdir} ${alnfile} ${bsset} ${fdname}
done

#Root the tree and extract the rooted trees into a file.

stepind=1
aftertreebuild ${curdir} ${fdname} ${stepind}

cp -r ${curdir}/${fdname}/bs-file ${curdir}/.bs-file
rm -rf Trees
mkdir Trees
mv ${curdir}/${fdname}/ ${curdir}/${fdname}_rooted/ ${curdir}/Trees/

#Reroot the bootstrapped trees
rm -rf bs-file_rooted
mkdir bs-file_rooted
for ufbfile in ${curdir}/.bs-file/*
do
  namepre=$(echo "${ufbfile}")
  name=$(basename ${namepre})
  nw_reroot ${ufbfile} ZF > ${curdir}/bs-file_rooted/${name}
done

rm -rf ${curdir}/.bs-file/

tstamp=$(date)
echo "
         Finished at ${tstamp}!
   ---------------------------------------------------------------"
#*************************************************************

#2.Generate permutation data
#*************************************************************
echo "
         Part 2. Generate permutation data"

python ${amasdir}/AMAS.py concat -i ${alnf}/*phy -f phylip -d dna -u phylip -t ${curdir}/AutoWZgn_concatenation.phy > /dev/null
rm -rf ${curdir}/partitions.txt
sed -i "s/ZF/ZF    /" AutoWZgn_concatenation.phy
sed -i "s/MW/MW    /" AutoWZgn_concatenation.phy
sed -i "s/CW/CW    /" AutoWZgn_concatenation.phy
sed -i "s/BC/BC    /" AutoWZgn_concatenation.phy
sed -i "s/GRW/GRW    /" AutoWZgn_concatenation.phy
sed -i "s/CRW/CRW    /" AutoWZgn_concatenation.phy

seednum=2
until [ $((seednum % 2)) -gt 0 ]
do
  seednum=${RANDOM}
done
fseqboot -sequence ${curdir}/AutoWZgn_concatenation.phy -outfile outfile -test o -seqtype d -rewriteformat p -reps 1000 -seed ${seednum} > /dev/null

taxanum=6
makedatasets ${curdir} ${taxanum} ${seednum}
rm -rf ${curdir}/outfile

tstamp=$(date)
echo "
         Finished at ${tstamp}!
   ---------------------------------------------------------------"
#*************************************************************

#3.Prepare files for permutation
#*************************************************************
echo "
         Part 3. File preparation for permutation analysis"

#Generate gene trees collection
fdname="iqMF-autoWZgn_br"
if [ -d "${curdir}/${fdname}" ]; then rm -rf ${curdir}/${fdname}; fi
mkdir ${curdir}/${fdname} ${curdir}/${fdname}/report ${curdir}/${fdname}/trivia

bsset="NoBS" #Don't need nodes support method.
for alnfile in ${curdir}/${alnf}/*phy
do
  IQTREE ${curdir} ${alnfile} ${bsset} ${fdname}
done

##Root the tree and extract the rooted trees into a file.
stepind=5
aftertreebuild ${curdir} ${fdname} ${stepind}

mv ${curdir}/${fdname}/ ${curdir}/${fdname}_rooted/ ${curdir}/Trees/  


#Get clusters information
rm -rf ${curdir}/Analysis
mkdir ${curdir}/Analysis
mv ${curdir}/gnTrees_collection.tre ${curdir}/impMsg.txt ${curdir}/get_cluster.R ${curdir}/Analysis/

cd ${curdir}/Analysis
Rscript get_cluster.R 2>&1 >/dev/null
rm -rf get_cluster.R

sexyn="n"
sexi=1
lnum=$(cat ${curdir}/Analysis/dend.order | wc -l)
i=1
while [ "${i}" -le "${lnum}" ]
do
  name=$(cat ${curdir}/Analysis/dend.order | sed -n "${i}p")
  nsmark=$(echo "${name}" | cut -c 1-2)
  if [ "${nsmark}" == "ns" ]; then 
     prefix=$(echo "${name}" | cut -d "_" -f 1)
     suffix=$(echo "${name}" | cut -d "_" -f 2)
     if [ "${suffix}" == "W" ] || [ "${suffix}" == "Z" ]; then 
        sexyn="y"
     fi
     if [ "${suffix}" == "W" ]; then
        if [ "${sexi}" -lt 10 ]; then seximd=$(echo "0${sexi}"); else seximd=${sexi}; fi
        echo "${name}" >> ${curdir}/Analysis/${seximd}-${prefix}.sex
        echo "${prefix}_Z" >> ${curdir}/Analysis/${seximd}-${prefix}.sex
        ((sexi++))
     fi
  else
     echo "${name}" >> ${curdir}/Analysis/A.auto
  fi
  ((i++))
done
rm -rf dend.order
if [ "${sexyn}" == "y" ]; then
   cat ${curdir}/Analysis/A.auto ${curdir}/Analysis/*sex > ${curdir}/Analysis/order
else
   mv ${curdir}/Analysis/A.auto ${curdir}/Analysis/order
fi
rm -rf ${curdir}/Analysis/A.auto ${curdir}/Analysis/*sex

rm -rf permutationGTRI_analysis
mkdir permutationGTRI_analysis
if [ ! -d clinfo ]; then mkdir clinfo; fi
mv *cl order clinfo/
mv clinfo/ permutationGTRI_analysis/

#Get partition file
pl_partition ${curdir} ${alnf}

#Manoevour the sequence files
cp -p ${curdir}/AutoWZgn_concatenation.phy ${curdir}/Analysis/permutationGTRI_analysis/
mv ${curdir}/permutation_data/*phy ${curdir}/Analysis/permutationGTRI_analysis/
rmdir ${curdir}/permutation_data

tstamp=$(date)
echo "
         Finished at ${tstamp}!
   ---------------------------------------------------------------"
#*************************************************************

#4.Permutation analysis
#*************************************************************
echo "
         Part 4. Permutation analysis"

cd ${curdir}/Analysis/permutationGTRI_analysis/

#Create partition-sorted file
clfi=0
rm -rf cl_sort
mkdir cl_sort
for file in clinfo/*cl
do  
  partSort ${file}
  filename=$(basename ${file})
  ((clfi++))
  clf[$clfi]=$(echo "${filename}_output")
  mv ./clinfo/*cl_output ./cl_sort/
done

#Execute the sequence splitting
ni=0
for file in *phy
do       
  ((ni++))
  name[$ni]=$(echo "${file}" | cut -d "." -f 1)
  mkdir ${name[$ni]}

  exei=1   
  while [ "${exei}" -le "${clfi}" ]     
  do
    python ${amasdir}/AMAS.py split -f phylip-int -d dna -i ${file} -l ./cl_sort/${clf[$exei]}  -u phylip > /dev/null
    kclnum=$(echo "${clf[$exei]}" | cut -d "." -f 1)
    mkdir ./${name[$ni]}/${kclnum}
    mv ./*cluster* ./${name[$ni]}/${kclnum}/
    ((exei++))
  done
done

#read -p $'\nContinue to Tree Construction? (y/n): ' totree
totree="y"
if [ "${totree}" = "y" ] || [ "${totree}" = "Y" ]; then
   
   #read -p $'\nSpecify the number of runs: ' runtime
   #if [ -z "${runtime}" ]; then runtime=2; fi
   runtime=1

   predir=$(pwd)
   tri=1
   while [ "${tri}" -le "${ni}" ]
   do
     trj=1
     while [ "${trj}" -le "${clfi}" ] 
     do    
       kclnum=$(echo "${clf[$trj]}" | cut -d "." -f 1)
       cd ${predir}/${name[$tri]}/${kclnum}/
       clnum=1
       for seqfile in *
       do         

         seqlen=$(sed -n "1p" ${seqfile} | cut -d " " -f 2)
         if [ "${seqlen}" -lt 10000 ]; then 
            nth="1"
         elif [ "${seqlen}" -lt 1000000 ]; then 
            nth="2"
         elif [ "${seqlen}" -lt 10000000 ]; then
            nth="3"
         else
            nth="4"
         fi

         if [ "${clnum}" -lt 10 ]; then
            teststring1=$(echo "cluster0${clnum}")
         else
            teststring1=$(echo "cluster${clnum}")
         fi

         teststring2=$(echo "${seqfile}" | grep "${teststring1}")
         if [ -n "${teststring2}" ]; then

            ri=0
            while [ "${ri}" -lt "${runtime}" ]
            do
              ((ri++))                  
              seqfname=$(echo "${seqfile}" | cut -d "." -f 1)
              iqtree -s ${seqfile} -m GTR+F+I -AICc -ninit 200 -ntop 50 -nt ${nth} -keep-ident -pre ${seqfname}.run > /dev/null
              ((clnum++))

              mkdir ${seqfname}_run${ri}

              for mvfile in *  #File sorting
              do
                mvfex=$(echo "${mvfile}" | cut -d "." -f 2)
                if [ "${mvfex}" = "run" ]; then
                   mv ${mvfile} ./${seqfname}_run${ri}/
                fi
              done             #File sorting finish
            done 

            rm -rf ${seqfile}  #For saving the storage room
            pickbestrun ${seqfname} ${ri} #Select the run with the highest log likelihood
            mkdir ${seqfname}_bestrun
            mv ${seqfname}_run${tnum}/* ${seqfname}_bestrun/ 

            di=1
            while [ "${di}" -le "${ri}" ]
            do
              rm -rf ${seqfname}_run${di}
              ((di++))
            done                          #Finish select the run with the highest log likelihood

         fi
       done
       ((trj++))
     done
 
     cd ${predir}
     mv ${name[$tri]} ${name[$tri]}_finish
     logLXtract ${name[$tri]} ${predir}
     ((tri++))
   done
   
   rm -rf ${predir}/*finish ${predir}/*phy

   #Extract log likelihood values from numerous files:
   mkdir ${predir}/value
   cd ${predir}/logLikelihood
   for file in *
   do
     name=$(echo "${file}")
     cat ${file} | cut -d $'\t' -f 2 > ${predir}/value/${name}
   done

   mv ${predir}/value/AutoWZgn_concatenation_logLhcollection ${predir}/value/AutoWZgn_concatenation_logLkhcollection

   #Generate logLkhforR.csv file
   mv ${curdir}/logLkhforR_make.R ${predir}/value/
   cd ${predir}/value
   Rscript logLkhforR_make.R 2>&1 >/dev/null
   rm -rf logLkhforR_make.R

   #File sorting
   mv ${predir}/value ${predir}/logLikelihood
   mv ${predir}/logLikelihood/value/logLkhcollection.csv ${predir}/logLikelihood/value/logLkhforR.csv ${predir}/
   rm -rf ${predir}/logLikelihood

   #Plot graph
   mv ${curdir}/graphs.R ${predir}/
   cd ${predir}
   Rscript graphs.R 2>&1 >/dev/null
   rm -rf graphs.R

   #File sorting
   rm -rf ${predir}/a.csv ${predir}/amd.csv ${predir}/b.csv ${predir}/bmd.csv ${predir}/cl_sort
fi

tstamp=$(date)
echo "
         Finished at ${tstamp}!
   ---------------------------------------------------------------"
#*************************************************************

#5.Post permutation analysis operation
#*************************************************************
echo "
         Part 5. Post permutation analysis operation"

#read -p "continue? " yn
yn="y"
if [ "${yn}" == "y" ]; then


if [ -e ${predir}/logLkhforR.csv ]; then
#Find potential number of clusters
   mv ${curdir}/find_potent_knum.R ${curdir}/AutoWZgn_concatenation.phy ${predir}/

   Rscript find_potent_knum.R 2>&1 >/dev/null
   rm -rf find_potent_knum.R
   pknumall=$(cat potent_knum)
   pknum=$(sed -n "1p" potent_knum)
   #read -p $"   
#Potential number of clusters are
#${pknumall}
#Continue with ${pknum} clusters? (y/n): " pkyn
   pkyn="y"
   if [ -z "${pkyn}" ] || [ "${pkyn}" == "Y" ] || [ "${pkyn}" == "y" ]; then 
      echo "
     ********************************************************
         Msg:The analysis continues with ${pknum} clusters
     ********************************************************"
   elif [ "${pkyn}" == "n" ] || [ "${pkyn}" == "N" ]; then
      read -p "Reset the number of clusters to: " pknumreset
      if [ -z "${pknumreset}" ] || [ "${pknumreset}" -eq "${pknumreset}" ] 1>/dev/null || [ "${pknumreset}" -gt "${clfi}" ]; then
         echo "
     *************************************************************************
        Msg:FAILED to reset...the analysis continues with ${pknum} clusters
     *************************************************************************"
      else
         pknum=${pknumreset}
         echo "
     ******************************************************************************************
        Msg:Number of clusters has been reset to ${pknum} and the analysis hence continues...
     ******************************************************************************************"
      fi
   else
      echo "
     *************************************************************************
        Msg:FAILED to reset...the analysis continues with ${pknum} clusters
     *************************************************************************"
   fi

   if [ -n ${pknum} ]; then
      #Create partition-sorted file
      rm -rf cl_sort
      mkdir cl_sort
      if [ "${pknum}" -gt 10 ]; then
         pknummd=$(echo "k${pknum}.cl")
      else
         pknummd=$(echo "k0${pknum}.cl")
      fi

      clfi=0
      for file in clinfo/${pknummd}
      do  
        partSort ${file}
        filename=$(basename ${file})
        ((clfi++))
        clf[$clfi]=$(echo "${filename}_output")
        mv ./clinfo/*cl_output ./cl_sort/
      done



      #Execute the sequence splitting
      ni=0
      for file in *phy
      do       
        ((ni++))
        name[$ni]=$(echo "${file}" | cut -d "." -f 1)
        mkdir ${name[$ni]}

        exei=1   
        while [ "${exei}" -le "${clfi}" ]     
        do
          python ${amasdir}/AMAS.py split -f phylip-int -d dna -i ${file} -l ./cl_sort/${clf[$exei]}  -u phylip > /dev/null
          kclnum=$(echo "${clf[$exei]}" | cut -d "." -f 1)
          mkdir ./${name[$ni]}/${kclnum}
          mv ./*cluster* ./${name[$ni]}/${kclnum}/
          ((exei++))
        done
      done

      #read -p $'\nContinue to Tree Construction? (y/n): ' totree
      totree="y"
      if [ "${totree}" = "y" ] || [ "${totree}" = "Y" ]; then
   
         #read -p $'\nSpecify the number of runs: ' runtime
         #if [ -z "${runtime}" ]; then runtime=2; fi
         runtime=2

         rm -rf ${curdir}/supergn_tree/
         mkdir ${curdir}/supergn_tree/
         tri=1
         while [ "${tri}" -le "${ni}" ]
         do
           trj=1
           while [ "${trj}" -le "${clfi}" ] 
           do    
             kclnum=$(echo "${clf[$trj]}" | cut -d "." -f 1)
             cd ${predir}/${name[$tri]}/${kclnum}/
             clnum=1
             for seqfile in *
             do         
               seqlen=$(sed -n "1p" ${seqfile} | cut -d " " -f 2)
               if [ "${seqlen}" -lt 10000 ]; then 
                  nth="1"
               elif [ "${seqlen}" -lt 1000000 ]; then 
                  nth="2"
               elif [ "${seqlen}" -lt 10000000 ]; then
                  nth="3"
               else
                  nth="4"
               fi

               if [ "${clnum}" -lt 10 ]; then
                  teststring1=$(echo "cluster0${clnum}")
               else
                  teststring1=$(echo "cluster${clnum}")
               fi

               teststring2=$(echo "${seqfile}" | grep "${teststring1}")
               if [ -n "${teststring2}" ]; then

                  ri=0
                  while [ "${ri}" -lt "${runtime}" ]
                  do
                    ((ri++))
                    if [ "${ri}" -gt 1 ]; then clnum=$((clnum-1)); fi                   
                    seqfname=$(echo "${seqfile}" | cut -d "." -f 1)
                    iqtree -s ${seqfile} -m TEST -AICc -ninit 200 -ntop 50 -nt ${nth} -keep-ident -pre ${seqfname}.run > /dev/null
                    ((clnum++))

                    mkdir ${seqfname}_run${ri}

                    for mvfile in *  #File sorting
                    do
                      mvfex=$(echo "${mvfile}" | cut -d "." -f 2)
                      if [ "${mvfex}" = "run" ]; then
                         mv ${mvfile} ./${seqfname}_run${ri}/
                      fi
                    done             #File sorting finish
                  done 

                  rm -rf ${seqfile}  #For saving the storage room
                  pickbestrun ${seqfname} ${ri} #Select the run with the highest log likelihood
                  mkdir ${seqfname}_bestrun
                  mv ${seqfname}_run${tnum}/* ${seqfname}_bestrun/ 
                  mv ${seqfname}_bestrun/*treefile ${curdir}/supergn_tree/

                  di=1
                  while [ "${di}" -le "${ri}" ]
                  do
                    rm -rf ${seqfname}_run${di}
                   ((di++))
                  done                          #Finish select the run with the highest log likelihood
               
               fi
             done
             ((trj++))
           done
 
           cd ${predir}
           mv ${name[$tri]} ${name[$tri]}_finish
          ((tri++))
        done
      
        rm -rf ${predir}/*phy
      
        stepind=9
        #Extract the supergn_tree
        aftertreebuild ${curdir} supergn_tree ${stepind}
        mv ${curdir}/supergn_tree* ${curdir}/Trees/
      fi
   fi
fi

tstamp=$(date)
echo "
         Finished at ${tstamp}!
   ---------------------------------------------------------------"

fi
#*************************************************************

#6.Cluster Assignment
#*************************************************************
echo "
         Part 6. Cluster assignment"

if [ -e ${curdir}/gnTrees_collection.tre ]; then
   #Distribute the supergene trees
   cd ${curdir}
   rm -rf supergn_trees
   mkdir supergn_trees   
   spi=1
   while [ "${spi}" -le "${pknum}" ]
   do
     name=$(cat gnTrees_collection.tre | sed -n "${spi}p" | cut -d " " -f 1)
     nwt=$(cat gnTrees_collection.tre | sed -n "${spi}p" | cut -d " " -f 2)     
     echo "${nwt}" > ./supergn_trees/${name}.tre
     ((spi++))
   done

   #Generate the tree collection file to calculate Geodesic distance
   rm -rf tree.Geodesic
   mkdir tree.Geodesic
   for sgt in ${curdir}/supergn_trees/*
   do
     sgtname=$(basename ${sgt} | cut -d "." -f 1)
     for bst in ${curdir}/bs-file_rooted/*
     do
       bstname=$(basename ${bst} | cut -d "." -f 1 | cut -d "-" -f 2)
       cat ${sgt} ${bst} > ${curdir}/tree.Geodesic/${bstname}-${sgtname}.tre
     done
   done
   
   #Calculate the Geodesic distances between each supergene tree and the BS trees 
   cd ${curdir}/tree.Geodesic/
   for tgfile in *tre
   do
     name=$(echo "${tgfile}" | cut -d "." -f 1)
     java -jar ${gtpdir}/gtp.jar -d -o ${name}.gd -r 0 ${tgfile}
     cat ${name}.gd | cut -d $'\t' -f 3 > ${name}.geod    
     clsort=$(echo "${tgfile}" | cut -d "-" -f 1)
     if [ ! -d "${curdir}/tree.Geodesic/${clsort}" ]; then mkdir ${curdir}/tree.Geodesic/${clsort}; fi
     mv ./*geod ./${clsort}/
     rm -rf ${tgfile} ${name}.gd
   done

   #Assign clusters for each BS tree
   for alnfile in ${curdir}/${alnf}/*phy
   do
     alnname=$(basename ${alnfile} | cut -d "." -f 1 | cut -d "-" -f 2)
     fnum=$(ls ${curdir}/tree.Geodesic/${alnname}/ | wc -l)     
     bsrep=$(cat ${curdir}/tree.Geodesic/${alnname}/${alnname}-cluster01.geod | wc -l)
     i=1
     while [ "${i}" -le "${bsrep}" ]
     do
       j=1
       while [ "${j}" -le "${fnum}" ]
       do
         if [ "${j}" -lt 10 ]; then jmd=$(echo "cluster0${j}"); else jmd=$(echo "cluster${j}"); fi       
         gdvpre=$(cat ${curdir}/tree.Geodesic/${alnname}/${alnname}-${jmd}.geod | sed -n "${i}p")
         gdv[$j]=$(echo "${gdvpre} * 100000000" | bc | cut -d "." -f 1)
         ((j++))
       done

       j=1
       gdvmin=${gdv[1]}
       min=1
       while [ "${j}" -le "${fnum}" ]
       do
         if [ "${gdv[$j]}" -lt "${gdvmin}" ]; then
            gdvmin=${gdv[$j]}
            min=${j}
         fi
         ((j++))
       done
       
       echo "${min}" >> ${curdir}/tree.Geodesic/${alnname}-cluster.assign
 
       ((i++))
     done
   done

   #Summary the information
   rm -rf ${curdir}/tree.Geodesic/Cluster_Assignment
   mkdir ${curdir}/tree.Geodesic/Cluster_Assignment
   
   rm -rf ${curdir}/tree.Geodesic/Cluster_Assignment.summary ${curdir}/tree.Geodesic/Cluster_Assignment_shortdata.summary

   odline=$(cat ${predir}/clinfo/order | wc -l)
   odi=1
   while [ "${odi}" -le "${odline}" ]
   do
     gnord=$(cat ${predir}/clinfo/order | sed -n "${odi}p")
     gdasn=$(echo "${curdir}/tree.Geodesic/${gnord}-cluster.assign")   
     clori=${odi}
     #clori=$(cat ${predir}/clinfo/${pknummd} | grep "${gnord}" | cut -d " " -f 2)
     if [ -z "${clori}" ]; then clori=${fnum}; fi
     if [ "${clori}" -lt 10 ]; then clori=$(echo "0${clori}"); fi
     name=$(echo -e "${clori}\t${gnord}")

     k=1
     while [ "${k}" -le "${fnum}" ]
     do
       kcount[$k]=0
       while read line
       do
         if [ "${line}" -eq "${k}" ]; then
            kcount[$k]=$((kcount[$k] + 1))
         fi
       done < ${gdasn}

       if [ -z "${kcount[$k]}" ]; then kcount[$k]=0; fi       

       ((k++))
     done

     #Produce short-data table output
     l=2
     ksum=${kcount[1]}
     while [ "${l}" -le "${fnum}" ]
     do
       ksum=$(echo -e "${ksum}\t${kcount[$l]}")
       ((l++))
     done
     echo -e "${name}\t${ksum}" >> ${curdir}/tree.Geodesic/Cluster_Assignment_shortdata.summary        

     #Produce long-data table output
     l=1
     while [ "${l}" -le "${fnum}" ]
     do
       echo -e "${name}\tcluster${l}\t${kcount[$l]}" >> ${curdir}/tree.Geodesic/Cluster_Assignment.summary
       ((l++))
     done
     
     #File sorting
     mv ${curdir}/tree.Geodesic/${gnord}-cluster.assign ${curdir}/tree.Geodesic/Cluster_Assignment

     ((odi++))
   done

   #File sorting   
   rm -rf ${curdir}/tree.Geodesic/Geodesic_distance
   mkdir ${curdir}/tree.Geodesic/Geodesic_distance
   for alnfile in ${curdir}/${alnf}/*phy
   do
     alnname=$(basename ${alnfile} | cut -d "." -f 1 | cut -d "-" -f 2)
     mv ${curdir}/tree.Geodesic/${alnname} ${curdir}/tree.Geodesic/Geodesic_distance/
   done
   if [ -d "${curdir}/untgn_aln" ]; then
      for alnfile in ${curdir}/untgn_aln/*phy
      do 
        alnname=$(basename ${alnfile} | cut -d "." -f 1 | cut -d "-" -f 2)
        mv ${curdir}/tree.Geodesic/${alnname} ${curdir}/tree.Geodesic/Geodesic_distance/
      done
   fi
   
   #Plot graphs for cluster assignment
   mv ${curdir}/clusterAssign.R ${curdir}/tree.Geodesic/
   cd ${curdir}/tree.Geodesic/
   Rscript clusterAssign.R 2>&1 >/dev/null
   rm -rf clusterAssign.R
   
fi

tstamp=$(date)
echo "
         Finished at ${tstamp}!
   ---------------------------------------------------------------"
#*************************************************************
rm -rf ${curdir}/gnTrees_collection.tre ${curdir}/supergn_trees 
mv ${curdir}/bs-file_rooted ${curdir}/Trees/

fi
