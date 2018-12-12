#!/bin/bash

#Functions
usage()
{
 echo "Usage: [-a] <Alignment folder>"
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

makedatasets()
{
#i=$((RANDOM/37+1))
i=1
#echo -e "The seed for the generation of the 1000-replicate permutation data is: ${3}.\nThe generated datasets(100 datasets in total) start from No.${i} dataset." > impMsg.txt
echo -e "The seed for the generation of the 1000-replicate permutation data is: ${3}." > impMsg.txt

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
   name=$(basename "${file}" | cut -d "." -f 1 | cut -d "-" -f 2)
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

logLXtract()
{
  clfiup=$(ls ${2}/*iqtree | wc -l)
  k=1
  while [ "${k}" -le "${clfiup}" ]
  do
    if [ "${k}" -lt 10 ]; then kmd=$(echo "0${k}"); else kmd=${k}; fi
    logLh=$(cat ${2}/*cluster${kmd}*iqtree | grep "Log-likelihood of the tree:" | cut -d " " -f 5)
    logLhmd=$(echo "${logLh} * 10000" | bc | cut -d "." -f 1 | cut -d "-" -f 2)
    ttllogLhmd=$((ttllogLhmd + logLhmd))
    ((k++))    
  done
  seqname=$(echo "${2}" | rev | cut -d "/" -f 2 | rev)
  echo -e "${seqname} k=${clfiup} -loglikelihood*10^4:\t${ttllogLhmd}" >> ${1}/${seqname}_logLhcollection
}

MinFinder()
{
 fnum=$(ls ${1} | wc -l)
 alnname=$(echo "${1}" | rev | cut -d "/" -f 1 | rev)
 bsreppre=$(cat ${1}/${alnname}-cluster01.geod | wc -l)
 bsrep=$((bsreppre - 1))
 i=1
 while [ "${i}" -le "${bsrep}" ]
 do
   j=1
   while [ "${j}" -le "${fnum}" ]
   do
      if [ "${j}" -lt 10 ]; then jmd=$(echo "cluster0${j}"); else jmd=$(echo "cluster${j}"); fi       
      geodvpre=$(cat ${1}/${alnname}-${jmd}.geod | sed -n "${i}p")
      geodv[$j]=$(echo "${geodvpre} * 100000000" | bc | cut -d "." -f 1)
      ((j++))
   done

   j=1
   geodvmin=${geodv[1]}
   min=1
   while [ "${j}" -le "${fnum}" ]
   do
     if [ "${geodv[$j]}" -lt "${geodvmin}" ]; then
        geodvmin=${geodv[$j]}
        min=${j}
     fi
     ((j++))
   done
       
   echo "${min}" >> ${2}/${alnname}-cluster.assign
 
   ((i++))
 done
}

gtpcal()
{
  name=$(echo "${3}" | cut -d "." -f 1)
  java -jar ${1}/gtp.jar -d -o ${name}.gd -r 0 ${2}/${3}
  cat ${name}.gd | cut -d $'\t' -f 3 > ${name}.geod
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
                  5.GTP               github.com/kgori/GTP
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
SECONDS=0

fdname="iqMF-autoWZgn_1000UFB"
if [ -d "${curdir}/${fdname}" ]; then rm -rf ${curdir}/${fdname}; fi
mkdir ${curdir}/${fdname} ${curdir}/${fdname}/report ${curdir}/${fdname}/trivia ${curdir}/${fdname}/bs-file

bsset="-bb 1000 -bnni -wbtl" #Perform Ultra Fast Bootstrapping with 1000 replicates.
cd ${curdir}/${alnf}/
parallel --no-notice -j 4 "iqtree -s {} -m TEST -AICc -ninit 200 -ntop 50 -nt 1 ${bsset}" ::: *.phy > /dev/null

#File sorting
mv ${curdir}/${alnf}/*ufboot ${curdir}/${fdname}/bs-file
mv ${curdir}/${alnf}/*iqtree ${curdir}/${fdname}/report/
mv ${curdir}/${alnf}/*treefile ${curdir}/${fdname}/
mv ${curdir}/${alnf}/*log ${curdir}/${alnf}/*gz ${curdir}/${alnf}/*bionj ${curdir}/${alnf}/*mldist ${curdir}/${alnf}/*contree ${curdir}/${alnf}/*splits.nex ${curdir}/${fdname}/trivia

#Extract the model information
for i in ${curdir}/${fdname}/report/*iqtree
do
  name=$(basename ${i} | cut -d "." -f 1-2)
  model=$(cat ${i} | grep "Best-fit model according to AICc:" | cut -d " " -f 6)
  echo -e "${name}\t${model}" >> ${curdir}/${alnf}/iq-modelset
done

#Root the tree and extract the rooted trees into a file.
stepind=1
aftertreebuild ${curdir} ${fdname} ${stepind}

cp -r ${curdir}/${fdname}/bs-file ${curdir}/

rm -rf ${curdir}/Trees
mkdir ${curdir}/Trees
mv ${curdir}/${fdname}/ ${curdir}/${fdname}_rooted/ ${curdir}/Trees/

#Reroot the bootstrapped trees
rm -rf ${curdir}/bs-file_rooted
mkdir ${curdir}/bs-file_rooted
for ufbfile in ${curdir}/bs-file/*
do
  namepre=$(echo "${ufbfile}")
  name=$(basename ${namepre})
  nw_reroot ${ufbfile} ZF > ${curdir}/bs-file_rooted/${name}
done

rm -rf ${curdir}/bs-file/
cd ${curdir}

duration=$SECONDS
echo "
         $(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
   ---------------------------------------------------------------"
#*************************************************************

#2.Generate permutation data
#*************************************************************
echo "
         Part 2. Generate permutation data"
SECONDS=0

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
fseqboot -sequence ${curdir}/AutoWZgn_concatenation.phy -outfile outfile -test o -seqtype d -rewriteformat p -reps 101 -seed ${seednum} > /dev/null

taxanum=6
makedatasets ${curdir} ${taxanum} ${seednum}
rm -rf ${curdir}/outfile

duration=$SECONDS
echo "
         $(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
   ---------------------------------------------------------------"
#*************************************************************

#3.Prepare files for permutation
#*************************************************************
echo "
         Part 3. File preparation for permutation analysis"
SECONDS=0

#Generate gene trees collection
fdname="iqMF-autoWZgn_br"
if [ -d "${curdir}/${fdname}" ]; then rm -rf ${curdir}/${fdname}; fi
mkdir ${curdir}/${fdname} ${curdir}/${fdname}/report ${curdir}/${fdname}/trivia

bsset="" #Don't need nodes support method.
cd ${curdir}/${alnf}/
parallel --no-notice -j 4 --colsep '\t' "iqtree -s {1} -m {2} -AICc -ninit 200 -ntop 50 -nt 1 ${bsset}" :::: iq-modelset > /dev/null
mv iq-modelset ${curdir}/

#File sorting
mv ${curdir}/${alnf}/*iqtree ${curdir}/${fdname}/report/
mv ${curdir}/${alnf}/*treefile ${curdir}/${fdname}/
mv ${curdir}/${alnf}/*log ${curdir}/${alnf}/*gz ${curdir}/${alnf}/*bionj ${curdir}/${alnf}/*mldist ${curdir}/${fdname}/trivia

##Root the tree and extract the rooted trees into a file.
stepind=5
aftertreebuild ${curdir} ${fdname} ${stepind}

for trefile in ${curdir}/${fdname}_rooted/*treefile
do
  name=$(basename ${trefile} | cut -d "." -f 1-2)  
  cat ${curdir}/bs-file_rooted/${name}.ufboot ${trefile} > ${curdir}/bs-file_rooted/${name}.ufboot-temp
  rm ${curdir}/bs-file_rooted/${name}.ufboot
  mv ${curdir}/bs-file_rooted/${name}.ufboot-temp ${curdir}/bs-file_rooted/${name}.ufboot
done

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

duration=$SECONDS
echo "
         $(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
   ---------------------------------------------------------------"
#*************************************************************

#4.Permutation analysis
#*************************************************************
echo "
         Part 4. Permutation analysis"
SECONDS=0

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
  rm -rf ${file}
done

#read -p $'\nContinue to Tree Construction? (y/n): ' totree
totree="y"
if [ "${totree}" = "y" ] || [ "${totree}" = "Y" ]; then

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
         directory=$(pwd)
         echo "${directory}/${seqfile}" >> ${curdir}/permtdata
         echo "${directory}" >> ${curdir}/permtdir-temp
       done
       ((trj++))
     done
     ((tri++))
   done
   
   cd ${curdir}
   parallel --no-notice -j 4 "iqtree -s {} -m GTR+F+I -AICc -ninit 200 -ntop 50 -nt 1" :::: permtdata > /dev/null
   rm -rf permtdata
   
   cat ${curdir}/permtdir-temp | sort -n | uniq > ${curdir}/permtdir
   rm -rf ${curdir}/permtdir-temp

   export -f logLXtract
   parallel --no-notice -j 4 logLXtract ::: ${predir} :::: permtdir > /dev/null
   rm -rf permtdir
   if [ ! -d "${predir}/logLikelihood" ]; then mkdir ${predir}/logLikelihood; fi
   mv ${predir}/*logLhcollection ${predir}/logLikelihood/ 
   rm -rf ${predir}/d*

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

duration=$SECONDS
echo "
         $(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
   ---------------------------------------------------------------"
#*************************************************************

#5.Post permutation analysis operation
#*************************************************************
echo "
         Part 5. Post permutation analysis operation"
SECONDS=0

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
   read -p $"   
Potential number of clusters are
${pknumall}
Continue with ${pknum} clusters? (y/n): " pkyn
   #pkyn="y"
   if [ -z "${pkyn}" ] || [ "${pkyn}" == "Y" ] || [ "${pkyn}" == "y" ]; then 
      echo "
     ********************************************************
         Msg:The analysis continues with ${pknum} clusters
     ********************************************************"
   elif [ "${pkyn}" == "n" ] || [ "${pkyn}" == "N" ]; then
      read -p "Reset the number of clusters to: " pknumreset
      if [ -n "${pknumreset}" ] && [ "${pknumreset}" -eq "${pknumreset}" ] 2>/dev/null && [ "${pknumreset}" -le "${clfi}" ]; then
         pknum=${pknumreset}
         echo "
     ******************************************************************************************
        Msg:Number of clusters has been reset to ${pknum} and the analysis hence continues...
     ******************************************************************************************"
      else
         echo "
     *************************************************************************
        Msg:FAILED to reset...the analysis continues with ${pknum} clusters
     *************************************************************************"
      fi
   else
      echo "
     *************************************************************************
        Msg:FAILED to reset...the analysis continues with ${pknum} clusters
     *************************************************************************"
   fi

   if [ -n ${pknum} ]; then

      if [ "${pknum}" -gt 10 ]; then
         pknummd=$(echo "k${pknum}")
      else
         pknummd=$(echo "k0${pknum}")
      fi

      rm -rf ${curdir}/supergn_tree
      mkdir ${curdir}/supergn_tree
      
      mv ${predir}/AutoWZgn_concatenation/${pknummd}/*treefile ${curdir}/supergn_tree

      stepind=9
      #Extract the supergn_tree
      aftertreebuild ${curdir} supergn_tree ${stepind}
      mv ${curdir}/supergn_tree* ${curdir}/Trees/
     
   fi
fi

duration=$SECONDS
echo "
         $(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
   ---------------------------------------------------------------"

fi
#*************************************************************

#6.Cluster Assignment
#*************************************************************
echo "
         Part 6. Cluster assignment"
SECONDS=0

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
       echo "${curdir}/tree.Geodesic/${bstname}" >> ${curdir}/tree.Geodesic/geodir-temp
       echo "${bstname}-${sgtname}.tre" >> ${curdir}/tree.Geodesic/treelist
     done
   done
   
   #Calculate the Geodesic distances between each supergene tree and the BS trees    
   cd ${curdir}/tree.Geodesic/
   export -f gtpcal
   parallel --no-notice -j 4  gtpcal ::: ${gtpdir} ::: ${curdir}/tree.Geodesic :::: treelist > /dev/null

   rm ${curdir}/tree.Geodesic/treelist
   rm -rf *tre *gd
   for kdfile in *geod
   do
     clsort=$(echo "${kdfile}" | cut -d "-" -f 1)
     if [ ! -d "${curdir}/tree.Geodesic/${clsort}" ]; then mkdir ${curdir}/tree.Geodesic/${clsort}; fi
     mv ${kdfile} ./${clsort}/
   done
     

   #Assign clusters for each BS tree
   export -f MinFinder
   cat ${curdir}/tree.Geodesic/geodir-temp | sort -n | uniq > ${curdir}/tree.Geodesic/geodir
   rm -rf ${curdir}/tree.Geodesic/geodir-temp
   cd ${curdir}/tree.Geodesic/
   parallel --no-notice -j 4 MinFinder :::: geodir ::: ${curdir}/tree.Geodesic > /dev/null

   #Summary the information
   rm -rf ${curdir}/tree.Geodesic/Cluster_Assignment
   mkdir ${curdir}/tree.Geodesic/Cluster_Assignment

   fnum=$(ls $(sed -n "1p" ${curdir}/tree.Geodesic/geodir) | wc -l)
   rm -rf ${curdir}/tree.Geodesic/geodir
   odline=$(cat ${predir}/clinfo/order | wc -l)
   odi=1
   while [ "${odi}" -le "${odline}" ]
   do
     gnord=$(cat ${predir}/clinfo/order | sed -n "${odi}p")
     geodasn=$(echo "${curdir}/tree.Geodesic/${gnord}-cluster.assign")   
     clori=${odi}
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
       done < ${geodasn}

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
   
   #Plot graphs for cluster assignment
   mv ${curdir}/clusterAssign.R ${curdir}/tree.Geodesic/
   cd ${curdir}/tree.Geodesic/
   Rscript clusterAssign.R ${pknum} 2>&1 >/dev/null
   rm -rf clusterAssign.R
   
fi

duration=$SECONDS
echo "
         $(($duration / 60)) minutes,$(($duration % 60)) seconds elapsed.
   ---------------------------------------------------------------"
#*************************************************************
rm -rf ${curdir}/supergn_trees 
mv ${curdir}/gnTrees_collection.tre ${curdir}/supergnTrees_collection.tre
mv ${curdir}/bs-file_rooted ${curdir}/supergnTrees_collection.tre ${curdir}/Trees/

tstamp=$(date)
echo "
   Analysis finished at ${tstamp}"
fi
