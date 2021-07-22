#!/bin/bash

while getopts c:d:x:t: flag
do
    case "${flag}" in
        c) contrasts=${OPTARG};;
        d) directory=${OPTARG};;
        x) cutsequence="${OPTARG}";;
        t) transcriptome=${OPTARG};;
    esac
done


# directory=/stor/work/Wilke/cdupai/Cory_RNAseq/USE_ME
setup=${directory}/setup_files
settemp=${setup}/temp
fastqs=${directory}/fastqs
remseq=${setup}/remove_sequences.fasta
# echo ${contrasts}
# echo $directory
# echo $cutsequence
# echo $transcriptome
contrasts=$setup/contrasts.tsv
# transcriptome=$setup/remove_sequences.fasta
# transcriptome="GENERATE"
# cutsequence="NO"
echo $directory
echo $cutsequence


firstread=()
secondread=()
paired_reads="FALSE"
combine="FALSE"

mkdir -p $settemp

twoconts=${setup}/paired_contrasts.tsv

rm -rf $twoconts

# Pair reads
  while IFS=$'\t' read -r name2 fp2 cond2 rep2 rp2 sb2
  do
    if [ $sb2 -gt 1 ];then
      combine="TRUE"
    fi
    
    if [ $rp2 -eq 2 ];then
      paired_reads="TRUE"
        while IFS=$'\t' read -r name1 fp1 cond1 rep1 rp1 sb1
        do
          if [ $cond1 == $cond2 ] && [ $rep1 -eq $rep2 ] && [ $sb1 -eq  $sb2 ] && [ $rp1 -eq 1 ];then
              firstread+=("$fp1")
              secondread+=("$fp2")
              # echo -e "${name1}\t${name2}\t${fp1}\t${fp2}\t${cond1}\t${rep1}" >> $twoconts
          fi
        done < $contrasts
    fi
  done < $contrasts

#Flexbar
if [ "$cutsequence" != "Continue" ];then
  echo "HERE"
  fpid=()
  trimdir=${fastqs}/trimmed
  mkdir -p $trimdir
  flexcont=$settemp/flex_contrasts.tsv
  rm -rf $flexcont

    while IFS=$'\t' read -r name fp cond rep rp sb
    do
      sing_flexbar="YES"
        # varout=$trimdir/${fp##*/}
        varout=${trimdir}/${name}
        # $trimdir/${cond}_${rep}_${rp}.fastq
        if [ $paired_reads == "TRUE" ];then
          for i in "${!firstread[@]}"; do
          # echo ${firstread[$i]}
            if [ ${firstread[$i]} == $fp ];then
            fp2=${secondread[$i]}
              varout2=${fp2##*/}
              varout2=${varout2%%.*fastq*}

              flexbar --reads $fp --reads2 ${secondread[$i]} --adapters $remseq \
              --adapter-trim-end ANY --target $varout &
              fpid+=($!)
              echo -e "${name}\t${varout}_1.fastq\t${cond}\t${rep}\t${rp}\t${sb}" >> $flexcont
              echo -e "${varout2}\t${varout}_2.fastq\t${cond}\t${rep}\t2\t${sb}" >> $flexcont

              sing_flexbar="NO"
            elif [ ${secondread[$i]} == $fp ];then
              sing_flexbar="NO"
            fi
          done
        fi
        # echo $sing_flexbar
        # echo $varout
        if [ $sing_flexbar == "YES" ];then
          flexbar --reads $fp --adapters $remseq \
              --adapter-trim-end ANY --target $varout &
              echo $fpid
              
              echo -e "${name}\t${varout}\t${cond}\t${rep}\t${rp}\t${sb}" >> $flexcont
        fi
    done < $contrasts
    contrasts=$flexcont
    echo $contrasts
    for f in ${fpid[*]}; do
      echo $f
      wait $f
    done
fi



#Combine files (LOOK HERE LATER)
if [ $combine == "TRUE" ];then
echo "Combine!"
  firstread=()
  secondread=()
  tempfastq=${fastqs}/temp
  combcont=${settemp}/temp_contrasts.tsv
  combcont2=${settemp}/temp_contrasts_temp.tsv
  rm -rf $tempfastq
  mkdir -p $tempfastq
  rm -rf $combcont

  while IFS=$'\t' read -r name fp cond rep rp sb
  do
    varout=$tempfastq/${cond}_${rep}_R${rp}.fastq
    zcat -f "${fp}" >> $varout
    
    echo -e "${varout}\t${cond}\t${rep}\t${rp}" >> $combcont
    
    if [ $rp -eq 2 ];then
        varout2=$tempfastq/${cond}_${rep}_R1.fastq
        firstread+=("$varout2")
        secondread+=("$varout")
    fi
  done < $contrasts
  sort $combcont | uniq > $combcont2
  mv -f $combcont2 $combcont
  rm -rf $combcont2
  contrasts=$combcont
fi
# 
#     rm ${fastqs}/*
#     mv -f $newcont ${setup}/contrasts.tsv
#     mv -f ${tempdir}/* ${fastqs}/
#     rm -rf ${tempdir}

# Generate transcriptome
# if [ $transcriptome == "GENERATE" ];then
#   echo $transcriptome
#   concated=$settemp/combined.fq
#   rm -rf $concated
#   while IFS=$'\t' read -r name fp cond rep rp
#   do
#     zcat -f "${fp}" >> $concated
#   done < $contrasts
# 
#   awk '(NR%4==2)' ${settemp}/combined.fq > ${settemp}/combined.fa
#   sort ${settemp}/combined.fa > ${settemp}/sorted.txt
#   uniq -cd ${settemp}/sorted.txt > ${settemp}/uniq.txt
#   # grep -E -v "^[[:space:]]*[1-9][[:space:]]|^[[:space:]]*[1-4][0-9][[:space:]]|^[[:space:]]*50" \
#   grep -E -v "^[[:space:]]*[1-9][[:space:]]|^[[:space:]]*[1][0-5][[:space:]]" \
#   ${settemp}/uniq.txt > ${settemp}/counts.txt
#   
#   awk '{print $2}' ${settemp}/counts.txt | \
#   awk '$1=">"$1"\n"$1' > ${setup}/transcriptome.fasta
#   
#   transcriptome=${setup}/transcriptome.fasta
# 
# fi
#transcriptome=${setup}/transcriptome.fasta
firstread=($(echo "${firstread[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
secondread=($(echo "${secondread[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

# Kallisto
kd=${directory}/kallisto
kcons=${setup}/deseq_contrasts.tsv
rm -rf $kcons
kallisto index --index $setup/kallisto.index $transcriptome
mkdir -p $kd/logs
kpid=()
while IFS=$'\t' read -r cond fp rep rp
do
  # echo $fp $cond $rep $rp
  # echo $fp
  # echo "--------------"
  for i in "${!firstread[@]}"; 
  do
    sing_kallisto="YES"
     
    # echo $i
    # echo ${firstread[$i]}
    if [ "$fp" == "${firstread[$i]}" ];then
      fp2=${secondread[$i]}
      echo $fp
      # echo $fp2
      # 
      # echo "|||||||||||||||||"
      
      kallisto quant --index $setup/kallisto.index --output $kd/${cond}_${rep} \
      $fp $fp2 &> $kd/logs/${cond}_${rep}_logs.txt &
      kpid+=($!)

              
      echo -e "${kd}/${cond}_${rep}/abundance.h5\t${cond}\t${rep}" >> $kcons
      sing_kallisto="NO"
    elif [ "$rp" == "2" ];then
        sing_kallisto="NO"
    fi
  done
  if [ sing_kallisto == "YES" ];then
    kallisto quant --index $setup/kallisto.index --output ${kd}/${cond}_${rep} \
    $fp --single --sd 1 --fragment-length 150 &> $kd/logs/${cond}_${rep}_logs.txt &
    kpid+=($!)

    echo -e "${kd}/${cond}_${rep}/abundance.h5\t${cond}\t${rep}" >> $kcons
  fi
done< $contrasts

# for k in ${kpid[*]}; do
#     echo $k
#     wait $k
# done

# DESeq
# Rscript DESEQ