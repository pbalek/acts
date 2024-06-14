#!/bin/bash

for fileCSV in $@
    do
    if [ -f ${fileCSV} ]
        then
        line=`sed 1d $fileCSV | head -1`
        nums=`echo $line | sed -e 's/,/\n/g' | tail -n +7 | grep -v 'e' | grep -v '\.' `
        for num in `echo $nums`
        do
            echo "${num}"
            sed -i "s/${num},/${num}.,/g" ${fileCSV}
        done

        fileROOT=`basename $fileCSV | sed 's/.csv$/.root/'`
        root -l -q "csv2root.C(\"${fileCSV}\", \"${fileROOT}\")"
    else
        echo "File ${fileCSV} does not exist"
    fi
done
