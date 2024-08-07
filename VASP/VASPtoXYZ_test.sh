#!/bin/sh
# lipai@mail.ustc.edu.cn
# creat xyz files using OUTCAR only
if [ $# = 0 ]; then
        out="OUTCAR"
        name=""
else
        out=$1
        name=$1"-"
fi
echo "Input filename:"  $out
rm *.temp

# total number of ions in the system
num_atom=`grep -m 1 "NIONS =" $out|awk '{print $12}'`
echo "Total number of ions: $num_atom"

# create temp files for writing xyz files
typenum=`grep -m 1 'ions per type' $out |head -1 |awk '{print NF}' `
typenum=$(($typenum-4)) # how many types of ions
for i in `seq $typenum`
do
        elename=`grep -m $i POTCAR $out |tail -1 |awk '{print $3}'`
        j=$(($i+4))
        elenum=`grep -m 1 "ions per type" $out |awk -v j=$j '{print $j}'`
        echo $elename $elenum
        for j in `seq $elenum`
        do
            echo $elename  >> type.temp
        done
done
grep -A 3 -m 1 "direct lattice vectors" $out \
|tail -3 |awk '{printf("%f %f %f ",$1,$2,$3)}' >primvec.temp
grep "energy  without " $out |awk '{print $4}' >energy.temp
awk '/POSITION/,/drift/{ if(NF==6) print $0 }' $out  > pos.temp

lines=`wc pos.temp|awk '{print $1}'`
num_str=`echo "$lines/$num_atom" |bc` # how many structures
echo "Number of structures: $num_str"
if [  -f all.xyz ]; then
    rm all.xyz
fi

for i in `seq $num_str`
do
        energy=`head -n $i energy.temp|tail -1`
        echo "$num_atom" >> str_$i.xyz
        echo -n "Lattice=\"" >> str_$i.xyz
        cat primvec.temp >> str_$i.xyz
        echo -n "\" Properties=species:S:1:pos:R:3:forces:R:3 " >> str_$i.xyz
        echo  "energy=$energy pbc=\"T T T\"" >> str_$i.xyz
        end=`echo "$i*$num_atom" |bc `
        head -n $end pos.temp|tail -n $num_atom >pos_i.temp
        paste type.temp pos_i.temp >> str_$i.xyz
        mv str_$i.xyz $name$i.xyz
        cat $name$i.xyz >>all.xyz
done

#rm *.temp
if [ ! -d "struc" ]; then
        mkdir struc
fi
mv *xyz struc
mv struc/all.xyz .
rm *.temp
