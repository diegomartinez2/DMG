#!/bin/bash

outcar=""
xmol="no"
moviemol="no"
rms="no"
evs="no"
toten="no"
ekin="no"
etotal="no"
maxforce="no"
frequency="no"
help="yes"
zipcom="cat"

poscar=""

options=$*

display_help()
{
     printf "
      Usage:

       $0 <options> <input file name>

      Synopsis:
         Extracts the geometries and the energies from the
         OUTCAR output file of VASP and displays them on screen.
         This screen output can be redirected to any file.

         If no input filename is give, the default name is \"OUTCAR\"


       Available options:
         (all options are case INSENSITIVE)

        -m   :  Extract the geom. in MOVIEMOL format for visualization.
        -x   :  Extract the geom. in XMOL format for visualization.
        -r   :  Extract the geom. in KS (our KS code) format for RMS.
        -fr  :  Extract the geom. in XMOL format for the frequency mode
                visualization. It also has frequencies.
                The movie for each frequency is written in the
                file \"movie.FREQUENCY_NUMBER.xyz\".
        -e   :  Extract the Eigenvalue spectrum. of all iterations.
        -b   :  Extract the Energies (TOTEN: binding energies) of all
                iterations.
        -k   :  Extract the Kinetic  Energy (EKIN) and temperatures of
                all iterations.
        -t   :  Extract the Total Energies (ETOTAL) of all iterations.
        -f   :  Show the Maximum Force after the last iteration.
        -bz2 :  Extract the data from a BZIPed file.
        -gz  :  Extract the data from a GZIPed file.

        -p <POSCAR>  :  Provide POSCAR file. Useful to get max forces with selective dynamics.

        -o  <input filename>:  Change the input file to this. This is
                               option normally not needed since even
                               without this option (i.e. -o) the input
                               new filename will be used.


   Copyleft  (C) 2004, Sajeev Chacko
   Department of Physics, University of Pune, Pune - 411 007, INDIA.

   Everyone is permitted to copy, change and  distribute verbatim
   copies of this license document.


"
}

set_value()
{
   moviemol=no
   xmol=no
   rms=no
   evs=no
   toten=no
   ekin=no
   etotal=no
   maxforce=no
   frequency=no
   help=no
   case $1 in
      moviemol )  moviemol=yes; return;;
      xmol )      xmol=yes; return;;
      rms )       rms=yes; return;;
      evs )       evs=yes; return;;
      toten )     toten=yes; return;;
      ekin )      ekin=yes; return;;
      etotal )    etotal=yes; return;;
      maxforce )  maxforce=yes; return;;
      frequency ) frequency=yes; return;;
      help )      help=yes; return;;
   esac
}

if [ $# -gt 0 ] ; then
   while [ $# -gt 0 ]
   do
      opt=$1
      case ${opt} in

         -m | -M ) set_value moviemol;;

         -x | -X ) set_value xmol;;

         -r | -R ) set_value rms;;

         -e | -E ) set_value evs;;

         -b | -B ) set_value toten;;

         -k | -K ) set_value ekin;;

         -t | -T ) set_value etotal;;

         -f | -F ) set_value maxforce;;

         -p | -P ) shift; poscar=$1;;

         -fr | -FR | -Fr | -fR ) set_value frequency;;

         -bz2 | -BZ2 | -Bz2 | -bZ2 ) zipcom="bunzip2 -c";;

         -gz | -GZ| -Gz| -gZ ) zipcom="gunzip -c";;

         -o | -O ) shift; outcar="$1";;
         -h ) display_help; exit;;
          * ) outcar="$1"
      esac
      shift
   done
fi


if [ "${help}" = "yes" ] ; then
   display_help
   exit
fi

if [ "${outcar}" = "" ] ; then
   if [ "${zipcom}" = "bunzip2 -c" ] ; then
      outcar=OUTCAR.bz2
   elif [ "${zipcom}" = "gunzip -c" ] ; then
      outcar=OUTCAR.gz
   else
      outcar=OUTCAR
   fi
fi


if [ ! -f "${outcar}" ] ; then
   echo ""
   echo "Error... File \"${outcar}\" does not exists..."
   echo ""
   exit
fi


natom="`${zipcom} ${outcar}  | head -5000 | grep \"  NIONS\"  | head -1 | awk '{print $NF}'`"
nbands="`${zipcom} ${outcar} | head -5000 | grep \"  NBANDS\" | head -1 | awk '{print $NF}'`"



#
#  Check for MOVIEMOL option...
#
if [ "${moviemol}" = "yes" ] ; then
   nlines=`\expr ${natom} + 2`
   nframes=`${zipcom} ${outcar} | grep "POSITION     " | wc -l | awk '{print $1}'`
   # if [ ${nframes} -gt 10000 ] ; then
   #    nframes=10000
   # fi
   echo "${nframes}"
   ${zipcom} ${outcar}                                              \
       |  grep -A ${nlines} "POSITION     "                         \
       |  grep -v "\-\-"                                            \
       |  awk '{printf("%9s %9s %9s   11\n",$1,$2,$3)}'             \
       |  sed -e"s/^ POSITION TOTAL-FORCE (eV\/Angst)   11/${natom}/"
fi


#
#  Check for XMOL option...
#
if [ "${xmol}" = "yes" ] ; then

   nlines=`\expr ${natom} + 2`
   ions_per_type=`${zipcom} ${outcar}|head -2000 | grep "ions per type"`
   ntype=`echo "${ions_per_type}" | wc -w | awk '{print $1-4}'`
   symbols=`${zipcom} ${outcar}|head -2000|grep POTCAR|head -${ntype}|awk '{printf("%s ",$3)}'`
   for ((i=1; i<=ntype; i++))
   do
      atomsymbol[i]=`echo "${symbols}" | awk -v n=${i} '{print $n}'`
   done

   for ((i=1; i<=ntype; i++))
   do
      ((j=4+i))
      atomtype[i]=`echo ${ions_per_type} | awk -v n=${j} '{print $n}'`;
   done

   ((itype=1))
   ((atm_count=0))
   ((iatm_begin=0))
   ((iatm_end=atomtype[itype]))
   for ((iatm=1; iatm<=natom; iatm++))
   do
      ((atm_count=atm_count+1))
      if [ ${iatm} -gt ${iatm_begin} -a ${iatm} -le ${iatm_end} ] ; then
        atom[iatm]=${atomsymbol[itype]}
      fi


      if [ ${atm_count} -eq ${atomtype[itype]} ] ; then
        ((atm_count=0))
        ((itype=itype+1))
        ((iatm_begin=iatm_end))
        ((iatm_end+=atomtype[itype]))
      fi
   done

   symbol_file="symbol_file"
   printf "" > symbol_file
   for ((iatm=1; iatm<=natom; iatm++))
   do
      echo "${atom[iatm]}"  >> ${symbol_file}
   done


   ${zipcom} ${outcar}                                              \
       |  grep -A ${nlines} "POSITION     "                         \
       |  grep -v "\-\-\-\-\-"                                      \
       |  awk '{printf("Na %9s %9s %9s\n",$1,$2,$3)}'               \
       |  sed -e"s/^Na        --/${natom}/"                        \
       |  sed -e"s/^Na  POSITION TOTAL-FORCE (eV\/Angst)//"         \
       |  awk -v symbol_file=${symbol_file} -v natom=${natom} '
                BEGIN {
                   for (iatm=1; iatm<=natom; iatm++) {
                     getline atom[$iatm] < symbol_file;
                     symbol[iatm]=atom[$iatm]
                   }
                   printf("%d\n\n", natom);
                }
                {
                  if ($1=="Na")
                    {
                      jatm++
                      printf("%2s  %9s  %9s  %9s\n", symbol[jatm],$2,$3,$4)
                    }
                  else if ($1==natom)
                    {
                      jatm=0
                      printf("%d\n\n", natom);
                    }
                }
          '

   /bin/rm  -f  ${symbol_file}

fi


#
#  Check for RMS option...
#
if [ "${rms}" = "yes" ] ; then
   nlines=`\expr ${natom} + 2`
   ${zipcom} ${outcar}                                              \
       |  grep -A ${nlines} "POSITION     "                         \
       |  grep -v "\-\-"                                            \
       |  awk '{printf("1 %9s %9s %9s 1\n",$1,$2,$3)}'      \
       |  sed -e"s/^1  POSITION TOTAL-FORCE (eV\/Angst) 1/1/"
fi


#
#  Check for EIGENVALUE option...
#
if [ "${evs}" = "yes" ] ; then
   nlines=`\expr ${nbands} + 0`
   ${zipcom} ${outcar}                                              \
       |  grep -A ${nlines} "^  band No\.  band energies     occupation"             \
       |  grep -v "\-\-"                                            \
       |  sed -e"s/    .0000/   0.0000/g"                           \
       |  sed -e"s/     / /g"                                       \
       |  sed -e"s/    2/ 2/g"                                      \
       |  sed -e"s/    1/ 1/g"                                      \
       |  sed -e"s/    0/ 0/g"                                      \
       |  sed -e"s/ -/   -/g"                                       \
       |  grep -v "band"                                            \
       |  sed -e"s/\.00000//g"
fi


#
#  Check for TOTEN (Total Binding Energy) option...
#
if [ "${toten}" = "yes" ] ; then
   # nlines=`\expr ${nbands} + 2`
   ${zipcom} ${outcar} | grep "  free  energy   TOTEN  =" | grep eV | awk '{print $5}'
fi


#
#  Check for EKIN (Kinetic Energy and Temperature) option...
#
if [ "${ekin}" = "yes" ] ; then
   # nlines=`\expr ${nbands} + 2`
   ${zipcom} ${outcar}|grep "  kinetic Energy EKIN   =  " |awk '{printf(" %s  %6.2f\n",$5,$7)}'
fi


#
#  Check for ETOTAL (Total Energy = P.E. + K.E) option...
#
if [ "${etotal}" = "yes" ] ; then
   # nlines=`\expr ${nbands} + 2`
   ${zipcom} ${outcar} | grep "  total energy   ETOTAL =" | awk '{print $5}'
fi


#
#  Check for MAX FORCE option...
#
if [ "${maxforce}" = "yes" ] ; then
   # printf "\n   Maximum FORCE is  :  "
   nlines=`\expr ${natom} + 2`

   poscar1="${poscar}"
   [ ! -f "${poscar1}" ] && poscar1=""

   ${zipcom} ${outcar}                                \
       |  grep "TOTAL-FORCE" -A ${nlines}             \
       |  grep -v "\-\-\-"                            \
       |  tail -${natom}                             \
       |  awk -v na=${natom} -v poscar="${poscar1}" ' BEGIN {
                   printf ("POSCAR:   %s\n", poscar)
                   if ( poscar != "" ) {
                       start=0
                       lineno=1
                       natoms=0
                       while ( (getline line < poscar) != 0) {
                           if (start==1) {
                               split (line, res, " ")
                               fxres[natoms] = tolower(res[4])
                               fyres[natoms] = tolower(res[5])
                               fzres[natoms] = tolower(res[6])
                               # printf ("AAA:   %3d    %s   %s   %s\n", natoms, res[4], res[5], res[6])
                               natoms++
                           }

                           line_lower=tolower(line)

                           if (lineno>=2 && (line_lower=="direct" || line_lower=="cart" || line_lower=="cartesian") ) {
                               start=1
                           }
                           else if (line==" ")  {
                               start++
                           }

                           lineno++
                       }
                   }

                   fxmax = -999999999999999999.0
                   fymax = -999999999999999999.0
                   fzmax = -999999999999999999.0
                   i = 0
               }
               {
                   fx = $4; if (fx<0) fx = -fx
                   fy = $5; if (fy<0) fy = -fy
                   fz = $6; if (fz<0) fz = -fz

                   if (fxres[i]!="f" && fxmax < fx) {fxmax=fx; coord="X"}
                   if (fyres[i]!="f" && fymax < fy) {fymax=fy; coord="Y"}
                   if (fzres[i]!="f" && fzmax < fz) {fzmax=fz; coord="Z"}

                   frms += fx*fx + fy*fy + fz*fz

                   i++

                }

                END {

                   fmax   =  fxmax
                   coord  =  "X"
                   if (fmax<fymax) {
                       fmax   =  fymax
                       coord  =  "Y"
                   }
                   if (fmax<fzmax) {
                       fmax   =  fzmax
                       coord  =  "Z"
                   }

                   frms = sqrt(frms/(3*na))

                   n = split (i, a, "")
                   if (a[n]==1) ext="st"
                   else if (a[n]==2) ext="nd"
                   else if (a[n]==3) ext="rd"
                   else if (a[n]>=4) ext="th"

                   printf ("\n")
                   printf ("    Maximum FORCE on %2s-coord of %4d%2s atom is: %10.6f", coord, i, ext, fmax)
                   printf ("   -    %10.6f   %10.6f   %10.6f\n", fxmax, fymax, fzmax)
                   printf ("    RMS     FORCE is                           : %10.6f\n", frms)
                 }'
       #|  sed -e"s/-/ /g"                             \
       #|  sort -n                                     \
       #|  tail -1
   echo ""
fi


#
#  Check for Frequency option...
#
if [ "${frequency}" = "yes" ] ; then
   ((nfreq=natom*3-6))
   ((natomp2=natom+2))
   ((nlines=natom+1))
   printf "\n   Extracting movie for frequency visualization"
   printf "\n   No of frequencies  =  ${nfreq}\n\n"

   tmpfile=chacko.vibration.`ps | grep bash | awk '{print $1}'`.tmp
   #while [ -e ${tmpfile} ]
   #do
   #   tmpfile="${tmpfile}.tmp"
   #done

   for ((freq=1; freq<=nfreq; freq++))
   do
      printf "   Extracting movie for frequency no.:   %3d ... " ${freq}
      ${zipcom} ${outcar} | grep -A ${nlines} " ${freq} f  = "    \
                    | tail -${natom}   > ${tmpfile}

      frequency=`${zipcom} ${outcar} | grep " ${freq} f  = "       \
                                     | tail -1               \
                                     | awk '{printf("%s THz    %s cm-1    %s meV\n",$4,$8,$10)}'`
      printf "" > movie.${freq}.xyz
      for ((i=-100; i<=100; i++))
      do
         scale=`echo "${i}*0.08" | bc -lq`
         awk -v n=${natom}  -v s=${scale} -v freq="${frequency}" '
              BEGIN {
                 printf("%d\n\n",n)
              }
              {
                 printf("Ga   %8.5f   %8.5f   %8.5f\n", $1+($4)*s, $2+($5)*s, $3+($6)*s)
              }' ${tmpfile}  >> movie.${freq}.xyz
      done
      for ((i=100; i>=-100; --i))
      do
         scale=`echo "${i}*0.08" | bc -lq`
         awk -v n=${natom}  -v s=${scale} -v freq="${frequency}" '
              BEGIN {
                 printf("%d\n\n",n)
              }
              {
                 printf("Ga   %8.5f   %8.5f   %8.5f\n", $1+($4)*s, $2+($5)*s, $3+($6)*s)
              }' ${tmpfile}  >> movie.${freq}.xyz
      done

      printf "   Done\n"
   done
   /bin/rm  -f  ${tmpfile}
fi
