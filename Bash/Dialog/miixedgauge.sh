#! /bin/sh
# $Id: mixedgauge,v 1.7 2010/01/13 10:20:03 tom Exp $

while true
do
dialog	--mixedgauge "DFTB+ Calculations." \
		0 0 10 \
		"Process 01"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-10000"|bc)"\
		"Process 02"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-1000"|bc)"\
		"Process 03"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-400"|bc)"\
		"Process 04"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"\
		"Process 05"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"\
		"Process 06"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"\
		"Process 07"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"\
		"Process 08"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"\
		"Process 09"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"\
		"Process 10"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"\
		"Process 11"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"\
		"Process 12"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"\
		"Process 13"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"\
		"Process 14"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"\
		"Process 15"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"\
		"Process 16"	"$(echo $(tail -300 /home/r2d2/DFTB+/dftb+mpi-negf.r4732.x86_64-linux/Device_G17.8G_S0.045/geo_end.xyz|grep Step|grep -o "[^ ]*$")"/-300"|bc)"

# break
sleep 5
done
