#!/usr/bin/env bash
#
#  ciclo_VASP_SSCHA.sh
#
#  Copyright 2023 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
# First copy the files to the cluster:
# scp [options] username1@source_host:directory1/filename1 username2@destination_host:directory2/filename2
scp -r vasp diegom@diegom@ekhi.cfm.ehu.es:~/vasp
# Second, run the VASP
# ssh -l ${USERNAME} ${HOSTNAME} "${SCRIPT}"
ssh -l diegom ekhi.cfm.ehu.es "cd ~/vasp; sbatch run.bash"
# Third, check if is still running (¿¿¿cómo hacerlo???)
if `squeue |grep diegom|awk -F, '{print $NF}'`==$codigo; then continua=True ; else continua=False; fi

# Forth, if finished retrieve the files back

scp -r diegom@diegom@ekhi.cfm.ehu.es:~/vasp vasp
