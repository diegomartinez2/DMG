#!/bin/bash
echo "🔬 ULTRA-FINA chimesFF (10 min)"
lmp_mpi_chimes -in min_ultra_fina.in > min_ultra.log
echo "🏆 RESULTADO:"
grep "FINAL ULTRA-FINA" min_ultra.log
echo "✅ Listo para Annealing!"
