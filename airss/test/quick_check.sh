#!/bin/bash
set -e

buildcell < quick_check.cell 2>&1 | grep -v " options " | grep -v "# on" | grep -v "# at" | grep -v "Total:" | grep -v "Timed:" | grep -v "Perc:" | awk '{$1=$1};1' > test.quick_check
airss.pl -pp3 -max 1 -seed quick_check
awk '{$1=$1};1' quick_check-*.res | grep -v "REM build" | grep -v " options " | grep -v 'REM Run started:' | grep -v 'TITL quick_check' | grep -v 'runtime:' | grep -v 'efficiency:' | grep -v -i version >> test.quick_check

# Test cabal conversions
echo "Testing cabal conversion..."
for suffix in cell shx res gulp cif psi4 xtl xyz poscar qe cp2k conf
do
	cabal res ${suffix} 0.1 < quick_check-*.res > converted.${suffix}
done

for suffix in cell shx xtl xyz poscar qe cp2k conf
do
	cabal ${suffix} res 0.1 < converted.${suffix} > converted.res
done
echo "cabal conversion tests finished"

rm *.minsep quick_check-*.* quick_check.cell.* converted.*
cat test.quick_check
diff test.quick_check benchmark.quick_check && exit 0 || exit 1
