# Running the Fortran experiments

To perform the experiments contained here you need to follow the following steps
1. Have compiled and installed PSBLAS and AMG4PSBLAS as described in the [README](/configure/README.md) under the configure folder.
2. Compiled the example executable in the fortransrc folder as described in the attached [README](/fortransrc/README.md).
3. Modify the file `katz_test.sh` that is currently configured to use the Toeplitz cluster of the University of Pisa to the
   desired settings and run it with
   ```bash
   sbatch katz_test.sh
   ```
4. Log files with results will be generated in the [katzolution](/runs/katzsolution) folder.
