### How to Compile
(1) Downloads 5 files: main.cpp, data.h, solver.h, solver.cpp, Makefile <br>
(2) enter the following command:
  ```
  make
  ```
  It will generate the executable file "legalizer" in the current directory<br>

If you want to remove it, please enter the following command:<br>
  ```
  make clean
  ```

### How to Run
enter the following command:<br>
```
./legalizer <first input file path> <second input file path> <third input file path> <output file path>
```

e.g.
```
./legalizer ./input/testcase1/architecture.txt ./input/testcase1/instance.txt ./input/testcase1/netlist.txt ./output/output1.txt
```

### Testcases
Decompress the public_cases.tar.gz to get the public testcases
```
tar zxvf public_cases.tar.gz
```
