# Variant calling workflow

In this notebook pipeline used for variant calling will be described. Pipeline will run mainly in batch 
command-line.

### Getting the files and the programs
Pipeline will be runing in its own virtual enviroment - to be sure, that none of the dependencies 
for this pipeline can create problems with already installed programs.
```bash
virtualenv exercise_1
source exercise_1/bin/activate
sudo apt install -y bwa
sudo apt install -y bcftools
```
Now we will get files, which we will be working with, locally.

```bash
```
