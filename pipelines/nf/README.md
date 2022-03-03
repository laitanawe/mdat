# Reproducible Research with Nextflow Pipelines
Nextflow is a workflow engine for data analysis pipelines, with strong focus on reproducibility, portability and scalability of research.

## Installation

In case you're trying to install Nextflow on a machine, please note that Nextflow requires java8 or later.
Here are the steps for Nextflow installation in your home directory.

Step1:
- Make a directory for nextflow in your home directory
```
mkdir -pv ~/opt/nextflowrm -rfv ~/opt/nextflow/*
```

Step2:
- cd into the directory that you just made for nextflow
```
cd ~/opt/nextflow
```

Step3:
- Download the nextflow installation and install it
```
curl -s https://get.nextflow.io| bash
```

Step4:
- Link the nextflow binary to your bin directory
```
ln -s "$PWD/nextflow" ~/bin/nextflow~/opt/nextflow/nextflow self-update
```

Step5:
- Check the version of Nextflow to confirm that it's installed
```
nextflow-v
```
