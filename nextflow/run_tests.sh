#!/bin/bash

# nf-test do not symlink executable python scripts and relative imports are lost
# when NF_TEST is found in the environment variables, nextflow will add the corresponding modules to the PYTHONPATH thanks to the beforeScript directive of each process
export NF_TEST=True

# Run test on processes
nf-test test tests/modules

# Run test on workflows
