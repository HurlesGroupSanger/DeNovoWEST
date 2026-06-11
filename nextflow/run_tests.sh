#!/bin/bash

# nf-test do not symlink executable python scripts and relative imports are lost
# when NF_TEST is found in the environment variables, nextflow will add the corresponding modules to the PYTHONPATH thanks to the beforeScript directive of each process
export NF_TEST=True
# ~/.nextflow may not be writable on some systems; use a local cache dir instead
export NXF_HOME="${NXF_HOME:-/tmp/nxf_home_$(whoami)}"

# Run test on processes
nf-test test nf_tests/modules

# Run test on workflows
