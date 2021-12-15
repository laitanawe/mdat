#!/usr/bin/env bash

mysleeptime=3
while [ 0 ]; do squeue -u $(whoami); echo; sleep $mysleeptime; done

## From the command line
## ./monitor_queue.sh
