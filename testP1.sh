#!/bin/bash

module add gridmathematica-9
  
INPUT_FILE="/auto/plzen1/home/yawa/DEMES/testPlzen2.m"

math-grid start

math -run <$INPUT_FILE
