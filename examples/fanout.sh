#!/bin/sh

../benchmark -t random -p r -o fanout -a 2 -d usertime -d maxflow --time-unit s -O 0:0000 -O 1:0001 -O 4:0010 -O 5:0011 -O 16:1000 -O 17:1001 -O 28:1110 -O 29:1111 fanout.rnd
