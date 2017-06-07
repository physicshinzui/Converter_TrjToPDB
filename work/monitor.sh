#!/bin/bash

for i in $(seq 1000); do
  grep Conf log1 | wc -l
  sleep 10s
done


