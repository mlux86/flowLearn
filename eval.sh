#!/bin/bash

mkdir -p results

echo "flowcap.bcell"
for i in 1 4 7 11 20; do
    echo $i
    time (R <<< "source('eval.R'); flEvalDataset('flowcap.bcell', $i, '/home/mlux/flowlearn.traindat/flowlearn.traindat.')" --no-save > R.log 2>&1)
done

echo "flowcap.tcell"
for i in 1 4 7 11 20; do
    echo $i
    time (R <<< "source('eval.R'); flEvalDataset('flowcap.tcell', $i, '/home/mlux/flowlearn.traindat/flowlearn.traindat.')" --no-save > R.log 2>&1)
done

echo "flowcap.DC"
for i in 1 4 7 11 20; do
    echo $i
    time (R <<< "source('eval.R'); flEvalDataset('flowcap.DC', $i, '/home/mlux/flowlearn.traindat/flowlearn.traindat.')" --no-save > R.log 2>&1)
done

echo "flowcap.treg"
for i in 1 4 7 11 20; do
    echo $i
    time (R <<< "source('eval.R'); flEvalDataset('flowcap.treg', $i, '/home/mlux/flowlearn.traindat/flowlearn.traindat.')" --no-save > R.log 2>&1)
done

echo "impc.bonemarrow"
for i in 1 2 5 10 50; do
    echo $i
    time (R <<< "source('eval.R'); flEvalDataset('impc.bonemarrow', $i, '/home/mlux/flowlearn.traindat/flowlearn.traindat.')" --no-save > R.log 2>&1)
done

