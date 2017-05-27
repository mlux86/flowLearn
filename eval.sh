#!/bin/bash

mkdir -p results

echo "impc.bonemarrow"
for i in 1 2 5 10 50; do
    echo "  $i"
    R <<< "source('eval.R'); flEvalDataset('impc.bonemarrow', $i, '/media/data/mlux/flowlearn.traindat/flowlearn.traindat.')" --no-save
done

echo "flowcap.bcell"
for i in 1 4 7 11 20; do
    echo "  $i"
    R <<< "source('eval.R'); flEvalDataset('flowcap.bcell', $i, '/media/data/mlux/flowlearn.traindat/flowlearn.traindat.')" --no-save
done

echo "flowcap.tcell"
for i in 1 4 7 11 20; do
    echo "  $i"
    R <<< "source('eval.R'); flEvalDataset('flowcap.tcell', $i, '/media/data/mlux/flowlearn.traindat/flowlearn.traindat.')" --no-save
done

echo "flowcap.DC"
for i in 1 4 7 11 20; do
    echo "  $i"
    R <<< "source('eval.R'); flEvalDataset('flowcap.DC', $i, '/media/data/mlux/flowlearn.traindat/flowlearn.traindat.')" --no-save
done

echo "flowcap.treg"
for i in 1 4 7 11 20; do
    echo "  $i"
    R <<< "source('eval.R'); flEvalDataset('flowcap.treg', $i, '/media/data/mlux/flowlearn.traindat/flowlearn.traindat.')" --no-save
done
