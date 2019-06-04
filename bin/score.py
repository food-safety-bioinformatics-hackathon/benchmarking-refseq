#!/usr/bin/env python
# This will load in a BAMCov output file and look at the coverage and mean depth
import sys
import os

results_file='/usr/users/QIB_fr017/lott/Downloads/bamcov/test.tsv' # TODO: Change this path to the results file name or specify the file name as a command line argument

print(sys.argv)
if len(sys.argv) > 1: results_file = sys.argv[1]

def getCoverageAndMeanDepth(filename):
    with open(filename) as f:
        lines = f.readlines()
        return getCoverageAndMeanDepthFromLine(lines[0])

def getCoverageAndMeanDepthFromLine(fileline):
        # Assume just one contpig so return the coverage and mean depth for that
        field = fileline.split('\t')

        return field[5],field[6]

sum_of_coverage_differences=0
coverage_differences_squared=0

sum_of_mean_depth_differences=0
mean_depth_differences_squared=0

with open(results_file) as f:
    lines = f.readlines()
    referenceLine=lines[len(lines)-1]
    print('Reference',referenceLine)

    true_ref_average_coverage, true_ref_mean_depth = getCoverageAndMeanDepthFromLine(referenceLine)

    for line in lines:
        print('Processing ',line)

        average_coverage, mean_depth = getCoverageAndMeanDepthFromLine(line)

        coverage_difference = float(true_ref_average_coverage)-float(average_coverage)
        mean_depth_difference = float(true_ref_mean_depth)-float(mean_depth)

        sum_of_coverage_differences+=abs(coverage_difference)
        coverage_differences_squared+=coverage_difference*coverage_difference

        sum_of_mean_depth_differences+=abs(mean_depth_difference)
        mean_depth_differences_squared+=mean_depth_difference*mean_depth_difference

        print('Difference in Coverage',coverage_difference)
        print('Difference in Mean Depth',mean_depth_difference)

print('\nSum of coverage differences',sum_of_coverage_differences)
print('Coverage differences squared',coverage_differences_squared)

print('Sum of mean differences',sum_of_mean_depth_differences)
print('Mean Depth differences squared',mean_depth_differences_squared)