# Coverage Visualization Script

## Overview

The coverage_visualisation.py script is a comprehensive Python tool designed to analyze and visualize sequencing coverage data from genome assemblies. It takes samtools coverage output and optional mapping statistics to generate detailed plots and summary statistics for assembly quality assessment.

## Features

### 📊 **Visual Analysis**
- **6-panel comprehensive visualization** including:
  - Coverage distribution histogram
  - Coverage vs contig length scatter plot
  - Cumulative coverage distribution
  - Contig length distribution
  - Coverage statistics summary panel
  - Quality & mapping metrics panel

### 📈 **Statistical Analysis**
- **Coverage metrics**: mean, median, weighted average, peak, minimum, standard deviation
- **Assembly quality metrics**: N50, total bases, contig counts, length distributions
- **Mapping statistics**: read counts, mapping rates, properly paired reads
- **Quality scores**: base quality and mapping quality analysis

### 📁 **Output Files**
- High-resolution PNG visualization plots
- Detailed summary statistics text files
- Comprehensive console output with key metrics

## Requirements

### Dependencies
```bash
pip install pandas matplotlib numpy
```

### Input Files
1. **Coverage file** (required): samtools coverage output in TSV format
2. **Flagstat file** (optional): samtools flagstat output for mapping statistics

## Usage

### Basic Usage
```bash
python coverage_visualisation.py <coverage_file> <output_plot.png>
```

### Full Usage with Options
```bash
python coverage_visualisation.py coverage.txt plot.png \
    --sample "SampleName" \
    --assembler "AssemblerName" \
    --flagstat flagstat.txt \
    --summary summary.txt
```

### Command Line Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `coverage_file` | ✅ | Path to samtools coverage output file |
| `output_file` | ✅ | Output path for visualization plot (PNG/PDF) |
| `--sample` | ❌ | Sample name for plot titles (default: "Sample") |
| `--assembler` | ❌ | Assembler name for plot titles (default: "Assembler") |
| `--flagstat` | ❌ | Path to samtools flagstat file for mapping stats |
| `--summary` | ❌ | Custom path for summary text file |

## Input File Formats

### 1. Coverage File (samtools coverage output)
Expected columns:
- `#rname`: Contig/reference name
- `startpos`: Start position (usually 1)
- `endpos`: End position (contig length)
- `numreads`: Number of mapped reads
- `covbases`: Number of covered bases
- `coverage`: Coverage percentage
- `meandepth`: Mean coverage depth
- `meanbaseq`: Mean base quality
- `meanmapq`: Mean mapping quality

### 2. Flagstat File (samtools flagstat output)
Expected format:
```
13479 + 0 in total (QC-passed reads + QC-failed reads)
2495 + 0 mapped (18.51% : N/A)
2340 + 0 properly paired (17.38% : N/A)
...
```

## Output Description

### 📊 **Visualization Plots (6-panel layout)**

#### **Panel 1: Coverage Distribution Histogram**
- Shows frequency distribution of coverage depths across contigs
- Vertical lines indicate mean, median, and peak coverage
- Helps identify coverage uniformity and outliers

#### **Panel 2: Coverage vs Contig Length Scatter Plot**
- X-axis: Contig length (log scale)
- Y-axis: Coverage depth
- Color: Coverage percentage (viridis colormap)
- Shows relationship between contig size and coverage

#### **Panel 3: Cumulative Coverage Distribution**
- Cumulative count of contigs vs coverage depth
- Reference lines for mean and median coverage
- Useful for understanding coverage distribution shape

#### **Panel 4: Contig Length Distribution**
- Histogram of contig lengths (log scale)
- Red line at 250bp threshold (common quality cutoff)
- Shows assembly fragmentation characteristics

#### **Panel 5: Coverage Statistics Panel**
- Text summary of key coverage metrics
- Assembly information and quality indicators
- Centered display with enhanced formatting

#### **Panel 6: Quality & Mapping Metrics Panel**
- Mapping statistics (if flagstat provided)
- Base and mapping quality scores
- Read distribution information

### 📋 **Summary Statistics File**

Contains all calculated metrics in text format:
```
Coverage Analysis Summary - SampleName (AssemblerName)
============================================================

sample: SampleName
assembler: AssemblerName
total_contigs: 1234
total_bases: 2500000
mean_coverage: 15.43
weighted_mean_coverage: 18.67
median_coverage: 12.34
n50: 45678
...
```

### 🖥️ **Console Output**

Real-time progress and summary information:
```
Loaded 1234 contigs from coverage.txt
Loaded mapping statistics from flagstat.txt

=== Coverage Summary for SampleName (AssemblerName) ===
Total contigs: 1,234
Contigs < 250bp: 156 (12.6%)
Total bases: 2,500,000
Average coverage: 15.43x
N50: 45,678 bp

=== Mapping Statistics ===
Total reads: 13,479
Mapped reads: 2,495
Mapping rate: 18.5%

Plot saved to: output_plot.png
Summary statistics saved to: output_summary.txt
```

## Key Metrics Explained

### **Coverage Metrics**
- **Mean Coverage**: Simple average of all contig coverages
- **Weighted Mean Coverage**: Average weighted by contig length (more representative)
- **Overall Coverage %**: Percentage of bases with any coverage
- **Peak Coverage**: Maximum coverage observed

### **Assembly Quality Metrics**
- **N50**: Length of shortest contig at 50% of total assembly length
- **Contigs < 250bp**: Count and percentage of small contigs (quality indicator)
- **Total Bases**: Sum of all contig lengths

### **Mapping Quality Metrics**
- **Mapping Rate**: Percentage of reads that mapped to assembly
- **Properly Paired**: Reads mapping in correct orientation and distance
- **Mean Base Quality**: Average quality score of mapped bases
- **Mean Mapping Quality**: Average confidence in read mapping positions

## Integration with Snakemake

This script is designed to work within Snakemake workflows:

```python
rule coverage_visualization:
    input:
        coverage="{sample}/coverage_stats.txt",
        flagstat="{sample}/flagstat.txt"
    output:
        plot="{sample}/coverage_plot.png",
        summary="{sample}/coverage_summary.txt"
    shell:
        """
        python workflow/scripts/coverage_visualisation.py {input.coverage} {output.plot} \
            --sample {wildcards.sample} \
            --assembler spades \
            --flagstat {input.flagstat} \
            --summary {output.summary}
        """
```

## Troubleshooting

### Common Issues

1. **Missing Dependencies**
   ```bash
   pip install pandas matplotlib numpy
   ```

2. **File Format Errors**
   - Ensure coverage file is tab-separated
   - Check that required columns are present
   - Verify file permissions and paths

3. **Memory Issues with Large Assemblies**
   - Script automatically closes plots to free memory
   - Consider filtering very small contigs if needed

4. **Plot Quality**
   - Output is saved at 300 DPI for publication quality
   - Use PNG for raster or PDF for vector graphics

### Performance Tips

- **Large assemblies**: Script handles millions of contigs efficiently
- **Memory usage**: Plots are closed automatically to conserve memory
- **File size**: Summary files are compact text format

## Example Use Cases

1. **Assembly Quality Assessment**: Compare coverage uniformity across different assemblers
2. **Sequencing Depth Analysis**: Evaluate if sequencing depth is sufficient
3. **Contamination Detection**: Identify contigs with unusual coverage patterns
4. **Assembly Optimization**: Guide parameter tuning for assembly tools