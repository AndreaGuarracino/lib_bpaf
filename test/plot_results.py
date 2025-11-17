#!/usr/bin/env python3
"""
Plot lib_bpaf compression test results from all_results.tsv

Usage:
    python3 plot_results.py [path_to_all_results.tsv]

Defaults to: ../test/bpaf_all_tests/all_results.tsv (relative to this script)
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path
import numpy as np

def format_bytes(bytes_val):
    """Convert bytes to human-readable format"""
    for unit in ['B', 'KB', 'MB', 'GB']:
        if bytes_val < 1024.0:
            return f"{bytes_val:.2f} {unit}"
        bytes_val /= 1024.0
    return f"{bytes_val:.2f} TB"

def plot_dataset_metrics(df, dataset_name, output_dir):
    """Create comprehensive plots for a single dataset"""

    # Filter data for this dataset
    data = df[df['dataset_name'] == dataset_name].copy()

    if len(data) == 0:
        print(f"Warning: No data found for dataset '{dataset_name}'")
        return

    # Sort by compression ratio (best first)
    data = data.sort_values('ratio_orig_to_bpaf', ascending=False)

    # Create strategy labels (combine first → second strategies)
    strategies = []
    for _, row in data.iterrows():
        first = row['strategy_first']
        second = row['strategy_second']
        strategies.append(f"{first}→{second}")

    # Create figure with 4 subplots (vertical stack)
    fig, axes = plt.subplots(4, 1, figsize=(72, 32))
    fig.suptitle(f'Compression Strategy Performance: {dataset_name}', fontsize=16, fontweight='bold', y=0.995)

    # Original file info
    orig_size = data['original_size_bytes'].iloc[0]
    num_records = data['num_records'].iloc[0]
    dataset_type = data['dataset_type'].iloc[0]

    fig.text(0.5, 0.975, f'Type: {dataset_type} | Records: {num_records:,} | Original Size: {format_bytes(orig_size)}',
             ha='center', fontsize=11, style='italic')

    # Color code by compression layer (based on requested strategy)
    colors = []
    for _, row in data.iterrows():
        strat = row['compression_strategy']
        if strat.endswith('-bgzip'):
            colors.append('#ff7f0e')  # orange for bgzip
        elif strat.endswith('-nocomp'):
            colors.append('#2ca02c')  # green for nocomp
        else:
            colors.append('#1f77b4')  # blue for zstd (default)

    # Plot 1: BPAF File Size
    ax1 = axes[0]
    bars1 = ax1.bar(range(len(strategies)), data['bpaf_size_bytes'], color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax1.set_xlabel('Compression Strategy', fontweight='bold')
    ax1.set_ylabel('BPAF File Size (bytes)', fontweight='bold')
    ax1.set_title('Final Compressed File Size', fontweight='bold')
    ax1.set_xticks(range(len(strategies)))
    ax1.set_xticklabels(strategies, rotation=90, ha='right', fontsize=8)
    ax1.grid(axis='y', alpha=0.3)
    ax1.ticklabel_format(axis='y', style='plain')

    # Annotate best (smallest)
    best_idx = data['bpaf_size_bytes'].idxmin()
    best_size = data.loc[best_idx, 'bpaf_size_bytes']
    best_first = data.loc[best_idx, 'strategy_first']
    best_second = data.loc[best_idx, 'strategy_second']
    ax1.text(0.98, 0.98, f'Best: {best_first}→{best_second}\n{format_bytes(best_size)}',
             transform=ax1.transAxes, ha='right', va='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), fontsize=9)

    # Plot 2: Compression Ratio (Original → BPAF)
    ax2 = axes[1]
    bars2 = ax2.bar(range(len(strategies)), data['ratio_orig_to_bpaf'], color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax2.set_xlabel('Compression Strategy', fontweight='bold')
    ax2.set_ylabel('Compression Ratio (x)', fontweight='bold')
    ax2.set_title('End-to-End Compression Ratio (Original → BPAF)', fontweight='bold')
    ax2.set_xticks(range(len(strategies)))
    ax2.set_xticklabels(strategies, rotation=90, ha='right', fontsize=8)
    ax2.grid(axis='y', alpha=0.3)

    # Annotate best (highest)
    best_idx = data['ratio_orig_to_bpaf'].idxmax()
    best_ratio = data.loc[best_idx, 'ratio_orig_to_bpaf']
    best_first = data.loc[best_idx, 'strategy_first']
    best_second = data.loc[best_idx, 'strategy_second']
    ax2.text(0.98, 0.98, f'Best: {best_first}→{best_second}\n{best_ratio:.2f}x',
             transform=ax2.transAxes, ha='right', va='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8), fontsize=9)

    # Plot 3: Seek Time (Mode B)
    ax3 = axes[2]
    bars3 = ax3.bar(range(len(strategies)), data['seek_mode_b_avg_us'], color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax3.set_xlabel('Compression Strategy', fontweight='bold')
    ax3.set_ylabel('Average Seek Time (μs)', fontweight='bold')
    ax3.set_title('Seek Performance (Mode B - Standalone Functions)', fontweight='bold')
    ax3.set_xticks(range(len(strategies)))
    ax3.set_xticklabels(strategies, rotation=90, ha='right', fontsize=8)
    ax3.grid(axis='y', alpha=0.3)

    # Annotate best (lowest)
    best_idx = data['seek_mode_b_avg_us'].idxmin()
    best_seek = data.loc[best_idx, 'seek_mode_b_avg_us']
    best_first = data.loc[best_idx, 'strategy_first']
    best_second = data.loc[best_idx, 'strategy_second']
    ax3.text(0.98, 0.98, f'Fastest: {best_first}→{best_second}\n{best_seek:.2f} μs',
             transform=ax3.transAxes, ha='right', va='top',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8), fontsize=9)

    # Plot 4: Compression Time
    ax4 = axes[3]
    bars4 = ax4.bar(range(len(strategies)), data['compression_runtime_sec'], color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
    ax4.set_xlabel('Compression Strategy', fontweight='bold')
    ax4.set_ylabel('Compression Time (seconds)', fontweight='bold')
    ax4.set_title('Compression Runtime', fontweight='bold')
    ax4.set_xticks(range(len(strategies)))
    ax4.set_xticklabels(strategies, rotation=90, ha='right', fontsize=8)
    ax4.grid(axis='y', alpha=0.3)

    # Annotate fastest
    best_idx = data['compression_runtime_sec'].idxmin()
    best_time = data.loc[best_idx, 'compression_runtime_sec']
    best_first = data.loc[best_idx, 'strategy_first']
    best_second = data.loc[best_idx, 'strategy_second']
    ax4.text(0.98, 0.98, f'Fastest: {best_first}→{best_second}\n{best_time:.2f} s',
             transform=ax4.transAxes, ha='right', va='top',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8), fontsize=9)

    # Add legend for compression layers
    zstd_patch = mpatches.Patch(color='#1f77b4', alpha=0.7, label='Zstd (default)')
    bgzip_patch = mpatches.Patch(color='#ff7f0e', alpha=0.7, label='BGZIP')
    nocomp_patch = mpatches.Patch(color='#2ca02c', alpha=0.7, label='No compression')
    fig.legend(handles=[zstd_patch, bgzip_patch, nocomp_patch],
               loc='lower center', ncol=3, frameon=True, fontsize=10,
               bbox_to_anchor=(0.5, -0.02))

    plt.subplots_adjust(left=0.05, right=0.98)
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])

    # Save plot
    output_file = output_dir / f'{dataset_name}_metrics.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")
    plt.close()


def generate_summary_table(df, output_dir):
    """Generate a summary markdown table showing top performers per dataset"""

    output_file = output_dir / 'summary_best_strategies.md'

    with open(output_file, 'w') as f:
        f.write("# Best Compression Strategies per Dataset\n\n")
        f.write("Summary of top-performing strategies across all metrics.\n\n")

        for dataset in df['dataset_name'].unique():
            data = df[df['dataset_name'] == dataset].copy()

            f.write(f"## {dataset}\n\n")
            f.write(f"- **Type**: {data['dataset_type'].iloc[0]}\n")
            f.write(f"- **Records**: {data['num_records'].iloc[0]:,}\n")
            f.write(f"- **Original Size**: {format_bytes(data['original_size_bytes'].iloc[0])}\n\n")

            # Best compression ratio
            best_ratio_idx = data['ratio_orig_to_bpaf'].idxmax()
            f.write(f"### Best Compression Ratio\n")
            first = data.loc[best_ratio_idx, 'strategy_first']
            second = data.loc[best_ratio_idx, 'strategy_second']
            f.write(f"- **Strategy**: {first}→{second}\n")
            f.write(f"- **Ratio**: {data.loc[best_ratio_idx, 'ratio_orig_to_bpaf']:.2f}x\n")
            f.write(f"- **BPAF Size**: {format_bytes(data.loc[best_ratio_idx, 'bpaf_size_bytes'])}\n\n")

            # Fastest seek
            best_seek_idx = data['seek_mode_b_avg_us'].idxmin()
            f.write(f"### Fastest Seek (Mode B)\n")
            first = data.loc[best_seek_idx, 'strategy_first']
            second = data.loc[best_seek_idx, 'strategy_second']
            f.write(f"- **Strategy**: {first}→{second}\n")
            f.write(f"- **Time**: {data.loc[best_seek_idx, 'seek_mode_b_avg_us']:.2f} μs\n")
            f.write(f"- **Ratio**: {data.loc[best_seek_idx, 'ratio_orig_to_bpaf']:.2f}x\n\n")

            # Fastest compression
            best_comp_idx = data['compression_runtime_sec'].idxmin()
            f.write(f"### Fastest Compression\n")
            first = data.loc[best_comp_idx, 'strategy_first']
            second = data.loc[best_comp_idx, 'strategy_second']
            f.write(f"- **Strategy**: {first}→{second}\n")
            f.write(f"- **Time**: {data.loc[best_comp_idx, 'compression_runtime_sec']:.2f} s\n")
            f.write(f"- **Ratio**: {data.loc[best_comp_idx, 'ratio_orig_to_bpaf']:.2f}x\n\n")

            # Best overall (balance of ratio and seek time)
            data['score'] = data['ratio_orig_to_bpaf'] / (data['seek_mode_b_avg_us'] / 10.0)
            best_overall_idx = data['score'].idxmax()
            f.write(f"### Best Overall Balance (Ratio/Seek)\n")
            first = data.loc[best_overall_idx, 'strategy_first']
            second = data.loc[best_overall_idx, 'strategy_second']
            f.write(f"- **Strategy**: {first}→{second}\n")
            f.write(f"- **Ratio**: {data.loc[best_overall_idx, 'ratio_orig_to_bpaf']:.2f}x\n")
            f.write(f"- **Seek**: {data.loc[best_overall_idx, 'seek_mode_b_avg_us']:.2f} μs\n\n")

            f.write("---\n\n")

    print(f"  ✓ Saved: {output_file}")


def main():
    # Determine input file
    if len(sys.argv) > 1:
        tsv_file = Path(sys.argv[1])
    else:
        # Default path relative to script location
        script_dir = Path(__file__).parent
        tsv_file = script_dir / 'bpaf_all_tests' / 'all_results.tsv'

    if not tsv_file.exists():
        print(f"Error: TSV file not found: {tsv_file}")
        print(f"\nUsage: {sys.argv[0]} [path_to_all_results.tsv]")
        sys.exit(1)

    print(f"Reading test results from: {tsv_file}")

    # Read TSV
    df = pd.read_csv(tsv_file, sep='\t')

    print(f"Loaded {len(df)} test results")
    print(f"Found {df['dataset_name'].nunique()} unique datasets")
    print()

    # Output directory (same as TSV)
    output_dir = tsv_file.parent

    # Generate plots for each dataset
    print("Generating plots...")
    for dataset in sorted(df['dataset_name'].unique()):
        print(f"  Processing: {dataset}")
        plot_dataset_metrics(df, dataset, output_dir)

    print()
    print("Generating summary table...")
    generate_summary_table(df, output_dir)

    print()
    print("=" * 60)
    print("All plots generated successfully!")
    print(f"Output directory: {output_dir}")
    print("=" * 60)


if __name__ == '__main__':
    main()
