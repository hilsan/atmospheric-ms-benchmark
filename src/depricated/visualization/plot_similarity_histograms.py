import argparse
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Function to collect the data from the files for each method
def collect_data_for_methods(directory, methods_files):
    data_dict = {}
    
    for method, file_name in methods_files.items():
        method_data = []
        for subdir in sorted(os.listdir(directory)):
            subdir_path = os.path.join(directory, subdir)
            file_path = os.path.join(subdir_path, file_name)
            
            if os.path.isdir(subdir_path) and os.path.exists(file_path):
                try:
                    df = pd.read_csv(file_path)
                    method_data.append(df)
                except Exception as e:
                    print(f"Could not read file {file_path}: {e}")

        if method_data:
            # Combine the data for each method across subdirectories
            data_dict[method] = pd.concat(method_data, ignore_index=True)
        else:
            print(f"No data found for {method}")

    return data_dict

# Plot the similarity histograms for each method
def plot_similarity_histograms(data_dict, output_dir):
    """Plot similarity score histograms for multiple methods."""
    metrics = ["Cosine", "Cosine_Weighted", "Tanimoto", "Dot_Product"]
    
    bins = np.linspace(0, 1000, 21)  # Define the bin edges for histograms
    colors = {"NEIMS": "#E059B0", "QCxMS": "#8CD3FA", "QCxMS2": "#A9CD6F"}

    # Loop through each metric and create a separate figure for each one
    for metric in metrics:
        fig, ax = plt.subplots(figsize=(10, 7.5), dpi=100)  # 4:3 aspect ratio (10, 7.5)

        # Plot data for each method
        for method, df in data_dict.items():
            if not df.empty:
                # Normalize the data by dividing by the total count of each method
                counts, _ = np.histogram(df[metric], bins=bins)
                normalized_counts = counts / counts.sum()  # Normalize to percentage of total
                
                # Plot the histogram bars for the method
                ax.bar(bins[:-1], normalized_counts, width=np.diff(bins), label=method, 
                       color=colors.get(method, "gray"), alpha=0.6)

        # Set title and labels
        ax.set_title(f"{metric} Similarity", fontsize=24)
        ax.set_xlabel("Similarity Score", fontsize=24)
        ax.set_ylabel("Frequency (Percentage)", fontsize=24)
        ax.set_xlim(0, 1000)
        ax.grid(True, alpha=0.3)

        # Set ticks for both axes
        ax.set_xticks(bins[::4])  # Set x-tick marks
        ax.set_xticklabels(np.round(bins[::4], 1), fontsize=24)
        ax.tick_params(axis='y', labelsize=24)  # Set y-tick labels size

        # Create legend with 24-point font, place it in the upper left corner, and make it transparent
        ax.legend(title="", loc="upper left", fontsize=24, frameon=False, bbox_to_anchor=(0, 1))

        # Adjust layout to prevent clipping
        plt.tight_layout()

        # Save the figure as a separate PNG file
        output_file = os.path.join(output_dir, f"{metric}_similarity_histogram.png")
        plt.savefig(output_file, dpi=300, transparent=True)
        plt.close(fig)  # Close the figure to free up memory

def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Plot similarity histograms.")
    
    # Argument for the directory where the subdirectories are located
    parser.add_argument("--dir", required=True, help="Parent directory containing the subdirectories.")
    
    # Argument for output directory
    parser.add_argument("--output", required=True, help="Output directory for the plot.")
    
    # Flags for enabling methods (NEIMS, QCxMS, QCxMS2)
    parser.add_argument("--qcxms_flag", action="store_true", help="Include QCxMS data.")
    parser.add_argument("--qcxms2_flag", action="store_true", help="Include QCxMS2 data.")
    parser.add_argument("--qcxms2_dft_flag", action="store_true", help="Include QCxMS2 data.")
    
    
    # Argument for specifying filenames
    parser.add_argument("--neims", default="exp2neims_similarity_scores.csv", help="Filename for NEIMS similarity scores.")
    parser.add_argument("--qcxms", default="exp2qcxms_similarity_scores.csv", help="Filename for QCxMS similarity scores.")
    parser.add_argument("--qcxms2", default="exp2qcxms2_similarity_scores.csv", help="Filename for QCxMS2 similarity scores.")
    parser.add_argument("--qcxms2_dft", default="exp2qcxms-dft_similarity_scores.csv", help="Filename for QCxMS2 similarity scores.")

    # Parse the arguments
    args = parser.parse_args()

    # Define the directory and methods to load
    directory = args.dir
    methods_files = {
        "NEIMS": args.neims
    }

    # Only include QCxMS and QCxMS2 if the flags are set
    if args.qcxms_flag:
        methods_files["QCxMS"] = args.qcxms
    if args.qcxms2_flag:
        methods_files["QCxMS2"] = args.qcxms2
    if args.qcxms2_dft_flag:
        methods_files["QCxMS2-dft"] = args.qcxms2_dft

    # Collect data
    data = collect_data_for_methods(directory, methods_files)

    # Plot the histograms and save them as separate files
    plot_similarity_histograms(data, args.output)

if __name__ == "__main__":
    main()
