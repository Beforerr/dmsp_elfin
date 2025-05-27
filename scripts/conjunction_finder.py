#!/usr/bin/env python3
"""
DMSP-ELFIN Conjunction Finder

This script identifies conjunction events between DMSP and ELFIN satellites
based on their MLT (Magnetic Local Time) and MLAT (Magnetic Latitude) values.

A conjunction is defined when:
- |MLT_DMSP - MLT_ELFIN| < d1
- |MLAT_DMSP - MLAT_ELFIN| < d2
- This condition persists for at least T_min

Author: Cascade AI
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import os
import argparse


class ConjunctionFinder:
    """Class to find conjunction events between DMSP and ELFIN satellites."""
    
    def __init__(self, mlt_threshold=0.5, mlat_threshold=2.0, min_duration_seconds=60):
        """
        Initialize the ConjunctionFinder with threshold parameters.
        
        Args:
            mlt_threshold (float): Maximum allowed difference in MLT (hours)
            mlat_threshold (float): Maximum allowed difference in MLAT (degrees)
            min_duration_seconds (int): Minimum duration for a conjunction event (seconds)
        """
        self.mlt_threshold = mlt_threshold
        self.mlat_threshold = mlat_threshold
        self.min_duration_seconds = min_duration_seconds
        
    def load_data(self, dmsp_file, elfin_file):
        """
        Load DMSP and ELFIN satellite data from CSV files.
        
        Expected columns:
        - timestamp: Time in ISO format or unix timestamp
        - mlt: Magnetic Local Time in hours (0-24)
        - mlat: Magnetic Latitude in degrees (-90 to 90)
        
        Args:
            dmsp_file (str): Path to DMSP data file
            elfin_file (str): Path to ELFIN data file
            
        Returns:
            tuple: (dmsp_df, elfin_df) pandas DataFrames
        """
        print(f"Loading DMSP data from {dmsp_file}")
        dmsp_df = pd.read_csv(dmsp_file)
        
        print(f"Loading ELFIN data from {elfin_file}")
        elfin_df = pd.read_csv(elfin_file)
        
        # Ensure timestamp column is datetime
        for df in [dmsp_df, elfin_df]:
            if 'timestamp' in df.columns:
                if pd.api.types.is_numeric_dtype(df['timestamp']):
                    # Convert unix timestamp to datetime
                    df['timestamp'] = pd.to_datetime(df['timestamp'], unit='s')
                else:
                    # Convert string timestamp to datetime
                    df['timestamp'] = pd.to_datetime(df['timestamp'])
            else:
                raise ValueError("Data must contain a 'timestamp' column")
            
            # Check for required columns
            for col in ['mlt', 'mlat']:
                if col not in df.columns:
                    raise ValueError(f"Data must contain a '{col}' column")
        
        return dmsp_df, elfin_df
    
    def interpolate_data(self, dmsp_df, elfin_df, freq='1s'):
        """
        Interpolate data to a common time grid for comparison.
        
        Args:
            dmsp_df (DataFrame): DMSP satellite data
            elfin_df (DataFrame): ELFIN satellite data
            freq (str): Frequency for interpolation (e.g., '1s' for 1 second)
            
        Returns:
            tuple: (dmsp_interp, elfin_interp) interpolated DataFrames
        """
        # Find overlapping time range
        start_time = max(dmsp_df['timestamp'].min(), elfin_df['timestamp'].min())
        end_time = min(dmsp_df['timestamp'].max(), elfin_df['timestamp'].max())
        
        if start_time >= end_time:
            raise ValueError("No overlapping time period between DMSP and ELFIN data")
        
        print(f"Interpolating data from {start_time} to {end_time} with frequency {freq}")
        
        # Create common time index
        common_time = pd.date_range(start=start_time, end=end_time, freq=freq)
        
        # Set timestamp as index for interpolation
        dmsp_df = dmsp_df.set_index('timestamp')
        elfin_df = elfin_df.set_index('timestamp')
        
        # Interpolate to common time grid
        dmsp_interp = dmsp_df.reindex(dmsp_df.index.union(common_time)).interpolate(method='time').reindex(common_time)
        elfin_interp = elfin_df.reindex(elfin_df.index.union(common_time)).interpolate(method='time').reindex(common_time)
        
        # Reset index to have timestamp as column
        dmsp_interp = dmsp_interp.reset_index().rename(columns={'index': 'timestamp'})
        elfin_interp = elfin_interp.reset_index().rename(columns={'index': 'timestamp'})
        
        return dmsp_interp, elfin_interp
    
    def handle_mlt_wraparound(self, mlt_diff):
        """
        Handle the wraparound case for MLT differences (e.g., 23.5 - 0.5 = 23 should be 1).
        
        Args:
            mlt_diff (float or array): MLT difference(s)
            
        Returns:
            float or array: Corrected MLT difference(s)
        """
        # MLT is in range [0, 24), so the maximum difference should be 12
        return np.minimum(np.abs(mlt_diff), 24 - np.abs(mlt_diff))
    
    def find_conjunctions(self, dmsp_df, elfin_df):
        """
        Find conjunction events between DMSP and ELFIN satellites.
        
        Args:
            dmsp_df (DataFrame): DMSP satellite data
            elfin_df (DataFrame): ELFIN satellite data
            
        Returns:
            DataFrame: Conjunction events with start time, end time, and duration
        """
        # Interpolate data to common time grid
        dmsp_interp, elfin_interp = self.interpolate_data(dmsp_df, elfin_df)
        
        # Calculate differences in MLT and MLAT
        mlt_diff = self.handle_mlt_wraparound(dmsp_interp['mlt'].values - elfin_interp['mlt'].values)
        mlat_diff = np.abs(dmsp_interp['mlat'].values - elfin_interp['mlat'].values)
        
        # Create a mask for conjunction conditions
        conjunction_mask = (mlt_diff < self.mlt_threshold) & (mlat_diff < self.mlat_threshold)
        
        # Find start and end indices of conjunction events
        state_changes = np.diff(np.concatenate(([0], conjunction_mask.astype(int), [0])))
        starts = np.where(state_changes == 1)[0]
        ends = np.where(state_changes == -1)[0] - 1
        
        # Calculate duration of each event
        durations = []
        valid_events = []
        
        for i, (start_idx, end_idx) in enumerate(zip(starts, ends)):
            start_time = dmsp_interp['timestamp'].iloc[start_idx]
            end_time = dmsp_interp['timestamp'].iloc[end_idx]
            duration = (end_time - start_time).total_seconds()
            
            if duration >= self.min_duration_seconds:
                valid_events.append({
                    'start_time': start_time,
                    'end_time': end_time,
                    'duration_seconds': duration,
                    'start_idx': start_idx,
                    'end_idx': end_idx
                })
                durations.append(duration)
        
        # Create DataFrame with conjunction events
        if valid_events:
            events_df = pd.DataFrame(valid_events)
            print(f"Found {len(events_df)} conjunction events")
            return events_df, dmsp_interp, elfin_interp
        else:
            print("No conjunction events found")
            return pd.DataFrame(), dmsp_interp, elfin_interp
    
    def plot_conjunction(self, event, dmsp_interp, elfin_interp, output_dir='plots'):
        """
        Plot a conjunction event.
        
        Args:
            event (Series): Row from conjunction events DataFrame
            dmsp_interp (DataFrame): Interpolated DMSP data
            elfin_interp (DataFrame): Interpolated ELFIN data
            output_dir (str): Directory to save plots
        """
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Extract event data
        start_idx = event['start_idx']
        end_idx = event['end_idx']
        
        # Add some margin for context
        margin = int(min(300, (end_idx - start_idx) * 0.5))  # 5 minutes or 50% of event duration
        plot_start = max(0, start_idx - margin)
        plot_end = min(len(dmsp_interp), end_idx + margin)
        
        # Extract time slice for plotting
        time_slice = slice(plot_start, plot_end)
        time = dmsp_interp['timestamp'].iloc[time_slice]
        
        # Create figure with subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
        
        # Plot MLT
        ax1.plot(time, dmsp_interp['mlt'].iloc[time_slice], 'b-', label='DMSP')
        ax1.plot(time, elfin_interp['mlt'].iloc[time_slice], 'r-', label='ELFIN')
        ax1.axvspan(event['start_time'], event['end_time'], alpha=0.2, color='green')
        ax1.set_ylabel('MLT (hours)')
        ax1.legend()
        ax1.grid(True)
        
        # Plot MLAT
        ax2.plot(time, dmsp_interp['mlat'].iloc[time_slice], 'b-', label='DMSP')
        ax2.plot(time, elfin_interp['mlat'].iloc[time_slice], 'r-', label='ELFIN')
        ax2.axvspan(event['start_time'], event['end_time'], alpha=0.2, color='green')
        ax2.set_ylabel('MLAT (degrees)')
        ax2.set_xlabel('Time (UTC)')
        ax2.legend()
        ax2.grid(True)
        
        # Format x-axis
        fig.autofmt_xdate()
        
        # Add title
        start_str = event['start_time'].strftime('%Y-%m-%d %H:%M:%S')
        duration_min = event['duration_seconds'] / 60
        fig.suptitle(f'DMSP-ELFIN Conjunction Event\nStart: {start_str}, Duration: {duration_min:.1f} minutes')
        
        # Save figure
        filename = f"conjunction_{event['start_time'].strftime('%Y%m%d_%H%M%S')}.png"
        filepath = os.path.join(output_dir, filename)
        plt.tight_layout()
        plt.savefig(filepath)
        plt.close()
        
        print(f"Saved plot to {filepath}")
        
    def analyze_conjunctions(self, events_df, dmsp_interp, elfin_interp):
        """
        Analyze conjunction events and generate summary statistics.
        
        Args:
            events_df (DataFrame): Conjunction events
            dmsp_interp (DataFrame): Interpolated DMSP data
            elfin_interp (DataFrame): Interpolated ELFIN data
            
        Returns:
            dict: Summary statistics
        """
        if events_df.empty:
            return {"count": 0}
        
        # Calculate statistics
        durations = events_df['duration_seconds'].values
        
        stats = {
            "count": len(events_df),
            "total_duration_seconds": np.sum(durations),
            "min_duration_seconds": np.min(durations),
            "max_duration_seconds": np.max(durations),
            "mean_duration_seconds": np.mean(durations),
            "median_duration_seconds": np.median(durations),
        }
        
        # Plot each conjunction event
        for i, event in events_df.iterrows():
            self.plot_conjunction(event, dmsp_interp, elfin_interp)
        
        return stats
    
    def save_results(self, events_df, stats, output_file='conjunction_events.csv', stats_file='conjunction_stats.json'):
        """
        Save conjunction events and statistics to files.
        
        Args:
            events_df (DataFrame): Conjunction events
            stats (dict): Summary statistics
            output_file (str): Path to save events
            stats_file (str): Path to save statistics
        """
        if not events_df.empty:
            events_df.to_csv(output_file, index=False)
            print(f"Saved {len(events_df)} conjunction events to {output_file}")
        
        # Save statistics as JSON
        import json
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        print(f"Saved statistics to {stats_file}")
        
        # Print summary
        print("\nConjunction Analysis Summary:")
        print(f"Found {stats.get('count', 0)} conjunction events")
        if stats.get('count', 0) > 0:
            print(f"Total duration: {stats['total_duration_seconds']/60:.1f} minutes")
            print(f"Mean duration: {stats['mean_duration_seconds']:.1f} seconds")
            print(f"Median duration: {stats['median_duration_seconds']:.1f} seconds")
            print(f"Min duration: {stats['min_duration_seconds']:.1f} seconds")
            print(f"Max duration: {stats['max_duration_seconds']:.1f} seconds")


def generate_sample_data(output_dir='.', duration_hours=24, sample_rate_seconds=10):
    """
    Generate sample DMSP and ELFIN data for testing.
    
    Args:
        output_dir (str): Directory to save sample data
        duration_hours (int): Duration of sample data in hours
        sample_rate_seconds (int): Time between samples in seconds
        
    Returns:
        tuple: (dmsp_file, elfin_file) paths to generated files
    """
    print(f"Generating sample data for {duration_hours} hours at {sample_rate_seconds}s sample rate")
    
    # Create time range
    start_time = datetime(2023, 1, 1, 0, 0, 0)
    end_time = start_time + timedelta(hours=duration_hours)
    timestamps = pd.date_range(start=start_time, end=end_time, freq=f'{sample_rate_seconds}s')
    
    # Number of samples
    n_samples = len(timestamps)
    
    # Generate DMSP data
    dmsp_data = {
        'timestamp': timestamps,
        'mlt': (12 + 6 * np.sin(np.linspace(0, 2*np.pi, n_samples))) % 24,  # MLT varies from 6 to 18
        'mlat': 60 + 20 * np.sin(np.linspace(0, 4*np.pi, n_samples))  # MLAT varies from 40 to 80
    }
    dmsp_df = pd.DataFrame(dmsp_data)
    
    # Generate ELFIN data with some overlap and some differences
    elfin_data = {
        'timestamp': timestamps,
        'mlt': (13 + 6 * np.sin(np.linspace(0.5, 2.5*np.pi, n_samples))) % 24,  # Slightly offset MLT
        'mlat': 58 + 20 * np.sin(np.linspace(0.2, 4.2*np.pi, n_samples))  # Slightly offset MLAT
    }
    elfin_df = pd.DataFrame(elfin_data)
    
    # Save to CSV
    dmsp_file = os.path.join(output_dir, 'dmsp_sample.csv')
    elfin_file = os.path.join(output_dir, 'elfin_sample.csv')
    
    dmsp_df.to_csv(dmsp_file, index=False)
    elfin_df.to_csv(elfin_file, index=False)
    
    print(f"Saved sample DMSP data to {dmsp_file}")
    print(f"Saved sample ELFIN data to {elfin_file}")
    
    return dmsp_file, elfin_file


def main():
    """Main function to run the conjunction finder."""
    parser = argparse.ArgumentParser(description='Find conjunction events between DMSP and ELFIN satellites')
    parser.add_argument('--dmsp', type=str, help='Path to DMSP data file')
    parser.add_argument('--elfin', type=str, help='Path to ELFIN data file')
    parser.add_argument('--mlt-threshold', type=float, default=0.5, help='Maximum allowed difference in MLT (hours)')
    parser.add_argument('--mlat-threshold', type=float, default=2.0, help='Maximum allowed difference in MLAT (degrees)')
    parser.add_argument('--min-duration', type=int, default=60, help='Minimum duration for a conjunction event (seconds)')
    parser.add_argument('--output', type=str, default='conjunction_events.csv', help='Output file for conjunction events')
    parser.add_argument('--generate-sample', action='store_true', help='Generate sample data for testing')
    parser.add_argument('--sample-duration', type=int, default=24, help='Duration of sample data in hours')
    
    args = parser.parse_args()
    
    # Generate sample data if requested
    if args.generate_sample:
        dmsp_file, elfin_file = generate_sample_data(duration_hours=args.sample_duration)
    else:
        dmsp_file, elfin_file = args.dmsp, args.elfin
    
    # Check if files are provided
    if not dmsp_file or not elfin_file:
        parser.error("Please provide paths to DMSP and ELFIN data files or use --generate-sample")
    
    # Initialize conjunction finder
    finder = ConjunctionFinder(
        mlt_threshold=args.mlt_threshold,
        mlat_threshold=args.mlat_threshold,
        min_duration_seconds=args.min_duration
    )
    
    # Load data
    dmsp_df, elfin_df = finder.load_data(dmsp_file, elfin_file)
    
    # Find conjunctions
    events_df, dmsp_interp, elfin_interp = finder.find_conjunctions(dmsp_df, elfin_df)
    
    # Analyze conjunctions
    stats = finder.analyze_conjunctions(events_df, dmsp_interp, elfin_interp)
    
    # Save results
    finder.save_results(events_df, stats, output_file=args.output)


if __name__ == "__main__":
    main()
